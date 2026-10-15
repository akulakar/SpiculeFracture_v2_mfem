#include "spicule/Setup.hpp"

#include <stdexcept>

using namespace mfem;

namespace spicule
{
    std::unique_ptr<mfem::ParMesh> BuildParMesh(const MeshParams &mesh_params, MPI_Comm comm)
    {
        mfem::Mesh serial(mesh_params.file.c_str(), 1, 1);
        serial.SetCurvature(mesh_params.order, false);
        // generate par mesh
        auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial);
        // refine par mesh
        for (int l = 0; l < mesh_params.refLevels; l++)
            pmesh->UniformRefinement();
        return pmesh;
    }

    int BasisFromName(const std::string &name)
    {
        if (name == "GaussLobatto")
            return BasisType::GaussLobatto;
        else if (name == "Positive")
            return BasisType::Positive;
        else if (name == "ClosedUniform")
            return BasisType::ClosedUniform;

        throw std::runtime_error("Unknown fe.basis: " + name);
    }

    int OrderingType(const std::string &ordering)
    {
        if (ordering == "byNODES")
            return Ordering::byNODES;
        else if (ordering == "byVDIM")
            return Ordering::byVDIM;
        else
            return Ordering::byNODES; // default is byNODES
    }

    void Fields::Update()
    {
        u = utemp;
        d = dtemp;
    }

    void Fields::Revert()
    {
        utemp = u;
        dtemp = d;
    }

    BoundaryConditions::BoundaryConditions(const std::vector<BCEntry> &bcs, ParFiniteElementSpace &disp)
    {
        ParMesh &mesh = *disp.GetParMesh();
        const int num_attr = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
        const int dim = mesh.Dimension();

        for (const BCEntry &bc : bcs)
        {
            Array<int> marker(num_attr);
            marker = 0;

            // Check if boundaries are valid
            for (int a : bc.attrs)
            {
                if (a < 0 || a >= num_attr)
                    throw std::runtime_error("bcs: boundary attribute " + std::to_string(a) + " not present in mesh");
                marker[a] = 1;
            }

            for (int comp : bc.components)
            {
                if (comp < 0 || comp >= dim)
                    throw std::runtime_error("bcs: component out of range");

                // Check if dirichlet dofs from each attribute
                Entry entry;
                disp.GetEssentialTrueDofs(marker, entry.ess_true_dofs, comp);
                entry.isLoad = bc.isNonZero;
                entry.value = bc.value;
                ess_tdofs_.Append(entry.ess_true_dofs);
                entries_.push_back(std::move(entry));
            }
        }
        ess_tdofs_.Sort();
        ess_tdofs_.Unique();
    }

    void BoundaryConditions::Apply(double load, ParGridFunction &u) const
    {
        u.GetTrueDofs(tv_);

        for (const Entry &entry : entries_)
        {
            const double disp_val = entry.isLoad ? load : entry.value;

            for (int i = 0; i < entry.ess_true_dofs.Size(); i++)
                tv_(entry.ess_true_dofs[i]) = disp_val;
        }

        u.SetFromTrueDofs(tv_);
    }

} // namespace spicule
