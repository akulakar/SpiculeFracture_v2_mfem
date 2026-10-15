#pragma once
// Mesh, finite element spaces, fields and boundary conditions.

#include <memory>
#include <mfem.hpp>
#include <vector>

#include "spicule/Config.hpp"

namespace spicule
{

    std::unique_ptr<mfem::ParMesh> BuildParMesh(const MeshParams &mesh_params, MPI_Comm comm);

    int BasisFromName(const std::string &name);
    int OrderingType(const std::string &ordering);

    struct FESpaces
    {
        std::unique_ptr<mfem::FiniteElementCollection> fec;
        std::unique_ptr<mfem::ParFiniteElementSpace> dispfespace;   // vector H1, byNODES
        std::unique_ptr<mfem::ParFiniteElementSpace> damagefespace; // scalar H1

        FESpaces(Config &ConfigFile, mfem::ParMesh *pmesh)
        {
            // H1 Lagrange FE collection.
            fec = std::make_unique<mfem::H1_FECollection>(ConfigFile.fe.order, pmesh->Dimension(), spicule::BasisFromName(ConfigFile.fe.basis));

            // displacement and damage finite element spaces, disp is vector valued space, damage is scalar. Both share same element order.
            dispfespace = std::make_unique<mfem::ParFiniteElementSpace>(pmesh, fec.get(), pmesh->Dimension(), spicule::OrderingType(ConfigFile.fe.ordering));
            damagefespace = std::make_unique<mfem::ParFiniteElementSpace>(pmesh, fec.get(), 1, spicule::OrderingType(ConfigFile.fe.ordering));

            pmesh->SetNodalFESpace(dispfespace.get());
            pmesh->EnsureNodes();
        }
    };

    // u, d: updated state. utemp, dtemp: temporary state read by the integrators.
    struct Fields
    {
        mfem::ParGridFunction u, d;
        mfem::ParGridFunction utemp, dtemp;

        Fields(FESpaces &fe_spaces) : u(fe_spaces.dispfespace.get()), d(fe_spaces.damagefespace.get()),
                                      utemp(fe_spaces.dispfespace.get()), dtemp(fe_spaces.damagefespace.get())
        {
            u = 0.0;
            utemp = 0.0;

            d = 0.0;
            dtemp = 0.0;
        };

        void Update(); // temp -> updated
        void Revert(); // accepted -> temp
    };

    // Dirichlet conditions from the bcs input array. Entries apply in input order,
    // so a later entry wins on shared dofs. Dof lists are computed once.
    class BoundaryConditions
    {

    private:
        struct Entry
        {
            mfem::Array<int> ess_true_dofs;
            bool isLoad;
            double value;
        };
        std::vector<Entry> entries_;
        mfem::Array<int> ess_tdofs_;
        mutable mfem::Vector tv_;

    public:
        BoundaryConditions(const std::vector<BCEntry> &bcs, mfem::ParFiniteElementSpace &disp);

        // Return essential true dofs
        const mfem::Array<int> &EssentialTrueDofs() const { return ess_tdofs_; }

        // Overwrites the dirichlet entries of u with the BC values (load-driven ones use load).
        void Apply(double load, mfem::ParGridFunction &u) const;
    };

} // namespace spicule
