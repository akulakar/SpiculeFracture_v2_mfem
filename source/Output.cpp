#include "spicule/Output.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace spicule
{
    OutputWriter::OutputWriter(const Config &config, mfem::ParMesh &mesh, Fields &fields, MPI_Comm comm)
        : mesh_(mesh), fields_(fields)
    {
        int myid, nranks;
        MPI_Comm_rank(comm, &myid);
        MPI_Comm_size(comm, &nranks);

        namespace fs = std::filesystem;
        const fs::path dir = fs::path(config.output.path) / config.testName;
        if (myid == 0)
            fs::create_directories(dir);
        MPI_Barrier(comm);

        const long long nElem = mesh.GetGlobalNE();
        const long long nDisp = fields.u.ParFESpace()->GlobalTrueVSize();
        const long long nDmg = fields.d.ParFESpace()->GlobalTrueVSize();
        if (myid == 0)
        {
            std::ofstream(dir / "input.json") << config.inputJson.dump(2) << "\n";
            std::ofstream det(dir / "SimulationDetails.txt");
            det << "testName\t" << config.testName << "\n"
                << "mpi_ranks\t" << nranks << "\n"
                << "Ey\t" << config.material.Ey << "\n"
                << "nu\t" << config.material.nu << "\n"
                << "mesh\t" << config.mesh.file << "\n"
                << "mesh_order\t" << config.mesh.order << "\n"
                << "ref_levels\t" << config.mesh.refLevels << "\n"
                << "fe_order\t" << config.fe.order << "\n"
                << "global_elements\t" << nElem << "\n"
                << "disp_dofs\t" << nDisp << "\n"
                << "damage_dofs\t" << nDmg << "\n";
        }

        const std::string &fmt = config.output.format;
        if (fmt == "adios2")
        {
            aos_ = std::make_unique<mfem::adios2stream>(
                (dir / "FractureMFEMoutput.bp").string(), mfem::adios2stream::openmode::out, comm, "BP5");
            aos_->SetRefinementLevel(0);
            aos_->SetParameter("FlushStepsCount", "1");
        }
        else if (fmt == "paraview")
        {
            pv_ = std::make_unique<mfem::ParaViewDataCollection>("FractureMFEMoutput", &mesh);
            pv_->SetPrefixPath(dir.string());
            pv_->SetHighOrderOutput(true);
            pv_->RegisterField("disp", &fields.u);
            pv_->RegisterField("damage", &fields.d);
        }
        else
            throw std::runtime_error("Unknown output.format: " + fmt);
    }

    void OutputWriter::Save(int cycle, double time)
    {
        if (aos_)
        {
            aos_->BeginStep();
            mesh_.Print(*aos_);
            aos_->SetTime(time);
            aos_->SetCycle(cycle);
            fields_.u.Save(*aos_, "disp", mfem::adios2stream::data_type::point_data);
            fields_.d.Save(*aos_, "damage", mfem::adios2stream::data_type::point_data);
            aos_->EndStep();
        }
        else
        {
            pv_->SetCycle(cycle);
            pv_->SetTime(time);
            pv_->Save();
        }
    }

} // namespace spicule
