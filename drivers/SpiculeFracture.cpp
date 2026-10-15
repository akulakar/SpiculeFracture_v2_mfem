// Regularized variational fracture driver. All choices come from the input JSON:
//   ./SpiculeFracture ../input/SpiculeFracture.json [../input/cases/override.json]
// The pipeline: Config -> Mesh -> FESpaces -> Fields/BCs -> Forms -> Solver -> Load loop -> Output.

#include <iostream>
#include <mfem.hpp>

#include "spicule/Config.hpp"
#include "spicule/Forms.hpp"
#include "spicule/Output.hpp"
#include "spicule/PhaseFieldSolver.hpp"
#include "spicule/Setup.hpp"

using namespace mfem;

int main(int argc, char *argv[])
{
    MPI_Comm comm = MPI_COMM_WORLD;
    Mpi::Init(argc, argv);
    Hypre::Init();
    Device device("cpu");

    if (argc < 2)
    {
        if (Mpi::Root())
            std::cerr << "Usage: " << argv[0] << " input.json [override.json]" << std::endl;
        return 1;
    }

    try
    {
        StopWatch timer;
        timer.Start();

        //----------------------------------------------------------------------------------------------------------------
        // Load input file
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Load input file and configure." << std::endl;

        // Load config from input file
        spicule::Config config = spicule::LoadConfig(argv[1], argc > 2 ? argv[2] : "");

        if (Mpi::Root())
            std::cout << "Run name: " << config.testName << std::endl;

        //----------------------------------------------------------------------------------------------------------------
        // Generate mesh
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Generate mesh." << std::endl;

        // Read serial mesh and convert to parallel mesh
        auto pmesh = BuildParMesh(config.mesh, comm);

        //----------------------------------------------------------------------------------------------------------------
        // Initialize finite element spaces
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize finite element spaces." << std::endl;

        // H1 Lagrange FE collection.
        // Lagrange finite elements for displacement and damage.
        spicule::FESpaces fe_spaces(config, pmesh.get());

        //----------------------------------------------------------------------------------------------------------------
        // Initialize grid functions for displacement and damage fields
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize grid functions." << std::endl;

        // Grid functions for displacement and damage fields
        spicule::Fields fields(fe_spaces);

        //----------------------------------------------------------------------------------------------------------------
        // Initialize bilinear and linear forms
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize bilinear and linear forms." << std::endl;

        // Bilinear and linear forms for the elasticity and damage equations
        spicule::Forms forms(config, fe_spaces, fields);
        forms.AddIntegrators(config);

        //----------------------------------------------------------------------------------------------------------------
        // Initialize boundary conditions class
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize boundary conditions." << std::endl;

        // Boundary conditions read attributes
        spicule::BoundaryConditions bcs(config.bcs, *(fe_spaces.dispfespace));

        //----------------------------------------------------------------------------------------------------------------
        // Initialize main phase field solver class
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize phase field fracture solver." << std::endl;

        spicule::PhaseFieldSolver phase_field_solver(config, fe_spaces, fields, forms, bcs, comm);

        //----------------------------------------------------------------------------------------------------------------
        // Initialize output writer
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize output." << std::endl;

        spicule::OutputWriter output(config, *pmesh, fields, comm);

        const auto num_disp_unknowns = fe_spaces.dispfespace->GlobalTrueVSize();
        const auto num_damage_unknowns = fe_spaces.damagefespace->GlobalTrueVSize();

        if (Mpi::Root())
            std::cout << num_disp_unknowns << " displacement dofs, "
                      << num_damage_unknowns << " damage dofs, Ey " << config.material.Ey << ", nu " << config.material.nu << std::endl;

        //----------------------------------------------------------------------------------------------------------------
        // Start staggered iterations
        //----------------------------------------------------------------------------------------------------------------

        // cycle number and time
        output.Save(0, 0.0);

        // Staggered solve loop
        const bool finished_successfully = phase_field_solver.RunStaggeredSolve([&](int step, double load)
                                                                                {   if(step % config.output.saveFrequency == 0)
                                                                                    output.Save(step, load); });

        //----------------------------------------------------------------------------------------------------------------
        // Final summary
        //----------------------------------------------------------------------------------------------------------------

        timer.Stop();

        if (!finished_successfully)
        {
            if (Mpi::Root())
                std::cerr << "Load stepping failed: increment collapsed without convergence" << std::endl;
            return 2;
        }

        if (Mpi::Root())
            std::cout << "Simulation completed in " << timer.RealTime() << " s" << std::endl;
    }

    catch (const std::exception &excep)
    {
        if (Mpi::Root())
            std::cerr << "Error: " << excep.what() << std::endl;
        return 1;
    }

    return 0;
}
