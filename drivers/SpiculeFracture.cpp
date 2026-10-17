// Regularized variational fracture driver. All choices come from the input JSON:
//   ./SpiculeFracture ../input/inputParameters/SpiculeFracture.json [override.json]
// Pipeline: input JSON -> mesh -> solver setup (FE spaces, fields, forms, BCs, linear solvers) -> load loop -> output.

#include <fstream>
#include <iostream>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "PhaseFieldFracture.hpp"

using namespace mfem;
using json = nlohmann::json;

// Read a JSON file; throws if it cannot be opened.
static json ReadJson(const std::string &path)
{
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Could not open the input file: " + path);
    json j;
    in >> j;
    return j;
}
// -------------------------------------------------------------------------------------------------------
// Small helpers
// -------------------------------------------------------------------------------------------------------
static int BasisFromName(const std::string &name)
{
    if (name == "GaussLobatto")
        return mfem::BasisType::GaussLobatto;
    if (name == "Positive")
        return mfem::BasisType::Positive;
    if (name == "ClosedUniform")
        return mfem::BasisType::ClosedUniform;
    throw std::runtime_error("Unknown Finite_Element.basis: " + name);
}

static int OrderingFromName(const std::string &name)
{
    if (name == "byNODES")
        return mfem::Ordering::byNODES;
    if (name == "byVDIM")
        return mfem::Ordering::byVDIM;
    throw std::runtime_error("Unknown Finite_Element.ordering: " + name + " (use byNODES or byVDIM)");
}

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
        // Load input file (an optional second file is merged on top of it)
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Load input file." << std::endl;

        json input = ReadJson(argv[1]);
        if (argc > 2)
            input.merge_patch(ReadJson(argv[2]));

        if (Mpi::Root())
            std::cout << "Run name: " << input.at("Test_Name").get<std::string>() << std::endl;

        //----------------------------------------------------------------------------------------------------------------
        // Generate mesh
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Generate mesh." << std::endl;

        const json &mesh_input = input.at("Mesh");
        const int mesh_order = mesh_input.value("order", 1);
        const int refinement_levels = mesh_input.value("refinement_levels", 0);

        // Read the serial mesh, partition it into parallel mesh, then free the serial mesh.
        std::unique_ptr<ParMesh> pmesh;
        {
            Mesh serial_mesh(mesh_input.at("file").get<std::string>().c_str(), 1, 1);
            serial_mesh.SetCurvature(mesh_order, false);
            pmesh = std::make_unique<ParMesh>(comm, serial_mesh);
        }
        for (int l = 0; l < refinement_levels; l++)
            pmesh->UniformRefinement();

        const int dim = pmesh->Dimension();

        //----------------------------------------------------------------------------------------------------------------
        // Initialize finite element collection and finite element spaces for displacement and damage fields.
        //----------------------------------------------------------------------------------------------------------------

        // H1 Lagrange collection and the vector displacement / scalar damage spaces (same element order).
        const json &fe = input.at("Finite_Element");
        const int order = fe.value("order", input.at("Mesh").value("order", 1));
        const int ordering = OrderingFromName(fe.value("ordering", "byNODES"));

        H1_FECollection fec(order, dim, BasisFromName(fe.value("basis", "GaussLobatto")));
        ParFiniteElementSpace dispfes(pmesh.get(), &fec, dim, ordering);
        ParFiniteElementSpace dmgfes(pmesh.get(), &fec, 1, ordering);

        pmesh->SetNodalFESpace(&dispfes);
        pmesh->EnsureNodes();

        //----------------------------------------------------------------------------------------------------------------
        // Initialize grid functions for displacement and damage fields.
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize fields." << std::endl;

        ParGridFunction u(&dispfes), utemp(&dispfes);
        ParGridFunction d(&dmgfes), dtemp(&dmgfes);

        u = 0.0;
        utemp = 0.0;
        d = 0.0;
        dtemp = 0.0;

        //----------------------------------------------------------------------------------------------------------------
        // Initialize bilinear and linear forms for elasticity and damage equations.
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize bilinear and linear forms." << std::endl;

        //----------------------------------------------------------------------------------------------------------------
        // Initialize coefficients

        const json &material = input.at("Material");
        const double Ey = material.at("youngs_modulus");
        const double nu = material.at("poisson_ratio");
        const double Gc = material.at("fracture_toughness");
        const double eps = material.at("length_scale");
        const double Keps = material.at("K_epsilon");

        const json &eigenstrain = input.at("Eigenstrain");
        const double beta = eigenstrain.value("swelling_coefficient", 0.0);
        const double deltaC = eigenstrain.value("Delta_concentration", 0.0);

        ConstantCoefficient Ey_coeff(Ey), nu_coeff(nu), volEig_coeff(beta * deltaC);
        ConstantCoefficient gOverEps_coeff(Gc / eps), gTimesEps_coeff(Gc * eps);

        //----------------------------------------------------------------------------------------------------------------
        // Bilinear form for elasticity
        ParBilinearForm kdisp(&dispfes);
        kdisp.AddDomainIntegrator(new mfemplus::DamagedIsotropicElasticityIntegrator(Ey_coeff, nu_coeff, Keps, dtemp, &dmgfes));

        //----------------------------------------------------------------------------------------------------------------
        // Bilinear forms for the damage equation. Split in two because only kdmg2 is reassembled every iteration.
        ParBilinearForm kdmg1(&dmgfes);
        kdmg1.AddDomainIntegrator(new mfem::MassIntegrator(gOverEps_coeff));
        kdmg1.AddDomainIntegrator(new mfem::DiffusionIntegrator(gTimesEps_coeff));

        ParBilinearForm kdmg2(&dmgfes);
        kdmg2.AddDomainIntegrator(new mfemplus::IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(Ey_coeff, nu_coeff, volEig_coeff, utemp, &dispfes));

        //----------------------------------------------------------------------------------------------------------------
        // Linear form for the elasticity equation
        ParLinearForm fdisp(&dispfes);
        fdisp.AddDomainIntegrator(new mfemplus::DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator(Ey_coeff, nu_coeff, volEig_coeff, Keps, dtemp, &dmgfes));

        //----------------------------------------------------------------------------------------------------------------
        // Linear form for the damage equation
        ParLinearForm fdmg(&dmgfes);
        fdmg.AddDomainIntegrator(new mfemplus::IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(Ey_coeff, nu_coeff, volEig_coeff, utemp, &dispfes));

        //----------------------------------------------------------------------------------------------------------------
        // Initialize boundary conditions object.
        //----------------------------------------------------------------------------------------------------------------

        if (Mpi::Root())
            std::cout << "Initialize boundary conditions." << std::endl;

        const json &bc = input.at("Boundary_Conditions");
        spicule::BoundaryConditions BCs(bc, *pmesh, dispfes);

        //----------------------------------------------------------------------------------------------------------------
        // Set up the phase field solver. Each step reads its own section of the input.
        //----------------------------------------------------------------------------------------------------------------
        spicule::PhaseFieldSolver solver(input, *pmesh, comm);

        if (Mpi::Root())
            std::cout << "Initialize phase field solver." << std::endl;

        solver.ReadFESpaces(&dispfes, &dmgfes);
        solver.ReadFields(&u, &utemp, &d, &dtemp);
        solver.ReadForms(&kdisp, &kdmg1, &kdmg2, &fdisp, &fdmg);
        solver.ReadBCs(&BCs);
        solver.InitializeSolvers(); // preconditioners and Krylov solvers

        //----------------------------------------------------------------------------------------------------------------
        // Output
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Initialize output." << std::endl;

        spicule::OutputWriter output(input, *pmesh, solver, comm);

        const auto num_disp_unknowns = solver.DispFESpace().GlobalTrueVSize();
        const auto num_damage_unknowns = solver.DamageFESpace().GlobalTrueVSize();

        if (Mpi::Root())
            std::cout << num_disp_unknowns << " displacement dofs, " << num_damage_unknowns << " damage dofs" << std::endl;

        const int save_frequency = std::max(1, input.value("Output", json::object()).value("save_frequency", 1));

        //----------------------------------------------------------------------------------------------------------------
        // Staggered solve over the load steps
        //----------------------------------------------------------------------------------------------------------------
        if (Mpi::Root())
            std::cout << "Starting staggered solves." << std::endl;

        output.Save(0, 0.0);

        const bool finished_successfully = solver.RunStaggeredSolve([&](int step, double load)
                                                                    {if (step % save_frequency == 0)
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
