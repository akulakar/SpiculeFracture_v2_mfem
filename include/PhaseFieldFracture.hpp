#pragma once
// Header-only classes for the regularized variational (AT2) fracture driver.
//
//   BoundaryConditions  Dirichlet dof lists from the input, applied to the working displacement.
//   PhaseFieldSolver    owns the FE spaces, fields, forms and linear solvers; runs the staggered load stepping.
//   OutputWriter        ADIOS2 BP5 or ParaView output plus the resolved input and run details.
//
// The solver is set up with small Initialize*() methods that the driver calls in this order:
//   InitializeFESpaces -> InitializeFields -> InitializeForms -> InitializeBCs -> InitializeSolvers
// Each method reads only its own section of the input JSON.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <mfem.hpp>
#include <mfemplus/SolidMechanics/mfemplus_solids.hpp>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <vector>

namespace spicule
{
    using json = nlohmann::json;

    //------------------------------------------------------------------------------------------------------------------
    // Dirichlet boundary conditions on displacement.
    //------------------------------------------------------------------------------------------------------------------
    // Built from Boundary_Conditions.dirichlet_attributes. Entries apply in input order, so a later entry wins on
    // shared dofs. Boundary attributes are 0-based indices into the mesh's boundary attributes (0 = mesh attribute 1).
    // A value of "nonzero_displacement" means the entry follows the load controller.
    class BoundaryConditions
    {
    private:
        struct Entry
        {
            mfem::Array<int> ess_true_dofs;
            bool isLoad = false;
            double value = 0.0;
        };
        std::vector<Entry> entries_;
        mfem::Array<int> ess_tdofs_;
        mutable mfem::Vector tv_;

        const json &dirichlet_attributes;

    public:
        double finalDisp_ = 0.0, minInc_ = 0.0, maxInc_ = 0.0;

    public:
        BoundaryConditions(const json &boundary_conditions, mfem::ParMesh &pmesh, mfem::ParFiniteElementSpace &disp) : dirichlet_attributes(boundary_conditions.at("dirichlet_attributes"))
        {

            const int num_attr = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
            const int dim = pmesh.Dimension();

            finalDisp_ = boundary_conditions.at("nonzero_displacement").at("final_displacement");
            minInc_ = boundary_conditions.at("nonzero_displacement").at("min_increment");
            maxInc_ = boundary_conditions.at("nonzero_displacement").at("max_increment");

            // Loop through each attribute
            for (const auto &bc : dirichlet_attributes)
            {
                const auto attributes = bc.at("attributes").get<std::vector<int>>();
                const auto components = bc.at("components").get<std::vector<int>>();

                bool is_load = false;
                double value = 0.0;

                if (bc.contains("value") && bc.at("value").is_string())
                {
                    if (bc.at("value").get<std::string>() != "nonzero_displacement")
                        throw std::runtime_error("Boundary condition value must be a number or \"nonzero_displacement\"");
                    is_load = true;
                }
                else
                    value = bc.value("value", 0.0);

                mfem::Array<int> marker(num_attr);
                marker = 0;
                for (int a : attributes)
                {
                    if (a < 0 || a >= num_attr)
                        throw std::runtime_error("Boundary attribute " + std::to_string(a) + " not present in mesh");
                    marker[a] = 1;
                }

                for (int comp : components)
                {
                    if (comp < 0 || comp >= dim)
                        throw std::runtime_error("Boundary condition component out of range");

                    Entry entry;
                    disp.GetEssentialTrueDofs(marker, entry.ess_true_dofs, comp);
                    entry.isLoad = is_load;
                    entry.value = value;
                    ess_tdofs_.Append(entry.ess_true_dofs);
                    entries_.push_back(std::move(entry));
                }
            }
            ess_tdofs_.Sort();
            ess_tdofs_.Unique();
        }

        const mfem::Array<int> &EssentialTrueDofs() const { return ess_tdofs_; }

        // Overwrites the Dirichlet entries of u with the BC values (load-driven entries use `load`).
        void Apply(double load, mfem::ParGridFunction &u) const
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
    };

    //------------------------------------------------------------------------------------------------------------------
    // Staggered phase-field fracture solver.
    //------------------------------------------------------------------------------------------------------------------
    //   bilinear:  Kdisp(d)  damaged isotropic elasticity (reads dtemp)
    //              Kdmg1     Gc/eps mass + Gc*eps diffusion (constant, assembled once)
    //              Kdmg2(u)  strain-energy driven reaction term (reads utemp)
    //   linear:    fdisp(d)  volumetric eigenstrain body force
    //              fdmg(u)   strain-energy driving force
    // u, d are the accepted fields; utemp, dtemp are the working copies the integrators read.
    class PhaseFieldSolver
    {
    protected:
        const json &input_;
        mfem::ParMesh &pmesh_;
        MPI_Comm comm_;

        // // Finite element spaces and fields
        // std::unique_ptr<mfem::FiniteElementCollection> fec_;
        // std::unique_ptr<mfem::ParFiniteElementSpace> dispfes_, dmgfes_;
        // std::unique_ptr<mfem::ParGridFunction> u_, utemp_, d_, dtemp_;

        // // Coefficients, forms
        // mfem::ConstantCoefficient Ey_, nu_, volEig_, gOverEps_, gTimesEps_;
        // std::unique_ptr<mfem::ParBilinearForm> kdisp_, kdmg1_, kdmg2_;
        // std::unique_ptr<mfem::ParLinearForm> fdisp_, fdmg_;

        // Finite element spaces and fields;
        mfem::ParFiniteElementSpace *dispfes_, *dmgfes_;
        mfem::ParGridFunction *u_, *utemp_, *d_, *dtemp_;

        // Coefficients, forms
        mfem::ParBilinearForm *kdisp_, *kdmg1_, *kdmg2_;
        mfem::ParLinearForm *fdisp_, *fdmg_;

        // Boundary conditions and load range
        BoundaryConditions *bcs_;

        // Persistent solvers. Krylov solvers may be mfem (CG, GMRES) or hypre (PCG, GMRES), chosen
        // independently for displacement and damage.

        std::unique_ptr<mfem::HypreBoomerAMG> dispPrec_, dmgPrec_;
        std::unique_ptr<mfem::Solver> dispKrylov_, dmgKrylov_;
        int dispMaxIter_ = 0, dmgMaxIter_ = 0;
        double staggerTol_ = 1.0e-6;
        int staggerMaxIter_ = 1000;

        // Linear system data. Matrices are swapped each iteration.
        std::unique_ptr<mfem::HypreParMatrix> Kdisp_, Kdmg1Mat_, Kdmg2Mat_, Kdmg_;
        std::unique_ptr<mfem::HypreParVector> Fdmg_;
        mfem::Vector u_new_, F_disp_, dOld_, dNew_;

        mfem::StopWatch disp_assembly_timer_, dmg_assembly_timer_, disp_solve_timer_, dmg_solve_timer_;

        // -------------------------------------------------------------------------------------------------------
        // Small helpers
        // -------------------------------------------------------------------------------------------------------

        template <class T>
        static void Require(const std::unique_ptr<T> &ptr, const char *what, const char *call_first)
        {
            if (!ptr)
                throw std::runtime_error(std::string(what) + " not initialized; call " + call_first + " first");
        }

        // Preconditioner from the input: none | amg | amg_systems | amg_elasticity.
        static std::unique_ptr<mfem::HypreBoomerAMG> MakePrec(const std::string &kind, mfem::ParFiniteElementSpace &fes)
        {
            if (kind == "none")
                return nullptr;

            auto amg = std::make_unique<mfem::HypreBoomerAMG>();
            amg->SetPrintLevel(0);
            if (kind == "amg_systems")
                amg->SetSystemsOptions(fes.GetVDim(), fes.GetOrdering() == mfem::Ordering::byNODES);
            else if (kind == "amg_elasticity")
            {
                if (fes.GetOrdering() != mfem::Ordering::byVDIM)
                    throw std::runtime_error("amg_elasticity requires Finite_Element.ordering = byVDIM");
                amg->SetElasticityOptions(&fes, true);
            }
            else if (kind == "amg")
            {
                amg->SetCoarsening(10); // HMIS
                amg->SetAggressiveCoarsening(2);
                amg->SetRelaxType(6); // symmetric GS
                amg->SetCycleType(1); // V
                amg->SetStrongThresholdR(0.25);
                amg->SetTol(0.0);
                amg->SetMaxIter(1);
            }
            else
                throw std::runtime_error("Unknown preconditioner: " + kind);
            return amg;
        }

        // Krylov solver from a solver block: library = mfem | hypre, type = pcg | gmres.
        // The stopping criteria differ: mfem solvers use `tolerance` relative to the initial preconditioned
        // residual, hypre solvers relative to ||b||.
        std::unique_ptr<mfem::Solver> MakeKrylov(const json &block, mfem::HypreBoomerAMG *prec) const
        {
            const std::string library = block.value("library", "mfem");
            const std::string type = block.value("type", "pcg");
            const double tol = block.value("tolerance", 1.0e-6);
            const int max_iter = block.value("max_iterations", 1000);

            const bool use_cg = (type == "pcg");
            if (!use_cg && type != "gmres")
                throw std::runtime_error("Unknown solver type: " + type + " (use pcg or gmres)");

            if (library == "mfem")
            {
                std::unique_ptr<mfem::IterativeSolver> solver;
                if (use_cg)
                    solver = std::make_unique<mfem::CGSolver>(comm_);
                else
                    solver = std::make_unique<mfem::GMRESSolver>(comm_);
                solver->SetRelTol(tol);
                solver->SetMaxIter(max_iter);
                solver->SetPrintLevel(0);
                if (prec)
                    solver->SetPreconditioner(*prec);
                return solver;
            }
            else if (library == "hypre")
            {
                if (use_cg)
                {
                    auto solver = std::make_unique<mfem::HyprePCG>(comm_);
                    solver->SetTol(tol);
                    solver->SetMaxIter(max_iter);
                    solver->SetPrintLevel(0);
                    if (prec)
                        solver->SetPreconditioner(*prec);
                    return solver;
                }
                auto solver = std::make_unique<mfem::HypreGMRES>(comm_);
                solver->SetTol(tol);
                solver->SetMaxIter(max_iter);
                solver->SetPrintLevel(0);
                if (prec)
                    solver->SetPreconditioner(*prec);
                return solver;
            }
            throw std::runtime_error("Unknown solver library: " + library + " (use mfem or hypre)");
        }

        // Hypre solvers do not report convergence; treat reaching max_iter as failure.
        // (Currently unused: the convergence checks in SolveDisplacement/SolveDamage are disabled.)
        static bool Converged(const mfem::Solver &solver, int max_iter)
        {
            if (auto *it = dynamic_cast<const mfem::IterativeSolver *>(&solver))
                return it->GetConverged();

            int iterations = 0;
            if (auto *pcg = dynamic_cast<const mfem::HyprePCG *>(&solver))
                pcg->GetNumIterations(iterations);
            else if (auto *gmres = dynamic_cast<const mfem::HypreGMRES *>(&solver))
                gmres->GetNumIterations(iterations);
            return iterations < max_iter;
        }

        // -------------------------------------------------------------------------------------------------------
        // Assembly methods
        // -------------------------------------------------------------------------------------------------------

        // Kdisp(dtemp), fdisp(dtemp) and the constrained system for utemp (essential values must already be in utemp).
        void AssembleDisplacement()
        {
            kdisp_->Update();
            kdisp_->Assemble(0);
            kdisp_->Finalize(0);

            fdisp_->Assemble();
            kdisp_->FormLinearSystem(bcs_->EssentialTrueDofs(), *utemp_, *fdisp_, *Kdisp_, u_new_, F_disp_);
        }

        // Kdmg1 + Kdmg2(utemp) and fdmg(utemp).
        void AssembleDamage()
        {
            kdmg2_->Update();
            kdmg2_->Assemble(0);
            kdmg2_->Finalize(0);

            Kdmg2Mat_.reset(kdmg2_->ParallelAssemble());
            Kdmg_.reset(mfem::Add(1.0, *Kdmg1Mat_, 1.0, *Kdmg2Mat_));

            fdmg_->Assemble();
            Fdmg_.reset(fdmg_->ParallelAssemble());
        }

        // -------------------------------------------------------------------------------------------------------
        // Solves
        // -------------------------------------------------------------------------------------------------------

        // Solve the elasticity equation: assemble, set up the solver and solve into utemp.
        bool SolveDisplacement()
        {
            disp_assembly_timer_.Start();
            AssembleDisplacement();
            disp_assembly_timer_.Stop();

            disp_solve_timer_.Start();
            if (dispPrec_)
                dispPrec_->SetOperator(*Kdisp_);

            dispKrylov_->SetOperator(*Kdisp_);
            dispKrylov_->Mult(F_disp_, u_new_);

            kdisp_->RecoverFEMSolution(u_new_, *fdisp_, *utemp_);
            disp_solve_timer_.Stop();

            // return Converged(*dispKrylov_, dispMaxIter_);
            return true;
        }

        bool SolveDamage()
        {
            dmg_assembly_timer_.Start();
            AssembleDamage();
            dmg_assembly_timer_.Stop();

            dmg_solve_timer_.Start();
            if (dmgPrec_)
                dmgPrec_->SetOperator(*Kdmg_);

            dmgKrylov_->SetOperator(*Kdmg_);
            dmgKrylov_->Mult(*Fdmg_, dNew_);
            dmg_solve_timer_.Stop();

            // return Converged(*dmgKrylov_, dmgMaxIter_);
            return true;
        }

        // One staggered solve at the given displacement condition.
        // damage irreversibity enforced via d = max(d_new, d_old).
        bool Stagger(double dispVal, int &iterations, double &error)
        {
            d_->GetTrueDofs(dOld_);
            dNew_ = dOld_; // initial guess for the damage solve
            error = 1.0;

            for (iterations = 1; iterations <= staggerMaxIter_; iterations++)
            {
                bcs_->Apply(dispVal, *utemp_);
                if (!SolveDisplacement())
                    return false;
                if (!SolveDamage())
                    return false;

                // Irreversibility; the error is the global max-norm of the damage update.
                double local = 0.0;
                for (int i = 0; i < dNew_.Size(); i++)
                {
                    const double dm = std::max(dNew_(i), dOld_(i));
                    local = std::max(local, dm - dOld_(i));
                    dOld_(i) = dm;
                }
                MPI_Allreduce(&local, &error, 1, MPI_DOUBLE, MPI_MAX, comm_);
                dtemp_->Distribute(dOld_);
                if (error <= staggerTol_)
                    return true;
            }
            iterations = staggerMaxIter_;
            return false;
        }

        void Accept() // temporary fields -> accepted fields
        {
            *u_ = *utemp_;
            *d_ = *dtemp_;
        }

        void Revert() // accepted fields -> temporary fields
        {
            *utemp_ = *u_;
            *dtemp_ = *d_;
        }

    public:
        PhaseFieldSolver(const json &input, mfem::ParMesh &pmesh, MPI_Comm comm)
            : input_(input), pmesh_(pmesh), comm_(comm) {}

        // Read fespaces.
        void ReadFESpaces(mfem::ParFiniteElementSpace *dispfespace_in, mfem::ParFiniteElementSpace *dmgfespace_in)
        {
            dispfes_ = dispfespace_in;
            dmgfes_ = dmgfespace_in;
        }

        // Read gridfunctions.
        void ReadFields(mfem::ParGridFunction *u_in, mfem::ParGridFunction *utemp_in, mfem::ParGridFunction *d_in, mfem::ParGridFunction *dtemp_in)
        {
            // Require(dispfes_, "FE spaces", "InitializeFESpaces()");
            u_ = u_in;
            utemp_ = utemp_in;
            d_ = d_in;
            dtemp_ = dtemp_in;
        }

        // Read forms
        void ReadForms(mfem::ParBilinearForm *kdisp_in, mfem::ParBilinearForm *kdmg1_in, mfem::ParBilinearForm *kdmg2_in, mfem::ParLinearForm *fdisp_in, mfem::ParLinearForm *fdmg_in)
        {
            // Require(utemp_, "Fields", "InitializeFields()");
            kdisp_ = kdisp_in;
            kdmg1_ = kdmg1_in;
            kdmg2_ = kdmg2_in;

            fdisp_ = fdisp_in;
            fdmg_ = fdmg_in;

            // Additionally, kdmg1 will be assembled since it only needs to be assembled once
            kdmg1_->Assemble(0);
            kdmg1_->Finalize(0);

            Kdmg1Mat_.reset(kdmg1_->ParallelAssemble());
        }

        // Dirichlet conditions and the displacement load range.
        void ReadBCs(BoundaryConditions *bcs_in)
        {
            // Require(dispfes_, "FE spaces", "InitializeFESpaces()");
            bcs_ = bcs_in;
        }

        // Preconditioners and Krylov solvers, chosen independently for displacement and damage.
        void InitializeSolvers()
        {
            // Require(kdisp_, "Forms", "InitializeForms()");
            // Require(bcs_, "Boundary conditions", "InitializeBCs()");

            const json staggered = input_.value("Staggered_Solve", json::object());
            staggerTol_ = staggered.value("tolerance", 1.0e-6);
            staggerMaxIter_ = staggered.value("max_iterations", 1000);

            const json &solvers = input_.at("Solvers");
            const json &disp_block = solvers.at("elasticity_solver");
            const json &dmg_block = solvers.at("damage_solver");

            dispMaxIter_ = disp_block.value("max_iterations", 1000);
            dmgMaxIter_ = dmg_block.value("max_iterations", 1000);

            dispPrec_ = MakePrec(disp_block.value("preconditioner", "amg_systems"), *dispfes_);
            dmgPrec_ = MakePrec(dmg_block.value("preconditioner", "amg"), *dmgfes_);

            dispKrylov_ = MakeKrylov(disp_block, dispPrec_.get());
            dmgKrylov_ = MakeKrylov(dmg_block, dmgPrec_.get());

            // warm start from the previous solution
            dispKrylov_->iterative_mode = true;
            dmgKrylov_->iterative_mode = true;

            Kdisp_ = std::make_unique<mfem::HypreParMatrix>();
        }

        // Accessors for the output writer
        mfem::ParFiniteElementSpace &DispFESpace() { return *dispfes_; }
        mfem::ParFiniteElementSpace &DamageFESpace() { return *dmgfes_; }
        mfem::ParGridFunction &Displacement() { return *u_; }
        mfem::ParGridFunction &Damage() { return *d_; }

        using SaveStep = std::function<void(int step, double load)>;

        // Runs load steps until the final displacement is reached; `onConvergence` is called after each converged
        // (staggered) step. Returns false if the load increment collapses without convergence.
        bool RunStaggeredSolve(const SaveStep &onConvergence)
        {
            // Require(Kdisp_, "Solvers", "InitializeSolvers()");

            const double incFloor = 1e-3 * bcs_->minInc_;
            double dispVal = 0.0, dispInc = std::min(2.0 * bcs_->minInc_, bcs_->maxInc_);
            int step = 0, substep = 1;

            while (dispVal < bcs_->finalDisp_ * (1.0 - 1e-12))
            {
                // Displacement increment for this step
                dispInc = std::min(dispInc, bcs_->finalDisp_ - dispVal);

                int iterations = 0;
                double damage_error = 0.0;

                // One stagger step: elasticity is first solved, then damage, until the damage update is small.
                const bool converged = Stagger(dispVal + dispInc, iterations, damage_error);

                if (converged)
                {
                    Accept(); // Updated grid functions
                    dispVal += dispInc;
                    step++;
                    substep = 1;

                    if (mfem::Mpi::Root())
                    {
                        std::cout << "Converged. Step " << step << " Substep " << substep << " iterations " << iterations << " dispVal " << dispVal << " dispInc " << dispInc << std::endl;
                    }

                    onConvergence(step, dispVal);

                    // Adapt the displacement increment to the number of stagger iterations.
                    if (step > 3)
                    {
                        if (iterations <= 1)
                            dispInc *= 2.0;
                        else if (iterations < 10)
                            dispInc *= 1.5;
                        else if (iterations < 30)
                            dispInc *= 1.15;
                        else if (iterations >= 2000)
                            dispInc /= 2.0;
                        else if (iterations >= 1000)
                            dispInc /= 1.5;
                        dispInc = std::min(dispInc, bcs_->maxInc_);
                    }
                }
                else
                {
                    // No convergence: restore the fields and retry with a smaller increment.
                    Revert();
                    substep++;

                    dispInc /= (substep < 4) ? 2.0 : (substep < 10 ? 4.0 : 10.0);

                    if (mfem::Mpi::Root())
                    {
                        std::cout << "reject at dispVal " << dispVal << " (damage_error " << damage_error << ", iterations " << iterations
                                  << "), dispInc -> " << dispInc << std::endl;
                    }

                    // Increment collapsed: give up.
                    if (dispInc < incFloor)
                        return false;
                }
            }

            if (mfem::Mpi::Root())
            {
                std::cout << "disp assembly time: " << disp_assembly_timer_.RealTime() << std::endl;
                std::cout << "disp solve time: " << disp_solve_timer_.RealTime() << std::endl;
                std::cout << "dmg assembly time: " << dmg_assembly_timer_.RealTime() << std::endl;
                std::cout << "dmg solve time: " << dmg_solve_timer_.RealTime() << std::endl;
            }
            return true;
        }
    };

    //------------------------------------------------------------------------------------------------------------------
    // Result writer: ADIOS2 BP5 (default) or ParaView, plus the resolved input and run details.
    //------------------------------------------------------------------------------------------------------------------
    class OutputWriter
    {
    private:
        mfem::ParMesh &mesh_;
        mfem::ParGridFunction &u_, &d_;
        std::unique_ptr<mfem::adios2stream> aos_;
        std::unique_ptr<mfem::ParaViewDataCollection> pv_;

    public:
        OutputWriter(const json &input, mfem::ParMesh &mesh, PhaseFieldSolver &solver, MPI_Comm comm)
            : mesh_(mesh), u_(solver.Displacement()), d_(solver.Damage())
        {
            int myid, nranks;
            MPI_Comm_rank(comm, &myid);
            MPI_Comm_size(comm, &nranks);

            const json output = input.value("Output", json::object());
            const std::string name = input.at("Test_Name");
            const std::string format = output.value("format", "adios2");

            namespace fs = std::filesystem;
            const fs::path dir = fs::path(output.value("path", "../results")) / name;
            if (myid == 0)
                fs::create_directories(dir);
            MPI_Barrier(comm);

            // Global sizes are collective calls, so every rank evaluates them.
            const long long num_elements = mesh.GetGlobalNE();
            const long long num_disp_dofs = solver.DispFESpace().GlobalTrueVSize();
            const long long num_dmg_dofs = solver.DamageFESpace().GlobalTrueVSize();

            if (myid == 0)
            {
                std::ofstream(dir / "input.json") << input.dump(2) << "\n";
                std::ofstream details(dir / "SimulationDetails.txt");
                details << "Test_Name\t" << name << "\n"
                        << "mpi_ranks\t" << nranks << "\n"
                        << "mesh\t" << input.at("Mesh").value("file", "") << "\n"
                        << "mesh_order\t" << input.at("Mesh").value("order", 1) << "\n"
                        << "refinement_levels\t" << input.at("Mesh").value("refinement_levels", 0) << "\n"
                        << "global_elements\t" << num_elements << "\n"
                        << "disp_dofs\t" << num_disp_dofs << "\n"
                        << "damage_dofs\t" << num_dmg_dofs << "\n";
            }

            if (format == "adios2")
            {
                aos_ = std::make_unique<mfem::adios2stream>((dir / "FractureMFEMoutput.bp").string(), mfem::adios2stream::openmode::out, comm, "BP5");
                aos_->SetRefinementLevel(0);
                aos_->SetParameter("FlushStepsCount", "1");
            }
            else if (format == "paraview")
            {
                pv_ = std::make_unique<mfem::ParaViewDataCollection>("FractureMFEMoutput", &mesh);
                pv_->SetPrefixPath(dir.string());
                pv_->SetHighOrderOutput(true);
                pv_->RegisterField("disp", &u_);
                pv_->RegisterField("damage", &d_);
            }
            else
                throw std::runtime_error("Unknown Output.format: " + format);
        }

        // Writes the accepted fields u and d.
        void Save(int cycle, double time)
        {
            if (aos_)
            {
                aos_->BeginStep();
                mesh_.Print(*aos_);
                aos_->SetTime(time);
                aos_->SetCycle(cycle);
                u_.Save(*aos_, "disp", mfem::adios2stream::data_type::point_data);
                d_.Save(*aos_, "damage", mfem::adios2stream::data_type::point_data);
                aos_->EndStep();
            }
            else
            {
                pv_->SetCycle(cycle);
                pv_->SetTime(time);
                pv_->Save();
            }
        }
    };

} // namespace spicule
