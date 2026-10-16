#include "spicule/PhaseFieldSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace mfem;

namespace spicule
{
    namespace
    {
        std::unique_ptr<HypreBoomerAMG> MakePrec(const std::string &kind, ParFiniteElementSpace *parfespace)
        {
            // no preconditioner
            if (kind == "none")
                return nullptr;

            auto amg = std::make_unique<HypreBoomerAMG>();
            amg->SetPrintLevel(0);
            if (kind == "amg_systems")
                amg->SetSystemsOptions(parfespace->GetVDim(), true); // byNODES
            else if (kind == "amg_elasticity")
                amg->SetElasticityOptions(parfespace); // byVDIM
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

        // Krylov solver from the input: library = mfem | hypre, type = pcg | gmres.
        // Note the stopping criteria differ: mfem solvers use `tol` relative to the initial
        // preconditioned residual, hypre solvers relative to ||b||.
        std::unique_ptr<Solver> MakeKrylov(const LinearSolverParams &params, HypreBoomerAMG *prec, MPI_Comm comm)
        {
            const bool use_cg = (params.type == "pcg");
            if (!use_cg && params.type != "gmres")
                throw std::runtime_error("Unknown solver type: " + params.type + " (use pcg or gmres)");

            if (params.library == "mfem")
            {
                std::unique_ptr<IterativeSolver> solver;
                if (use_cg)
                    solver = std::make_unique<CGSolver>(comm);
                else
                    solver = std::make_unique<GMRESSolver>(comm);

                solver->SetRelTol(params.tol);
                solver->SetMaxIter(params.maxIter);
                solver->SetPrintLevel(0);
                solver->iterative_mode = true;
                if (prec)
                    solver->SetPreconditioner(*prec);
                return solver;
            }
            else if (params.library == "hypre")
            {
                if (use_cg)
                {
                    auto solver = std::make_unique<HyprePCG>(comm);
                    solver->SetTol(params.tol);
                    solver->SetMaxIter(params.maxIter);
                    solver->SetPrintLevel(0);
                    if (prec)
                        solver->SetPreconditioner(*prec);
                    return solver;
                }
                else
                {
                    auto solver = std::make_unique<HypreGMRES>(comm);
                    solver->SetTol(params.tol);
                    solver->SetMaxIter(params.maxIter);
                    solver->SetPrintLevel(0);
                    if (prec)
                        solver->SetPreconditioner(*prec);
                    return solver;
                }
            }
            throw std::runtime_error("Unknown solver library: " + params.library + " (use mfem or hypre)");
        }

        // Hypre solvers do not report convergence; treat reaching max_iter as failure.
        bool Converged(const Solver &solver, int maxIter)
        {
            if (auto *it = dynamic_cast<const IterativeSolver *>(&solver))
                return it->GetConverged();

            int iterations = 0;
            if (auto *pcg = dynamic_cast<const HyprePCG *>(&solver))
                pcg->GetNumIterations(iterations);
            else if (auto *gmres = dynamic_cast<const HypreGMRES *>(&solver))
                gmres->GetNumIterations(iterations);
            return iterations < maxIter;
        }
    } // unnamed namespace

    // Constructor
    PhaseFieldSolver::PhaseFieldSolver(const Config &config, FESpaces &spaces, Fields &fields, Forms &forms, const BoundaryConditions &bcs, MPI_Comm comm)
        : config_(config), spaces_(spaces), fields_(fields), forms_(forms), bcs_(bcs), comm_(comm)
    {
        // Preconditioners
        dispPrec_ = MakePrec(config.dispSolver.preconditioner, spaces_.dispfespace.get());
        dmgPrec_ = MakePrec(config.dmgSolver.preconditioner, spaces_.damagefespace.get());

        // Krylov solvers, chosen independently for displacement and damage
        dispKrylov_ = MakeKrylov(config.dispSolver, dispPrec_.get(), comm);
        dmgKrylov_ = MakeKrylov(config.dmgSolver, dmgPrec_.get(), comm);

        // warm start from the previous solution
        dispKrylov_->iterative_mode = true;
        dmgKrylov_->iterative_mode = true;

        // hypre par matrix
        Kdisp_ = std::make_unique<HypreParMatrix>();
    }

    // Solve elasticity equation. Set up stiffness matrix, linear system and solve.
    bool PhaseFieldSolver::SolveDisplacement()
    {
        disp_assembly_timer.Start();
        forms_.AssembleDisplacement(bcs_.EssentialTrueDofs(), *Kdisp_, u_new, F_disp);
        disp_assembly_timer.Stop();
        // Kdisp_ = K_op.As<HypreParMatrix>();

        disp_solve_timer.Start();
        if (dispPrec_)
            dispPrec_->SetOperator(*Kdisp_);

        // set operator and mult
        dispKrylov_->SetOperator(*Kdisp_);
        dispKrylov_->Mult(F_disp, u_new);

        forms_.RecoverDisplacement(u_new);
        disp_solve_timer.Stop();

        // return Converged(*dispKrylov_, config_.dispSolver.maxIter);
        return 1;
    }

    bool PhaseFieldSolver::SolveDamage(Vector &dNew)
    {
        dmg_assembly_timer.Start();
        Kdmg_ = forms_.AssembleDamageMatrix();
        Fdmg_ = forms_.AssembleDamageRhs();
        dmg_assembly_timer.Stop();

        dmg_solve_timer.Start();
        if (dmgPrec_)
            dmgPrec_->SetOperator(*Kdmg_);

        // set operator
        dmgKrylov_->SetOperator(*Kdmg_);

        dmgKrylov_->Mult(*Fdmg_, dNew);
        dmg_solve_timer.Stop();

        // return Converged(*dmgKrylov_, config_.dmgSolver.maxIter);
        return 1;
    }

    bool PhaseFieldSolver::Stagger(double dispVal, int &iterations, double &error)
    {
        fields_.d.GetTrueDofs(dOld);
        dNew = dOld; // initial guess for the damage solve
        error = 1.0;

        for (iterations = 1; iterations <= config_.staggered.maxIter; iterations++)
        {
            bcs_.Apply(dispVal, fields_.utemp);
            if (!SolveDisplacement())
                return false;
            if (!SolveDamage(dNew))
                return false;

            // Irreversibility: d = max(d_new, d_old). Error is the global max-norm of the update.
            double local = 0.0;
            for (int i = 0; i < dNew.Size(); i++)
            {
                const double dm = std::max(dNew(i), dOld(i));
                local = std::max(local, dm - dOld(i));
                dOld(i) = dm;
            }
            MPI_Allreduce(&local, &error, 1, MPI_DOUBLE, MPI_MAX, comm_);
            fields_.dtemp.Distribute(dOld);
            if (error <= config_.staggered.tol)
                return true;
        }
        iterations = config_.staggered.maxIter;
        return false;
    }

    bool PhaseFieldSolver::RunStaggeredSolve(const SaveStep &onConvergence)
    {
        const LoadParams &load_params = config_.load;
        const double incFloor = 1e-3 * load_params.minInc;
        double dispVal = 0.0, dispInc = std::min(2.0 * load_params.minInc, load_params.maxInc);

        int step = 0, substep = 1;

        // While dispval is less the final disp
        while (dispVal < load_params.finalDisp * (1.0 - 1e-12))
        {
            // Determine displacement increment for this step
            dispInc = std::min(dispInc, load_params.finalDisp - dispVal);

            int iterations = 0;
            double damage_error = 0.0;

            // One stagger step. Solves elasticity problem and then the damage problem
            const bool converged = Stagger(dispVal + dispInc, iterations, damage_error);

            // Update fields and save if converged.
            if (converged)
            {
                fields_.Update();
                dispVal += dispInc;
                step++;
                substep = 1;

                if (Mpi::Root())
                {
                    std::cout << "Converged. Step " << step << " Substep " << substep << " iterations " << iterations << std::endl;
                    std::cout << " dispVal " << dispVal << " dispInc " << dispInc << std::endl;
                }

                // Save step
                onConvergence(step, dispVal);

                // Some conditions to adaptively determine dispVal steps.
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
                    dispInc = std::min(dispInc, load_params.maxInc);
                }
            }
            else
            {
                // No convergence. Revert fields to try new dispInc.
                fields_.Revert();
                substep++;

                dispInc /= (substep < 4) ? 2.0 : (substep < 10 ? 4.0 : 10.0);

                if (Mpi::Root())
                {
                    std::cout << "reject at dispVal " << dispVal << " (damage_error " << damage_error << ", iterations " << iterations
                              << "), dispInc -> " << dispInc << std::endl;
                }

                // If dispInc has reduced too much, simulation needs to be killed, so return false.
                if (dispInc < incFloor)
                    return false;
            }
        }
        if (Mpi::Root())
        {
            std::cout << "disp assembly time: " << disp_assembly_timer.RealTime() << std::endl;
            std::cout << "disp solve time: " << disp_solve_timer.RealTime() << std::endl;
            std::cout << "dmg assembly time: " << dmg_assembly_timer.RealTime() << std::endl;
            std::cout << "dmg solve time: " << dmg_solve_timer.RealTime() << std::endl;
        }
        return true;
    }

} // namespace spicule
