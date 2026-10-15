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
        std::unique_ptr<IterativeSolver> MakeKrylov(const LinearSolverParams &params, MPI_Comm comm)
        {
            std::unique_ptr<IterativeSolver> solver;

            if (params.type == "pcg")
                solver = std::make_unique<CGSolver>(comm);
            else if (params.type == "gmres")
                solver = std::make_unique<GMRESSolver>(comm);
            else
                throw std::runtime_error("Unknown solver type: " + params.type);

            solver->SetRelTol(params.tol);
            solver->SetMaxIter(params.maxIter);
            solver->SetPrintLevel(0);

            return solver;
        }

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
                amg->SetElasticityOptions(parfespace); // byNODES
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
    } // unnamed namespace

    // Constructor
    PhaseFieldSolver::PhaseFieldSolver(const Config &config, FESpaces &spaces, Fields &fields, Forms &forms, const BoundaryConditions &bcs, MPI_Comm comm)
        : config_(config), spaces_(spaces), fields_(fields), forms_(forms), bcs_(bcs), comm_(comm)
    {
        // Preconditioners
        dispPrec_ = MakePrec(config.dispSolver.preconditioner, spaces_.dispfespace.get());
        dmgPrec_ = MakePrec(config.dmgSolver.preconditioner, spaces_.dispfespace.get());

        // Krylov solvers
        dispKrylov_ = MakeKrylov(config.dispSolver, comm);
        dmgKrylov_ = MakeKrylov(config.dmgSolver, comm);

        if (dispPrec_)
            dispKrylov_->SetPreconditioner(*dispPrec_);
        if (dmgPrec_)
            dmgKrylov_->SetPreconditioner(*dmgPrec_);

        // warm start from the previous displacement
        dispKrylov_->iterative_mode = true;
        dmgKrylov_->iterative_mode = true;
    }

    // Solve elasticity equation. Set up stiffness matrix, linear system and solve.
    bool PhaseFieldSolver::SolveDisplacement()
    {
        forms_.AssembleDisplacement(bcs_.EssentialTrueDofs(), K_op, u_new, F_disp);

        K = K_op.As<HypreParMatrix>();

        if (dispPrec_)
            dispPrec_->SetOperator(*K);

        // set operator and mult
        dispKrylov_->SetOperator(*K);
        dispKrylov_->Mult(F_disp, u_new);

        forms_.RecoverDisplacement(u_new);

        return dispKrylov_->GetConverged();
    }

    bool PhaseFieldSolver::SolveDamage(Vector &dNew)
    {
        Kdmg_ = forms_.AssembleDamageMatrix();
        Fdmg_ = forms_.AssembleDamageRhs();

        if (dmgPrec_)
            dmgPrec_->SetOperator(*Kdmg_);

        // set operator
        dmgKrylov_->SetOperator(*Kdmg_);

        dmgKrylov_->Mult(*Fdmg_, dNew);

        return dmgKrylov_->GetConverged();
    }

    bool PhaseFieldSolver::Stagger(double dispVal, int &iterations, double &error)
    {
        fields_.d.GetTrueDofs(dOld);
        dNew.SetSize(dOld.Size());
        dNew = 0.0;
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
        return true;
    }

} // namespace spicule
