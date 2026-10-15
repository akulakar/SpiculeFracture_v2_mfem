#pragma once
// Staggered (alternate minimisation) solver with adaptive displacement-load stepping.

#include <functional>
#include <memory>
#include <mfem.hpp>

#include "spicule/Config.hpp"
#include "spicule/Forms.hpp"
#include "spicule/Setup.hpp"

namespace spicule
{
    class PhaseFieldSolver
    {
    private:
        const Config &config_;
        FESpaces &spaces_;
        Fields &fields_;
        Forms &forms_;
        const BoundaryConditions &bcs_;
        MPI_Comm comm_;

        // Persistent solvers
        std::unique_ptr<mfem::HypreBoomerAMG> dispPrec_, dmgPrec_;
        std::unique_ptr<mfem::IterativeSolver> dispKrylov_, dmgKrylov_;

        // Operators are swapped each iteration
        std::unique_ptr<mfem::HypreParMatrix> Kdmg_;
        std::unique_ptr<mfem::HypreParVector> Fdmg_;

        // Elasticity
        mfem::HypreParMatrix *K;
        mfem::OperatorHandle K_op;
        mfem::Vector u_new, F_disp;

        // Damage
        mfem::Vector dOld, dNew;

    public:
        PhaseFieldSolver(const Config &cfg, FESpaces &spaces, Fields &fields, Forms &forms, const BoundaryConditions &bcs, MPI_Comm comm);

        using SaveStep = std::function<void(int step, double load)>;

        // Main method that runs load steps until final load step is reached. converged is true after each converged (staggered) step. Returns false if no convergence.
        bool RunStaggeredSolve(const SaveStep &onConvergence);

    private:
        // One staggered solve at displacement_value; irreversible damage via d = max(d_new, d_old).
        bool Stagger(double displacement, int &iters, double &err);
        bool SolveDisplacement();
        bool SolveDamage(mfem::Vector &dNew);
    };

} // namespace spicule
