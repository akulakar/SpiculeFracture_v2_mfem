#pragma once
// Regularized variational fracture (with AT2) bilinear and linear forms built from mfemplus fracture integrators.
//
//   bilinear:  Kdisp(d)   damaged isotropic elasticity            (reads dtemp)
//              Kdmg1      Gc/eps mass + Gc*eps diffusion          (constant, assembled once)
//              Kdmg2(u)/Kdmg2(quad)   strain-energy driven reaction term      (reads utemp or quadrature point energy)
//
//   linear:    fdisp(d)   volumetric eigenstrain body force
//              fdmg(u)/fdmg(quad)    strain-energy driving force

#include <memory>
#include <mfem.hpp>

#include "spicule/Config.hpp"
#include "spicule/Setup.hpp"

namespace spicule
{
    class Forms
    {

    private:
        FESpaces &fe_spaces_;
        Fields &fields_;
        mfem::ConstantCoefficient Ey_, nu_, volEig_, gOverEps_, gTimesEps_;
        // mfem::FunctionCoefficient volEigFun_;

        std::unique_ptr<mfem::ParBilinearForm> kdisp_, kdmg1_, kdmg2_;
        std::unique_ptr<mfem::ParLinearForm> fdisp_, fdmg_;
        std::unique_ptr<mfem::HypreParMatrix> Kdmg1Mat_, Kdmg2Mat_;

    public:
        Forms(const Config &config, FESpaces &fe_spaces, Fields &fields) : fe_spaces_(fe_spaces), fields_(fields),
                                                                           Ey_(config.material.Ey), nu_(config.material.nu),
                                                                           volEig_(config.eigenstrain.beta * config.eigenstrain.deltaC),
                                                                           gOverEps_(config.material.Gc / config.material.eps),
                                                                           gTimesEps_(config.material.Gc * config.material.eps) {};

        void AddIntegrators(const Config &config);
        // Assembles Kdisp and fdisp for the current dtemp and forms the lifted linear system
        // for utemp (essential values must already be in utemp).
        // with vector class
        void AssembleDisplacement(const mfem::Array<int> &essTdofs, mfem::OperatorHandle &A, mfem::Vector &X, mfem::Vector &B);
        void RecoverDisplacement(const mfem::Vector &X);
        // withe hypre par matrix
        void AssembleDisplacement(const mfem::Array<int> &essTdofs, mfem::HypreParMatrix &K, mfem::Vector &X, mfem::Vector &B);

        // Kdmg1 + Kdmg2(utemp) and fdmg(utemp). Caller owns the results.
        std::unique_ptr<mfem::HypreParMatrix> AssembleDamageMatrix();
        std::unique_ptr<mfem::HypreParVector> AssembleDamageRhs();
    };

} // namespace spicule
