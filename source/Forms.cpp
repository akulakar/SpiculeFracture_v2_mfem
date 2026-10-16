// some comments

#include "spicule/Forms.hpp"
#include <mfemplus/SolidMechanics/mfemplus_solids.hpp>

using namespace mfem;

namespace spicule
{
    void Forms::AddIntegrators(const Config &config)
    {
        ParFiniteElementSpace *dispFes = fe_spaces_.dispfespace.get();
        ParFiniteElementSpace *dmgFes = fe_spaces_.damagefespace.get();
        const double Keps = config.material.Kepsilon;

        // Bilinear form for elasticity
        kdisp_ = std::make_unique<ParBilinearForm>(dispFes);
        kdisp_->AddDomainIntegrator(new mfemplus::DamagedIsotropicElasticityIntegrator(Ey_, nu_, Keps, fields_.dtemp, dmgFes));

        // Bilinear form for damage equation. Split into to because only one is reassembled every iteration
        kdmg1_ = std::make_unique<ParBilinearForm>(dmgFes);
        kdmg1_->AddDomainIntegrator(new MassIntegrator(gOverEps_));
        kdmg1_->AddDomainIntegrator(new DiffusionIntegrator(gTimesEps_));

        kdmg2_ = std::make_unique<ParBilinearForm>(dmgFes);
        kdmg2_->AddDomainIntegrator(new mfemplus::IsotropicStrainEnergyVolumetricEigenStrainDamageIntegrator(Ey_, nu_, volEig_, fields_.utemp, dispFes));

        // Linear form for elasticity equation.
        fdisp_ = std::make_unique<ParLinearForm>(dispFes);
        fdisp_->AddDomainIntegrator(new mfemplus::DamagedIsotropicVolumetricEigenstrainBodyForceLFIntegrator(Ey_, nu_, volEig_, Keps, fields_.dtemp, dmgFes));

        // Linear form for damage equation.
        fdmg_ = std::make_unique<ParLinearForm>(dmgFes);
        fdmg_->AddDomainIntegrator(new mfemplus::IsotropicStrainEnergyVolumetricEigenStrainDamageLFIntegrator(Ey_, nu_, volEig_, fields_.utemp, dispFes));

        // Constant part of the damage operator.
        kdmg1_->Assemble(0);
        kdmg1_->Finalize(0);
        // Assemble once.
        Kdmg1Mat_.reset(kdmg1_->ParallelAssemble());
    }

    void Forms::AssembleDisplacement(const Array<int> &essTdofs, OperatorHandle &A, Vector &X, Vector &B)
    {
        kdisp_->Update();
        kdisp_->Assemble(0);
        kdisp_->Finalize(0);

        fdisp_->Assemble();
        // Boundary conditions need to be projected on utemp
        kdisp_->FormLinearSystem(essTdofs, fields_.utemp, *fdisp_, A, X, B);
    }

    void Forms::AssembleDisplacement(const Array<int> &essTdofs, HypreParMatrix &K, Vector &X, Vector &B)
    {
        kdisp_->Update();
        kdisp_->Assemble(0);
        kdisp_->Finalize(0);

        fdisp_->Assemble();
        // Boundary conditions need to be projected on utemp
        kdisp_->FormLinearSystem(essTdofs, fields_.utemp, *fdisp_, K, X, B);
    }

    void Forms::RecoverDisplacement(const Vector &X)
    {
        kdisp_->RecoverFEMSolution(X, *fdisp_, fields_.utemp);
    }

    std::unique_ptr<HypreParMatrix> Forms::AssembleDamageMatrix()
    {
        kdmg2_->Update();
        kdmg2_->Assemble(0);
        kdmg2_->Finalize(0);
        Kdmg2Mat_.reset(kdmg2_->ParallelAssemble());
        return std::unique_ptr<HypreParMatrix>(Add(1.0, *Kdmg1Mat_, 1.0, *Kdmg2Mat_));
    }

    std::unique_ptr<HypreParVector> Forms::AssembleDamageRhs()
    {
        fdmg_->Assemble();
        return std::unique_ptr<HypreParVector>(fdmg_->ParallelAssemble());
    }

} // namespace spicule
