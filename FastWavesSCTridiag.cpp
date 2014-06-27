#include "StencilFramework.h"
#include "DycoreConfiguration.h"
#include "FastWavesSCTridiag.h"

enum
{
    a, b, c, rhs, tmp, bet, y
};

namespace FastWavesSCTridiagStages
{
    template <typename TEnv>
    struct Forward
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, y  ) // Output
        STAGE_PARAMETER(FullDomain, a  ) // Input
        STAGE_PARAMETER(FullDomain, b  ) // Input
        STAGE_PARAMETER(FullDomain, c  ) // Input
        STAGE_PARAMETER(FullDomain, rhs) // Input
        STAGE_PARAMETER(FullDomain, tmp) // Output buffer
        STAGE_PARAMETER(FullDomain, bet) // Output buffer

        __ACC__
        static void Do(Context ctx, KMinimumCenter)
        {
            ctx[bet::Center()] = ctx[b::Center()];
            ctx[y::Center()] = ctx[rhs::Center()] / ctx[bet::Center()];
        }

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[tmp::Center()] = ctx[c::At(kminus1)] / ctx[bet::Center()];
            ctx[bet::Center()] = ctx[b::Center()] - ctx[a::Center()] * ctx[tmp::Center()];
            ctx[y::Center()] = (ctx[rhs::Center()] - ctx[a::Center()] * ctx[y::At(kminus1)]) / ctx[bet::Center()];
        }
    };

    template <typename TEnv>
    struct Backward
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, y  ) // Output
        STAGE_PARAMETER(FullDomain, tmp) // Input buffer

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[y::Center()] -= ctx[tmp::At(kplus1)] * ctx[y::At(kplus1)];
        }
    };

} // End namespace FastWavesSCTridiagStages

FastWavesSCTridiag::FastWavesSCTridiag() {}
FastWavesSCTridiag::~FastWavesSCTridiag() {}

void FastWavesSCTridiag::Init(const DycoreConfiguration& dycoreConfiguration, DycoreRepository& dycoreRepository, TimeIntegratorRepository& timeIntegratorRepository, FastWavesSCRepository& fastWavesSCRepository)
{
    using namespace FastWavesSCTridiagStages;

    StencilCompiler::Build(
        stencil_,
        "FastWavesSCTridiag",
        dycoreRepository.calculationDomain(),
        StencilConfiguration<Real, FastWavesSCTridiagBlockSize>(),
        pack_parameters(
            /* output fields */
            Param<y, cInOut>(dycoreRepository.w_in()),
            /* input fields */
            Param<a  , cIn>(fastWavesSCRepository.lgs_a()  ),
            Param<b  , cIn>(fastWavesSCRepository.lgs_b()  ),
            Param<c  , cIn>(fastWavesSCRepository.lgs_c()  ),
            Param<rhs, cIn>(fastWavesSCRepository.lgs_rhs())
        ),
        define_temporaries(
            StencilBuffer<tmp, Real, KRange<FullDomain,0,1> >(),
            StageVariable<bet, Real>()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<Forward, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,1> >()
                )
            )
            ,
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<Backward, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}

void FastWavesSCTridiag::Apply()
{
    stencil_.Apply();
}

