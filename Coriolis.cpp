#include <cassert>
#include "StencilFramework.h"
#include "FiniteDifferenceFunctions.h"
#include "Coriolis.h"

// define functions

/**
* @struct CoriolisForce
* Function stencil calculating the Coriolis force using the following parameters:
* 1. fc constant Coriolis force factor
* 2. vel velocity used to calculate the force
*/
template<typename TEnv>
struct CoriolisForce
{
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, fc) 
    FUNCTION_PARAMETER(1, vel)

    __ACC__
    static T Do(Context ctx)
    {
        return ctx[fc::Center()] * ctx[vel::Center()];
    }
};

// define parameter enum
enum 
{ 
    utens, vtens, 
    u, v, fc 
};

namespace CoriolisStages
{
    /**
    * @struct USlowTensStage
    * Stage calculating the u tendency update due to Coriolis force
    */
    template<typename TEnv>
    struct USlowTensStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, utens)
        STAGE_PARAMETER(FullDomain, v)
        STAGE_PARAMETER(FullDomain, fc)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[utens::Center()] += 
                ctx[Call<Average>::With(jminus1, Call<CoriolisForce>::With(fc::Center(), Call<Average>::With(iplus1, v::Center())))];
        }
    };

    /**
    * @struct VSlowTensStage
    * Stage calculating the v tendency update due to Coriolis force
    */
    template<typename TEnv>
    struct VSlowTensStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, vtens)
        STAGE_PARAMETER(FullDomain, u)
        STAGE_PARAMETER(FullDomain, fc)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[vtens::Center()] -= 
                ctx[Call<Average>::With(iminus1, Call<CoriolisForce>::With(fc::Center(), Call<Average>::With(jplus1, u::Center()))) ];
        }
    };
}

Coriolis::Coriolis() {}
Coriolis::~Coriolis() {}

void Coriolis::Init(DycoreRepository& dycoreRepository)
{
    using namespace CoriolisStages;
    // initialize the stencil
    StencilCompiler::Build(
        stencil_,
        "Coriolis",
        dycoreRepository.calculationDomain(),
        StencilConfiguration<Real, CoriolisBlockSize>(),
        pack_parameters(
            /* output fields */
            Param<utens, cInOut>(dycoreRepository.utens()),
            Param<vtens, cInOut>(dycoreRepository.vtens()),
            /* input fields */
            Param<u, cIn>(dycoreRepository.u_nnow()),
            Param<v, cIn>(dycoreRepository.v_nnow()),
            Param<fc, cIn>(dycoreRepository.fc())
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<USlowTensStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<VSlowTensStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}

void Coriolis::Do() 
{ 
    stencil_.Apply(); 
}
    
Stencil& Coriolis::stencil()
{
    return stencil_;
}
  
