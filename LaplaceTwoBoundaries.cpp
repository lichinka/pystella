/*
 * 05. LaplaceTwoBoundaries
 *
 * This tutorial shows the solution of the heat equation by using a time stepping
 * technique. The boundary conditions are enforced through the HaloUpdater
 * framework and consist of an heat inflow, an heat outflow and two isolating walls
 * (i.e. two value boundary condition and two zero-gradient boundary conditions).
 */

// Include the framework
#include "StencilFramework.h"
#include "HaloUpdateFramework.h"

// Other includes
#include "stencilsUtils.h"
#include <iostream>
#include <iomanip>
#include "check.h"

/*******************************
 * IMPLEMENTATION OF THE STAGE *
 *******************************/

// Enum declaring the parameter
enum { dataIn, dataOut };

// Define stages inside a namespace to avoid conflicts with other stencils
namespace LaplaceTwoBoundariesStages
{
    template<typename TEnv>
    struct LaplaceStage 
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, dataIn)
        STAGE_PARAMETER(FullDomain, dataOut)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[dataOut::Center()] = ctx[dataIn::Center()] + 0.2 * (
                    -4.*ctx[dataIn::Center()]
                         + ctx[dataIn::At(iplus1)] + ctx[dataIn::At(iminus1)]
                         + ctx[dataIn::At(jplus1)] + ctx[dataIn::At(jminus1)]
              );
        }
    };
} 

/*****************************
 * DEFINITION OF THE STENCIL *
 *****************************/

class LaplaceStencil
{
public:
    LaplaceStencil() : westValue_(0.), eastValue_(0.) { }

    void Init(IJKRealField& fieldIn, IJKRealField& fieldOut,
              Real westValue, Real eastValue)
    {
        using LaplaceTwoBoundariesStages::LaplaceStage;

        StencilCompiler::Build(
            stencil_,
            "LaplaceStencil",
            fieldIn.calculationDomain(),
            StencilConfiguration<Real, BlockSize<4, 4> >(),
            pack_parameters(
                Param<dataIn, cIn>(fieldIn),
                Param<dataOut, cInOut>(fieldOut)
            ),
            concatenate_sweeps(
                define_sweep<cKIncrement>(
                    define_stages(
                        StencilStage<
                            LaplaceStage,
                            IJRange<cComplete,0,0,0,0>,
                            KRangeFullDomain
                        >()
                    )
                )
            )
        );


        // Temporary objects for boundary conditions:
        // * inner and outer boundaries
        // * BoundaryCondition
        IJBoundary innerBoundary, outerBoundary;
        
        // We need 3 different boundary condition objects:
        //   At the west: a value boundary condition that we use to set a fix value at the west boundary
        //   At the east: a value boundary condition that we use to set a fix value at the east boundary
        //   At the North/South: a zero gradient boundary condition to describe an isolation condition
        // Declare these 3 different objects with the corresponding type:
        //   valueBoundaryWest, valueBoundaryEast and isolationBoundary
        IJKRealValueBoundaryCondition valueBoundaryWest, valueBoundaryEast;
        IJKRealZeroGradientBoundaryCondition isolationBoundary;


        // Configuration initialization:
        // * non-periodic boundary condition (2*false)
        // * condition applied on the global boundary (4*true)
        haloConfiguration_.Init(false, false, true, true, true, true);

        // Initialization of the halo updater
        haloUpdate_.Init("haloUpdate", haloConfiguration_);

        // We have to initialize the 3 boundary conditions we need to use in this example.
        // Remember that for ValueBoundaryConditions the value used to fill the boundaries,
        // has to exist until the end of the computation, therefore we have to store this value
        // in a property of our class.
        // First set the westValue_ and eastValue_ to values passed by argument
        // to this Init method.
        // Then initialize the two value boundary conditions (with the corresponding values) and the
        // zero gradient boundary condition (which represents an isolation condition).
        westValue_ = westValue;
        valueBoundaryWest.Init(fieldIn, westValue_);

        eastValue_ = eastValue;
        valueBoundaryEast.Init(fieldIn, eastValue_);

        isolationBoundary.Init(fieldIn);

        // We need to describe with inner & outer boundaries the west boundary condition.
        // Initialize innerBoundary and outerBoundar with the proper offsets to describe a halo of 1 line
        //   at the west boundary.
        // Add the job with the valueBoundaryWest and the current inner & outer boundaries to the haloUpdate.
        innerBoundary.Init(0,0,0,0);
        outerBoundary.Init(-1,0,0,0);
        haloUpdate_.AddJob(valueBoundaryWest, innerBoundary, outerBoundary);

        // We need to describe with inner & outer boundaries the east boundary condition.
        // Initialize innerBoundary and outerBoundar with the proper offsets to describe a halo of 1 line
        //   at the east boundary.
        // Add the job with the valueBoundaryEast and the current inner & outer boundaries to the haloUpdate.
        innerBoundary.Init(0,0,0,0);
        outerBoundary.Init(0,1,0,0);
        haloUpdate_.AddJob(valueBoundaryEast, innerBoundary, outerBoundary);

        // We need to describe with inner & outer boundaries the isolation condition.
        // Initialize innerBoundary and outerBoundar with the proper offsets to describe a halo of 1 line
        //   at the north and south boundary.
        // Add the job with the isolationBoundary and the current inner & outer boundaries to the haloUpdate.
        innerBoundary.Init(0,0,0,0);
        outerBoundary.Init(0,0,-1,1);
        haloUpdate_.AddJob(isolationBoundary, innerBoundary, outerBoundary);

    }

    void Do() 
    {
        // Before applying the stencil, apply the halo exchange plus boundary conditions
        // using the halo update object.
        haloUpdate_.Apply();

        // Stencil computation
        stencil_.Apply();
    }

private:
    Stencil stencil_;

    HaloUpdateConfiguration haloConfiguration_;
    IJKRealHaloUpdate haloUpdate_;
    Real westValue_, eastValue_;
};

/*****************
 * MAIN FUNCTION *
 *****************/

TEST(Tutorial3, LaplaceTwoBoundaryEvolution)
{
    // Set number of time steps
    const int timeSteps = 20000;

    // Boundary condition values
    const Real westValue = 0.0, eastValue = 1;

    // instantiate a swap data field
    SwapDataField<IJKRealField> data;

    // initialize the sizes of the domain of our field
    IJKSize domain;
    domain.Init(100,60,1);
    data.in().Init("InitialField", domain, KBoundary());
    // initialize values of the field to 0.5
    for (int i = 0; i < data.in().calculationDomain().iSize(); ++i) {
        for (int j = 0; j < data.in().calculationDomain().jSize(); ++j) {
            data.in()(i, j, 0) = 0.5;
        }
    }

    // initialize also the output field
    data.out().Init("myfieldOut", data.in().calculationDomain(), KBoundary());

    // Initialize the movie
    Movie movie;
    movie.Init("Heat2D");

    // Generate the laplace object
    LaplaceStencil laplace;
    laplace.Init(data.in(), data.out(), westValue, eastValue);

    // Non-periodic boundary condition configuration, boundary condition applied
    // on the four sides of each plane

    // Initialize the reference
    LaplaceTwoBoundariesReference reference;
    reference.RegisterInitialCondition(data.in(), westValue, eastValue);


    // Perform time steps
    for(int t = 0; t <= timeSteps; ++t) 
    {

        // Apply the boundary condition and the Laplace stencil
        laplace.Do();

        // Store frame every 20 steps
        if(t%20 == 0)
        {
            std::cout << "Step "
                      << std::setw(5) << t << "/" << timeSteps << std::endl;
            movie.AddImage(data.out());
        }

        // Check reference
        ASSERT_TRUE(check(data.out(), reference, 0.00001));

        // Swap fields
        data.Swap();
    }

    movie.Finalize();

}
