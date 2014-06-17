import numpy as np



# --------------------------------------------------------------
# DEFINITION of the Coriolis stencil object

coriolis = Stencil ( )

#
# define the U stencil stage in DSL
#
coriolis.addStage ("uStageDo", """
    Kernel_2D_Begin (Field_2D_out utens, Field_2D_in v, float fc)
        utens[0,0] += ( fc * ((v[0,0] + v[1,0]) / 2.0) +
                        fc * ((v[0,-1] + v[1,-1]) / 2.0)
                      ) / 2.0;
    Kernel_2D_End
""")

#
# define the V stencil stage in DSL
#
coriolis.addStage ("vStageDo", """
    Kernel_2D_Begin (Field_2D_out vtens, Field_2D_in u, float fc)
        vtens[0,0] += ( fc * ((u[0,1] + u[0,0]) / 2.0) +
                        fc * ((u[-1,0] + u[-1,1]) / 2.0)
                      ) / 2.0;
    Kernel_2D_End
""")


#
# the output of the U stage is used as input of the V stage
#
#coriolis.addLoop ("""Concat_Sweep_K_Increment_Begin
#                        uStageDo;
#                        vStageDo;
#                     Concat_Sweep_K_Increment_End
#                   """)

#
# these loops do not share any data during stencil execution
#
coriolis.addLoop ("""Sweep_K_Increment_Begin
                        uStageDo;
                        vStageDo;
                     Sweep_K_Increment_End
                  """)


# --------------------------------------------------------------
# USAGE of the Coriolis stencil object defined above

#
# the calculation domain on which the stencil will be applied
#
calculationDomain = IJKSize (8, 8, 2)

#
# no boundaries in the K dimension
#
kBoundary = KBoundary (0, 0)

#
# data-field definitions
#
u = IJKRealField (calculationDomain, kBoundary)
v = IJKRealField (calculationDomain, kBoundary)
utens = IJKRealField (calculationDomain, kBoundary)
vtens = IJKRealField (calculationDomain, kBoundary)

#
# put some values in the data fields
#
map (random.random ( ), u)
map (random.random ( ), v)
map (random.random ( ) * 10.0, utens)
map (random.random ( ) * 10.0, vtens)

#
# print out the initial state
#
print ("Initial state")
print (utens)
print (vtens)

#
# apply the stencil in 3 time steps
#
for step in xrange (3):
    #
    # apply the stencil with the field and scalar variables:
    #   the names (i.e., dictionary keys) should match those used 
    #   in the `stage' functions at stage level, otherwise a
    #   StencilCompilationException would occur at runtime
    #
    coriolis.apply (fields={'u': u,
                            'v': v,
                            'utens': utens,
                            'vtens': vtens},
                    scalars={'fc': 3.5})

    #
    # print the situation after each step
    #
    print ("State after time step", step)
    print (utens)
    print (vtens)

