import random



# --------------------------------------------------------------
# DEFINITION of the Coriolis stencil object

coriolis = Stencil ( )

#
# add the U stage to the Coriolis stencil
#
uSlowTensStage = coriolis.addStage ( )

@uSlowTensStage.attachDo 
def uStageDo (utens, v, fc):
    """
    The 'Do' function of the U stage, with the Coriolis force directly applied:

       utens   a 3D numpy array, representing ???;
       v       a 3D numpy array, representing ???;
       fc      a scalar representing the force.-
    """
    res = fc * average (v, v.iplus1)
    res += fc * average (v.jminus1, v.jminus1.iplus1)
    utens += res / 2.0

#
# add the V stage to the Coriolis stencil
#
vSlowTensStage = coriolis.addStage ( )

@vSlowTensStage.attachDo
def vStageDo (vtens, u, fc):
    """
    The 'Do' function of the V stage, with the Coriolis force defined
    as a private function:

       vtens   a 3D numpy array, representing ???;
       u       a 3D numpy array, representing ???;
       fc      a scalar, representing the force.-
    """
    def coriolisForce (frc, vel):
        """
        Calculates the Coriolis force:

            fc   constant Coriolis force factor;
            vel  velocity used to calculte the force.

        """
        return frc * vel

    res = coriolisForce (fc, average (u.jplus1, u))
    res += coriolisForce (fc, average (u.iminus1, u.iminus1.jplus1))
    vtens += res / 2.0

#
# the output of the 'uSlowTensStage' is used as input of the 'vSlowTensStage'
#
#coriolis.addKLoop (sweep='kIncrement', (uSlowTensStage,
#                                        vSlowTensStage))

#
# these loops do not share any data whithin the stencil execution
#
coriolis.addKLoop (sweep='kIncrement', (uSlowTensStage))
coriolis.addKLoop (sweep='kIncrement', (vSlowTensStage))


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
    #   in the `Do' functions at stage level, otherwise a
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

