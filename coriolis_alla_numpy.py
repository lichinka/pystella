import numpy as np

from scipy.ndimage import convolve



# --------------------------------------------------------------
# DEFINITION of the Coriolis stencil object

coriolis = Stencil ( )


@coriolis.addStage
def uStageDo (utens, v, fc):
    """
    The 'Do' function of the U stage:

       utens   a 3D numpy array, representing ???;
       v       a 3D numpy array, representing ???;
       fc      a scalar representing the force.-
    """
    #
    # neighborhood for v.center and v.iplus1
    #
    neigh = np.array ([[0.0, 0.0, 0.0],
                       [0.0,  fc,  fc],
                       [0.0, 0.0, 0.0])
    #
    # the 'mode' and 'cval' parameters define the boundary conditions
    #
    res = (convolve (v, neigh, mode='constant', cval=0.0)) / 2.0

    #
    # neighborhood for v.jminus1 and v.jminus1.iplus1
    #
    neigh = np.array ([[0.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0],
                       [0.0,  fc,  fc])
    #
    # the 'mode' and 'cval' parameters define the boundary conditions
    #
    res += (convolve (v, neigh, mode='constant', cval=0.0)) / 2.0
    utens += res / 2.0


@coriolis.addStage
def vStageDo (vtens, u, fc):
    """
    The 'Do' function of the V stage:

       vtens   a 3D numpy array, representing ???;
       u       a 3D numpy array, representing ???;
       fc      a scalar, representing the force.-
    """
    #
    # neighborhood for u.jplus1 and u.center
    #
    neigh = np.array ([[0.0,  fc, 0.0],
                       [0.0,  fc, 0.0],
                       [0.0, 0.0, 0.0])
    #
    # the 'mode' and 'cval' parameters define the boundary conditions
    #
    res = (convolve (u, neigh, mode='constant', cval=0.0)) / 2.0

    #
    # neighborhood for u.iminus1 and u.iminus1.jplus1
    #
    neigh = np.array ([[ fc, 0.0, 0.0],
                       [ fc, 0.0, 0.0],
                       [0.0, 0.0, 0.0])
    #
    # the 'mode' and 'cval' parameters define the boundary conditions
    #
    res += ((convolve (u, neigh, mode='constant', cval=0.0)) / 2.0) * fc
    vtens += res / 2.0

#
# the output of the 'uSlowTensStage' is used as input of the 'vSlowTensStage'
#
#coriolis.addKLoop (sweep='kIncrement', (uStageDo,
#                                        vStageDo))

#
# these loops do not share any data during stencil execution
#
coriolis.addKLoop (sweep='kIncrement', (uStageDo))
coriolis.addKLoop (sweep='kIncrement', (vStageDo))



# --------------------------------------------------------------
# USAGE of the Coriolis stencil object defined above

#
# the calculation domain on which the stencil will be applied
#
calculationDomain = (8, 8, 2)

#
# boundaries in the K dimension: support is function-specific
#

#
# data-field definitions
#
u = np.random.random_integers (0, 10, calculationDomain)
v = np.random.random_integers (0, 10, calculationDomain)
utens = np.random.random_integers (0, 10, calculationDomain)
vtens = np.random.random_integers (0, 10, calculationDomain)

#
# put some values in the data fields: done along definition
#

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

