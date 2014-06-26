from stencil_kernel import *
import sys
import numpy
import math

width = int(sys.argv[2])
height = int(sys.argv[3])
image_in = open(sys.argv[1], 'rb')
stdev_d = 3
stdev_s = 70
radius = stdev_d * 3



# --------------------------------------------------------------
# DEFINITION of the Coriolis stencil object

class CoriolisKernel (StencilKernel):
    """
    Class definition of the Coriolis stencil.-
    """
    #
    # neighborhood definitions for the UStage
    #
    _neighborhoods    = []
    _neighborhoods[0] = [(0, 0), (1, 0)]
    _neighborhoods[1] = [(0, -1), (1, -1)]

    def _USlowTensStage (self, center, in_v, in_fc):
        """
        The 'Do' function of the U stage.
        This stage uses neighborhood definitions.-
        """
        return ( in_fc * np.average (in_v.neigh (center, 0)) +
                 in_fc * np.average (in_v.neigh (center, 1))
               ) / 2.0


    def _VSlowTensStage (self, center, in_u, in_fc):
        """
        The 'Do' function of the V stage.
        This stage uses enumerated neighbors.-
        """
        return ( in_fc * np.average (in_v.neigh (center, ((0, 0),
                                                          (0, 1)))) +
                 in_fc * np.average (in_v.neigh (center, ((-1, 0),
                                                          (-1, 1))))
               ) / 2.0


    def kernel (self, in_u, in_v, in_fc, out_utens, out_vtens):
        """
        This stencil comprises two independent stages.-
        """
        out_utens.define_sweep ('cKIncrement')
        for p in out_utens.interior_points ( ):
            out_utens[p] += self._USlowTensStage (p, in_v, in_fc) 

        out_vtens.define_sweep ('cKIncrement')
        for p in out_vtens.interior_points ( ):
            out_vtens[p] -= self._VSlowTensStage (p, in_u, in_fc)


# --------------------------------------------------------------
# USAGE of the Coriolis stencil object defined above

kernel               = CoriolisKernel ( )
kernel.should_unroll = False
kernel.pure_python   = False

#
# the calculation domain on which the stencil will be applied
#
calculationDomain = (8, 8, 2)

#
# boundaries definition in all three dimensions
#
boundaries = {'i': (0, 0),
              'j': (0, 0),
              'k': (0, 0)}
#
# data-field definitions
#
u = IJKRealField (np.random.random_floats (0, 10, calculationDomain))
u.set_boundaries (**boundaries)

v = IJKRealField (np.random.random_floats (0, 10, calculationDomain))
v.set_boundaries (**boundaries)

utens = IJKRealField (np.random.random_floats (0, 10, calculationDomain))
utens.set_boundaries (**boundaries)

vtens = IJKRealField (np.random.random_floats (0, 10, calculationDomain))
vtens.set_boundaries (**boundaries)

#
# print out the initial state
#
print ("Initial state")
print (utens)
print (vtens)

#
# apply the stencil
#
kernel.kernel (u,
               v,
               3.5,
               utens,
               vtens)
#
# print the final state
#
print ("State after initial step")
print (utens)
print (vtens)

