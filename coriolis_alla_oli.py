from stencil_kernel import *




# --------------------------------------------------------------
# DEFINITION of the Coriolis stencil object
#
class CoriolisKernel (StencilKernel):
    """
    Class definition of the Coriolis stencil.-
    """
    #
    # neighborhood definitions for the UStage
    #
    _neighborhoods    = {}
    _neighborhoods['first'] = [(0, 0), (1, 0)]
    _neighborhoods['1'] = [(0, -1), (1, -1)]

    def _USlowTensStage (self, ctr, in_v, in_fc):
        """
        The 'Do' function of the U stage.
        This stage uses neighborhood definitions.-
        """
        return ( in_fc * np.average (in_v.neigh (ctr, 'first')) +
                 in_fc * np.average (in_v.neigh (ctr, '1'))
               ) / 2.0

    def _VSlowTensStage (self, ctr, in_u, in_fc):
        """
        The 'Do' function of the V stage.
        This stage uses enumerated neighbors.-
        """
        return ( in_fc * np.average (in_v[ctr[0, 0]], in_v[ctr[0, 1]]) +
                 in_fc * np.average (in_v[ctr[-1, 0]], in_v[ctr[-1, 1]])
               ) / 2.0

    def kernel (self, in_u, in_v, in_fc, out_utens, out_vtens):
        """
        This stencil comprises two independent stages.-
        """
        for p in out_utens.interior_points (sweep='cKIncrement'):
            out_utens[p] += self._USlowTensStage (p, in_v, in_fc) 

        for p in out_vtens.interior_points (sweep='cKIncrement'):
            out_vtens[p] -= self._VSlowTensStage (p, in_u, in_fc)



# --------------------------------------------------------------
# USAGE of the Coriolis stencil object defined above
#
kernel = CoriolisKernel ( )

kernel.compilation.should_unroll = False
kernel.compilation.backend       = 'cxx'

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
u = IJKRealField (np.random.random (0, 10, calculationDomain))
u.set_boundaries (**boundaries)

v = IJKRealField (np.random.random (0, 10, calculationDomain))
v.set_boundaries (**boundaries)

utens = IJKRealField (np.random.random (0, 10, calculationDomain))
utens.set_boundaries (**boundaries)

vtens = IJKRealField (np.random.random (0, 10, calculationDomain))
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
# print out the final state
#
print ("State after initial step")
print (utens)
print (vtens)

