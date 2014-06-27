from stencil_kernel import *




# --------------------------------------------------------------
# DEFINITION of the FastWavesSCTridiag stencil object
#
class FastWavesSCTridiag (StencilKernel):
    """
    Class definition of the FastWavesSCTridiag stencil.-
    """
    def __init__ (self, domain):
        #
        # call the parent's constructor
        #
        super.__init__ (self)

        #
        # output fields
        #
        self.y   = IJKRealField (domain)
        self.tmp = IJKRealField (domain)
        self.bet = IJKRealField (domain)


    def _ForwardStage (self, ctr, in_b, in_rhs):
        """
        The 'Do' function of the Forward stage, applied over the KMinimunCenter.-
        """
        self.bet[ctr] = in_b[ctr]
        self.y[ctr]   = in_rhs[ctr] / self.bet[ctr]


    def _ForwardStageFull (self, ctr, in_a, in_b, in_c, in_rhs):
        """
        The 'Do' function of the Forward stage, applied over the FullDomain.-
        """
        self.tmp[ctr] = in_c[ctr[0, 0, -1]] / self.bet[ctr]
        self.bet[ctr] = in_b[ctr] - in_a[ctr] * self.tmp[ctr]
        self.y[ctr]   = (in_rhs[ctr] - 
                         in_a[ctr] * self.y[ctr[0, 0, -1]] / self.bet[ctr]


    def _BackwardStage (self, ctr, in_tmp):
        """
        The 'Do' function of the Backward stage.-
        """
        self.y[ctr] -= in_tmp[ctr[0, 0, 1]] * self.y[ctr[0, 0, 1]]


    def kernel (self, in_a, in_b, in_c, in_rhs):
        """
        This stencil comprises two independent stages.-
        """
        for p in self.y.interior_points (sweep='cKIncrement',
                                         height='KMinimumCenter'):
            self._ForwardStage (p, in_b, in_rhs) 
                                         
        for p in self.y.interior_points (sweep='cKIncrement'):
            self._ForwardStageFull (p, in_a, in_b, in_c, in_rhs) 
        #
        # using `self.tmp' as an input field creates a dependency among the loops
        #
        for p in self.y.interior_points (sweep='cKDecrement'):
            self._BackwardStage (p, self.tmp)



# --------------------------------------------------------------
# USAGE of the FastWavesSCTridiag object defined above
#

#
 the calculation domain on which the stencil will be applied
#
calculationDomain = (8, 8, 2)

kernel = FastWavesSCTridiag (calculationDomain)

kernel.compilation.should_unroll = False
kernel.compilation.backend       = None     # runs in pure Python (default)

#
# boundaries definition in all three dimensions
#
boundaries = {'i': (0, 0),
              'j': (0, 0),
              'k': (0, 1)}
#
# data-field definitions
#
a = IJKRealField (np.random.random (0, 10, calculationDomain))
a.set_boundaries (**boundaries)

b = IJKRealField (np.random.random (0, 10, calculationDomain))
b.set_boundaries (**boundaries)

c = IJKRealField (np.random.random (0, 10, calculationDomain))
c.set_boundaries (**boundaries)

rhs = IJKRealField (np.random.random (0, 10, calculationDomain))
rhs.set_boundaries (**boundaries)

#
# boundaries for the output fields
#
kernel.y.set_boundaries   (**boundaries)
kernel.bet.set_boundaries (**boundaries)
kernel.tmp.set_boundaries (**boundaries)

#
# print out the initial state
#
print ("Initial state")
print (kernel.y)

#
# apply the stencil
#
kernel.kernel (a,
               b,
               c,
               rhs)
#
# print out the final state
#
print ("State after initial step")
print (self.y)

