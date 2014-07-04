# -------------------------------------------------------------------
#
# Exposing Stella functionality to Python
#
# -------------------------------------------------------------------
import numpy as np

from stella.data import SwapField
from stella.stencil import LaplaceStencil



# set number of time steps to run the simulation
timeSteps = 20000

# boundary condition values
westValue = 0.0
eastValue = 1.0

# initialize the sizes of the domain of our fields
fieldDomain = np.array ((100, 60, 1))

# instantiate a swap field, containing input and output data fields
data = SwapField (name='MyDataField',
                  domain=fieldDomain)

# initialize values of the input field to 0.5
data.input.storage = np.array ([0.5] * in_data.domain.size)

#// Initialize the movie
#Movie movie;
#movie.Init("Heat2D");

# instantiate the Laplace stencil object, implemented in Stella
laplace = LaplaceStencil (data.input, data.output, westValue, eastValue)

# non-periodic boundary condition configuration, boundary condition applied
# on the four sides of each plane

# perform the simulation
for t in range (timeSteps):
    # apply the boundary condition and the Laplace stencil
    laplace.Do ( )

    # create a movie frame for every 20 steps
    if (t % 20) == 0:
        print ("Step %d/%d" % (t, timeSteps))
        #movie.AddImage(data.out());

    # swap the input and output data fields
    data.swap ( );

# create a movie from the intermediate frames
#movie.finalize ( );

