"""Generic functions"""

import numpy as np
from numpy import zeros

@np.vectorize
def delta_lorentzian(x, eta):
    """The lorentzian representation of a delta function."""
    return (eta / np.pi) / (x ** 2 + eta ** 2)

