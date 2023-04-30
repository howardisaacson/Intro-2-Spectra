
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt 
  
## generate the data and plot it for an ideal normal curve
  
import numpy as np
import matplotlib.pyplot as plt

def gaussian_curve(width, center, amplitude):
    """Generate a Gaussian curve with given width, center location, and amplitude."""
    return amplitude * np.exp(-((center) / width) ** 2)

def repeating_gaussian_curves(x, num_repeats, width, amplitude):
    """Generate a repeating set of Gaussian curves with given number of repeats, width, center, and amplitude."""
    x = np.asarray(x)
    for i in range(int(x[0]), int(x[-1]), int(len(x) / num_repeats)):
        y = gaussian_curve(width, i, amplitude)
        plt.plot(x, y)


# Example usage
x = np.linspace(-10, 10, 1000)
y = repeating_gaussian_curves(x, num_repeats=3, width=2, amplitude=1)
plt.show()

    


