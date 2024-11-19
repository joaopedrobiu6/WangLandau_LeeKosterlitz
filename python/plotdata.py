import numpy as np
import matplotlib.pyplot as plt

def read_data(data):
    """Read the data from the file."""
    return np.loadtxt(data)

g = read_data('../pottscpp/g.txt')

# shape if (n, 2)
# remove all the points where the first column is 1
g = g[g[:, 1] > 1]
plt.plot(g[:, 0], g[:, 1])
plt.show()


