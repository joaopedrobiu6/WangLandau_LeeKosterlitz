import numpy as np
import matplotlib.pyplot as plt

def read_data(data):
    """Read the data from the file."""
    return np.loadtxt(data)

g = read_data('../pottscpp/lng.txt')

# shape if (n, 2)
# remove all the points where the first column is 1
g = g[g[:, 1] > 1]

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(g[:, 0], g[:, 1])
plt.xlabel("Energy")
plt.ylabel(r"$\Omega(E)$")
plt.title("Density of States vs Energy")
plt.subplot(1, 2, 2)
plt.plot(g[:, 0], g[:, 1])
plt.xlabel("Energy")
plt.ylabel(r"log($\Omega(E)$)")
plt.title("Log Density of States vs Energy")
plt.show()


