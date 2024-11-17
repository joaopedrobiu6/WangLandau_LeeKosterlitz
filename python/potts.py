import numpy as np
import matplotlib.pyplot as plt
class Potts:
    """Class for simulating a Potts model and Wang-Landau method for DOS."""

    def __init__(self, T, L, q, J=1):
        self.T = T        # Temperature
        self.L = L        # Linear size of lattice
        self.q = q+1       # Number of states
        self.lattice = np.random.randint(1, self.q, (self.L, self.L))  # Random initial state
        self.E = self.energy(self.lattice)  # Initial energy
        self.J = J

    def energy(self, lattice):
        """Compute the energy of the lattice."""
        E = 0
        for i in range(self.L):
            for j in range(self.L):
                s0 = self.lattice[i, j]
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    x = (i + dx) % self.L
                    y = (j + dy) % self.L
                    s1 = self.lattice[x, y]
                    E -= 1 if s1 == s0 else 0
        return (E/2.)
    
    def energy_limits(self, L, J):
        """Compute the energy limits of the lattice."""
        print(f"E_min: {-J*2*L**2} \t E_max: 0")
        return -J*2*self.L**2, 0
    
    def deltaE(self, x, y, s0):
        dE = 0
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            x1 = (x + dx) % self.L
            y1 = (y + dy) % self.L
            s1 = self.s[x1, y1]
            print(f"s0: {s0} \t s1: {s1}",)

            dE -= 1 if s1 == s0 else 0
        
        return dE
        
    def isFlat(self, hist, f):
        """Check if the histogram is flat."""
        # remove entries with zero counts
        hist_new = {k: v for k, v in hist.items() if v > 0}
        values = list(hist_new.values())
        if np.min(values) - 0.8*np.mean(values) > 0: # and np.max(values) < 1.2*np.mean(values):
            return True, np.min(values) - 0.8*np.mean(values)
        else:
            return False, np.min(values) - 0.8*np.mean(values)
    
    def WangLandau(self, N_sweeps):
        """Wang-Landau algorithm to estimate the density of states."""
        print("Wang-Landau")
        print(f"Lattice:\n{self.lattice}")
        
        # Initialize the density of states and histogram
        E_min, E_max = self.energy_limits(self.L, self.J)
        energy_bins = np.linspace(E_min, E_max, -E_min+1)

        g = {E: 1.0 for E in energy_bins}
        hist = {E: 0 for E in energy_bins}

        f = np.e
        while f-1 > 1e-4:
            for ii in range(N_sweeps):
                current_energy = self.energy(self.lattice)
                # choose a site at random and propose a new state
                x, y = np.random.randint(0, self.L), np.random.randint(0, self.L)
                s0 = self.lattice[x, y]
                
                s1 = np.random.randint(1, self.q)
                self.lattice[x, y] = s1
                new_energy = self.energy(self.lattice)

                if (np.random.rand() < min(1, g[int(current_energy)]/g[int(new_energy)])):
                    g[int(new_energy)] *= f
                    hist[int(new_energy)] += 1
                else:
                    self.lattice[x, y] = s0
                    g[int(current_energy)] *= f
                    hist[int(current_energy)] += 1
                if ii % 50 == 0:
                    isFlat, flatness = self.isFlat(hist, f)
                 
                    if isFlat:
                        f = np.sqrt(f)
                        hist = {E: 0 for E in energy_bins}
                        
        g_x = np.array(list(g.keys()))
        g_y = np.array(list(g.values()))
        return g_x, g_y

if __name__ == "__main__":
    P = Potts(T=1, L=3, q=8)
    g_x, g_y = P.WangLandau(1000)

    plt.plot(g_x, g_y/np.max(g_y), label="Density of states")
    plt.legend()
    plt.show()
