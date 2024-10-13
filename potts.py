import numpy as np
import matplotlib.pyplot as plt

class Potts:
    """Class for simulating a Potts model and Wang-Landau method for DOS."""

    def __init__(self, T=1.0, L=2, q=2):
        self.T = T        # Temperature
        self.L = L        # Linear size of lattice
        self.q = q        # Number of states
        self.reset()

    def mcmc(self,n):
        """Metropolis Monte Carlo moves"""
        self.metropolis(n)

    def metropolis(self, n)



