import numpy as np
import random
import math
import os

# Constants
L = 10  # Lattice size (LxL)
q = 10  # Number of spin states
mc_steps = 100000  # Number of Monte Carlo steps
n_bins = 100  # Number of bins for the histogram
ener_min = -100  # Minimum energy
ener_max = 100  # Maximum energy
bin_width = (ener_max - ener_min) / n_bins  # Bin width
abscisse = np.linspace(-2, 2, n_bins)  # Abscissa values for histogram
abscissem = np.linspace(-1, 1, n_bins)  # Abscissa values for magnetization histogram

# Initializing histograms
histo = np.zeros(n_bins)
g_e = np.zeros(n_bins)
g_em = np.zeros(n_bins)
histom = np.zeros(n_bins)

# Generate qlist for spin states
qlist = np.arange(1, q + 1)

# Helper functions
def energy(spin, L, q):
    ener = 0
    for i in range(L):
        for j in range(L):
            ip = (i + 1) % L
            im = (i - 1 + L) % L
            jp = (j + 1) % L
            jm = (j - 1 + L) % L
            if spin[ip][j] == spin[i][j]:
                ener -= 1
            if spin[im][j] == spin[i][j]:
                ener -= 1
            if spin[i][jp] == spin[i][j]:
                ener -= 1
            if spin[i][jm] == spin[i][j]:
                ener -= 1
    return ener / 2.0  # Normalize by dividing by 2

def invert(spin, L, q, qlist, s, t, b):
    spin[s][t] = qlist[b]
    return energy(spin, L, q)

def magnetization(spin, L, q):
    summ = 0.0
    for i in range(L):
        for j in range(L):
            p = 1 if spin[i][j] == 1 else 0
            summ += q * p - 1
    return summ / (L * L * (q - 1))

def histogram_check(histo, n_bins, mcs):
    summ = 0.0
    cc1 = 0
    for u in range(n_bins):
        if histo[u] > 0:
            cc1 += 1
            summ += histo[u]
    optimal = summ / float(cc1)
    proc_optimal = optimal / 5.0

    count1 = 0
    for u in range(n_bins):
        if histo[u] > 0:
            if (optimal - proc_optimal > histo[u]) or (histo[u] - optimal < proc_optimal):
                count1 += 1

    return 1 if cc1 == count1 else 10

def ener_in_range(ener_min, ener_max, spin, L, q, qlist):
    e1 = energy(spin, L, q)
    while e1 < ener_min or e1 > ener_max:
        for i in range(L):
            for j in range(L):
                spin[i][j] = random.choice(qlist)
        e1 = energy(spin, L, q)

# Initialize lattice
spin = np.random.choice(qlist, size=(L, L))

# Main simulation loop
f = 1.0  # Weighting factor

# Create output directories
os.makedirs("data", exist_ok=True)

# Start the simulation
print(f"LATTICE {L}x{L}")
print(f"{ener_min} < range < {ener_max}")

# Initialize histogram and energy range
ener_in_range(ener_min, ener_max, spin, L, q, qlist)

# Weighting factor loop
for i in range(10):  # nb_f = 10
    print(f"Weighting factor loop: {i + 1}")
    
    # Monte Carlo loop
    for mcs in range(1, mc_steps + 1):
        for _ in range(L * L):  # One sweep (L * L iterations)
            e1 = energy(spin, L, q)
            m1 = magnetization(spin, L, q)
            
            # Select a random spin and assign it a random state from qlist
            s = random.randint(0, L - 1)
            t = random.randint(0, L - 1)
            b = random.randint(0, q - 1)  # Random spin state
            
            e2 = invert(spin, L, q, qlist, s, t, b)
            
            if ener_min <= e2 <= ener_max:
                k = min(n_bins - 1, max(0, int((2.0 + e1) / bin_width)))
                r = min(n_bins - 1, max(0, int((2.0 + e2) / bin_width)))
                g = min(n_bins - 1, max(0, int((1.0 - m1) / (bin_width / 2.0))))
                h = min(n_bins - 1, max(0, int((1.0 - m1) / (bin_width / 2.0))))

                prob = 1.0 if g_e[r] <= g_e[k] else math.exp(g_e[k] - g_e[r])
                
                if random.random() <= prob:
                    histo[r] += 1
                    g_e[r] += math.log(f)
                    histom[h] += 1
                    g_em[h] += m1
                else:
                    histo[k] += 1
                    g_e[k] += math.log(f)
                    histom[g] += 1
                    g_em[g] += m1
        
        # Check histogram flatness every 10,000 steps
        if mcs % 10000 == 0:
            print(f"Monte Carlo step: {mcs}")
            u = histogram_check(histo, n_bins, mcs)
            if u < 5:
                print(f"Histogram is flat at step: {mcs}")
                break  # Exit if the histogram is flat enough
    
    # Update weighting factor
    f = math.sqrt(f)
    if i < 9:  # Reset histogram for the next iteration
        histo.fill(0)

# Writing the results to files
def write_to_file(filename, data, header=""):
    with open(filename, 'w') as outFile:
        outFile.write(header)
        for i in range(n_bins):
            if data[i] > 0:
                outFile.write(f"{abscisse[i]} {data[i]}\n")

write_to_file("data/histogramq10.dat", histo)
write_to_file("data/density_of_stateq10.dat", g_e, header="Density of States\n")
write_to_file("data/DOSmagq10.dat", g_em, header="Magnetization\n")

print("Simulation finished!")
