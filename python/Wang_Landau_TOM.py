import numpy as np
from numpy.random import rand as rnd
import time
import matplotlib.pyplot as plt
import json

def calEnergy(lattice,J):
    # Energy of a 2D Ising lattice
    E_N = 0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            WF = lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L] + lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L] + lattice[(i - 1) % L,(j - 1) % L] + lattice[(i + 1) % L,(j + 1) % L]
            E_N -= J * WF * S  # Each neighbor gives environment energy
    return int(E_N / 2.)  # Counted twice

def Energy_list(L,J) :
    energy_list = [-3 * L ** 2, -3 * L ** 2 + 12, -3 * L ** 2 + 20]

    # Incrémenter de 4 jusqu'à L**2
    current_value = energy_list[-1] + 4
    while current_value <= L ** 2:
        energy_list.append(current_value)
        current_value += 4

    energy_list = [J * energy for energy in energy_list]
    print(energy_list)
    return energy_list

# Définition des constantes
J = 1
L = 16
N = L * L
MCMC_Sweep = 20000000
histo_flat = 0.1

# Initialisation de la lattice
lattice = np.random.choice([-1, 1], size=(L, L))  # Crée une lattice au hasard
E_system = calEnergy(lattice,J)  # Calcule l'énergie initiale du système
print(f'Start at {E_system}')
time.sleep(2)

# Définir la plage des énergies visitables
Energies = Energy_list(L,J)

lng = {E: 0 for E in Energies}  # Densité d'états initialisée à 1.0
H = {E: 0 for E in Energies}  # Histogramme
f = np.e  # Facteur de modification
lnf = np.log(f)

it = 0
for it in range(MCMC_Sweep) :
    n = int(np.random.rand() * N)  # The site to flip
    (i, j) = (int(n % L), int(n / L))  # The coordinates of the site
    S = lattice[i, j]  # its spin
    WF = lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L] + lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L] + lattice[(i - 1) % L,(j - 1) % L] + lattice[(i + 1) % L,(j + 1) % L]
    E_new = E_system + 2 * J * S * WF

    if E_new not in Energies :
        print("Erreur")

    lnP = lng[E_system] - lng[E_new]# Utilisation de log pour la probabilité
    if lnP > np.log(rnd()):
            lattice[i, j] = -S  # Accepter le changement
            E_system = E_new  # Mettre à jour l'énergie du système

    # Mettre à jour l'histogramme et la densité d'états
    H[E_system] += 1  # Mettre à jour l'histogramme
    lng[E_system] += lnf  # Mettre à jour la densité d'états


    # Vérifier la condition de flatness tous les 100 itérations
    if it % 100 == 0:
        mean_H = np.mean(list(H.values()))
        min_H = min(H.values())
        #print(f"Minimum H: {min_H}, Average H: {mean_H * histo_flat}")

        if min_H > mean_H * histo_flat:  # Condition de flatness
            # Normaliser l'histogramme
            H = {E: 0 for E in Energies}  # Réinitialiser l'histogramme
            lnf /= 8  # Réduire le facteur de modification
            print(f"Step in {it} iterations, f = {np.exp(lnf)}")

print(lng)

# Sauvegarder le dictionnaire lng dans un fichier JSON
with open("lng_resultsF.json", "w") as outfile:
    json.dump(lng, outfile)
    print("Densité d'états sauvegardée dans lng_results.json")


# Tracer la densité d'états
plt.plot(Energies, list(lng.values()), '-o', label='lng(E)')
plt.xlabel('Energy')
plt.ylabel('Density of States lng(E)')
plt.legend(loc='best')
plt.title('Densité d\'états de l\'Ising 2D')
plt.show()

#https://www.wellesu.com/10.1007/s11467-007-0026-3