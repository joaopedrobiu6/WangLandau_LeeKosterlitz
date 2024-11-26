#include <iostream>
#include <vector>
#include <set>
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()

// Function to calculate energy
double energy(const std::vector<std::vector<int>>& lattice, int L, double J) {
    double E = 0.0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int s0 = lattice[i][j];
            for (auto [dx, dy] : {std::pair{1, 0}, std::pair{0, 1}}) {
                int x = (i + dx) % L;
                int y = (j + dy) % L;
                int s1 = lattice[x][y];
                E -= (s0 == s1) ? 1 : 0;
            }
        }
    }
    return J * E;
}

// Function to generate a random lattice
std::vector<std::vector<int>> generateRandomLattice(int L, int n) {
    std::vector<std::vector<int>> lattice(L, std::vector<int>(L));
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            lattice[i][j] = rand() % n + 1; // Random value between 1 and n
        }
    }
    return lattice;
}

int main() {
    int L = 2;        // Lattice size
    int n = 8;         // Number of possible states (e.g., values 1 to 8)
    double J = 1.0;    // Coupling constant
    int samples = 1000000; // Number of random samples

    std::set<double> energyLevels;
    srand(static_cast<unsigned>(time(0))); // Seed random number generator

    // Randomly sample matrices and compute energies
    for (int i = 0; i < samples; ++i) {
        std::vector<std::vector<int>> lattice = generateRandomLattice(L, n);
        double E = energy(lattice, L, J);
        energyLevels.insert(E);
    }

    // Output the distinct energy levels
    std::cout << "Found " << energyLevels.size() << " distinct energy levels.\n";
    for (double E : energyLevels) {
        std::cout << E << "\n";
    }

    return 0;
}
