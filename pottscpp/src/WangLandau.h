#include "iostream"
#include "vector"
#include "math.h"
#include "random"
#include "map"
#include "algorithm"
#include "fstream"
#include "Lattice.h"

bool isFlat(const std::map<int, int> &hist, double h_tol);

std::map<int, double> WangLandauPotts(PottsLattice lat, int MC_N, int q, double f_tol, double h_tol, bool NoLog = false);

// std::map<int, double> WangLandauIsing(IsingLattice lat, int MC_N);