#include "iostream"
#include "vector"
#include "math.h"
#include "random"
#include "map"
#include "algorithm" 
#include "fstream"
#include "Lattice.h"


bool isFlat(const std::map<int, int>& hist, double f, double& deviation);

std::map<int, double> WangLandauPotts(std::vector<std::vector<int>> lat, int MC_N, int q);

std::map<int, double> WangLandauIsing(IsingLattice lat, int MC_N);