#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 3;
    float J = 1.0;
    int q = 8;

    PottsLattice potts(L, q, J);
    potts.PrintLattice();
    std::cout << "Energy = " << potts.Potts_Energy() << std::endl;

    std::map<int, double> g = WangLandauPotts(potts, 1000, q);
    save_data(g, "g.txt");
    return 0;

//     IsingLattice ising(L, J);
//     ising.PrintLattice();
//     std::cout << "Energy = " << ising.Ising_Energy() << std::endl;
    
//     auto g = WangLandauIsing(ising, 10000000);
//     save_data(g, "g.txt");
//     return 0;
}