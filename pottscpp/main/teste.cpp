#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 17;
    float J = 1.0;

    IsingLattice ising(L, J);
    ising.PrintLattice();
    std::cout << "Energy = " << ising.Ising_Energy() << std::endl;
    
    auto g = WangLandauIsing(ising, 100000);
    save_data(g, "g.txt");
    return 0;
}