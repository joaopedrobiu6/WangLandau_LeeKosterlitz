#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 12;
    float J = 1.0;
    int q = 5;
    double f_tol = 1.0e-8, h_tol = 0.80;
    bool NoLog = false;
    int sampleInterval = 10000, f_factor = 2;
    std::string filename = GetFilename(L, J, q);

    PottsLattice potts(L, q, J);
    potts.PrintLattice();
    std::cout << "Energy = " << potts.Potts_Energy() << std::endl;

    std::map<int, double> lng = WangLandauPotts(potts, 10000000, q, f_tol, h_tol, NoLog, sampleInterval, f_factor, filename);
    return 0;
}