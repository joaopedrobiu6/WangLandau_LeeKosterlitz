#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 12;
    float J = 1.0;
    int q = 10;
    double f_tol = 1.0e-6, h_tol = 0.80;
    bool NoLog = false;

    PottsLattice potts(L, q, J);
    potts.PrintLattice();
    std::cout << "Energy = " << potts.Potts_Energy() << std::endl;

    std::map<int, double> lng = WangLandauPotts(potts, 1000, q, f_tol, h_tol, NoLog);
    save_data(lng, "lng.txt");
    return 0;
}