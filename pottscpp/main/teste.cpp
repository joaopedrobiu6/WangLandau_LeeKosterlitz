#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 20;
    float J = 1.0;
    int q = 8;
    double f_tol = 1.0e-6, h_tol = 0.80;
    bool NoLog = false;
    int sampleInterval = 5000, f_factor = 5;
    std::string filename = "L20_J1_q8.txt";

    PottsLattice potts(L, q, J);
    potts.PrintLattice();
    std::cout << "Energy = " << potts.Potts_Energy() << std::endl;

    std::map<int, double> lng = WangLandauPotts(potts, 10000, q, f_tol, h_tol, NoLog, sampleInterval, f_factor, filename);
    save_data(lng, "lng.txt");
    return 0;
}