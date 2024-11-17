#include "Lattice.h"

// implementation of the PottsLattice class

PottsLattice::PottsLattice(int L, int q, float J) : L(L), q(q), J(J)
{
    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        for (int j = 0; j < L; j++)
        {
            float value = 1 + rand() % q;
            row.push_back(value);
        }
        lattice.push_back(row);
    }
}

double PottsLattice::Potts_Energy()
{
    int L = lattice.size();
    double E = 0.0;


    // maybe we can iterate over only half of the lattice to avoid double counting
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            int s0 = lattice[i][j];
            for (auto [dx, dy] : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1}, std::pair{0, 1}})
            {
                int x = (i + dx + L) % L;
                int y = (j + dy + L) % L;
                int s1 = lattice[x][y];
                E -= (s1 == s0) ? 1 : 0;
            }
        }
    }
    return J * E * 0.5;
}

std::pair<float, float> PottsLattice::Energy_Limit()
{
    float E_min = -2 * L * L * J;
    float E_max = 2* L * L * J;
    return std::make_pair(E_min, E_max);
}

void PottsLattice::PrintLattice()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            std::cout << lattice[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}

// implementation of the IsingLattice class

IsingLattice::IsingLattice(int L, float J) : L(L), J(J)
{

    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        for (int j = 0; j < L; j++)
        {
            float value = -1 + rand() % 2;
            if (value == 0)
            {
                value = 1;
            }
            row.push_back(value);
        }
        lattice.push_back(row);
    }
}

double IsingLattice::Ising_Energy()
{
    int L = lattice.size();
    double E = 0.0;

    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            int s0 = lattice[i][j];
            int WF = lattice[(i + 1) % L][j] + lattice[i][(j + 1) % L] + lattice[(i - 1 + L) % L][j] + lattice[i][(j - 1 + L) % L];

            E -= J * WF * s0;
        }
    }
    return E * 0.5;
}

std::pair<float, float> IsingLattice::Energy_Limit()
{
    float E_min = -2 * L * L * J;
    float E_max = 2 * L * L * J;
    return std::make_pair(E_min, E_max);
}

void IsingLattice::PrintLattice()
{
    std::cout << "Printing lattice" << std::endl;
    std::cout << "L = " << L << std::endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            std::cout << lattice[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}