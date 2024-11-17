#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <map>
#include <algorithm> 
#include <fstream>

#include "utils.h"

// g++ -std=c++20 -o potts potts.cpp

std::vector<std::vector<int>> lattice(int L, int q)
{
    std::vector<std::vector<int>> lat;
    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        for (int j = 0; j < L; j++)
        {
            float value = 1 + rand() % q;
            row.push_back(value);
        }
        lat.push_back(row);
    }
    return lat;
}

double Potts_Energy(const std::vector<std::vector<int>> &lattice)
{
    int L = lattice.size(); // Assuming lattice is L x L
    double E = 0.0;

    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            int s0 = lattice[i][j];
            // Check the four neighbors with periodic boundaries
            for (auto [dx, dy] : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1}, std::pair{0, 1}})
            {
                int x = (i + dx + L) % L; // Wrap around using modulo L
                int y = (j + dy + L) % L;
                int s1 = lattice[x][y];
                E -= (s1 == s0) ? 1 : 0;
            }
        }
    }

    return E / 2.0;
}

std::pair<float, float> Energy_Limit(int L, float J)
{
    float E_min = -2 * L * L * J;
    float E_max = 0;
    return std::make_pair(E_min, E_max);
}

bool isFlat(const std::map<int, int>& hist, double f, double& deviation)
{
    // Remove entries with zero counts
    std::vector<int> values;
    for (const auto& entry : hist)
    {
        if (entry.second > 0)
        {
            values.push_back(entry.second);
        }
    }

    if (values.empty())
    {
        deviation = 0;
        return false;
    }

    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double min_value = *std::min_element(values.begin(), values.end());

    deviation = min_value - 0.8 * mean;
    return true;
}

std::map<int, double> WangLandauPotts(std::vector<std::vector<int>> lat, int MC_N, int q)
{
    int L = lat.size();
    std::pair<float, float> E_limit = Energy_Limit(L, 1.0);
    float E_min = E_limit.first;
    float E_max = E_limit.second;
    float num_bins = -E_min + 1;

    std::map<int, double> g; // Density of states g(E)
    std::map<int, int> hist; // Histogram H(E)

    for (int E = E_min; E <= E_max; E += 1)
    { // Using increments of 2 for energy bins
        g[E] = 1.0;
        hist[E] = 0;
    }

    double f = std::exp(1.0); // Factor to multiply g(E) by
    while (f - 1 > 1e-4)
    {   
        #pragma omp parallel
        for (int i = 0; i < MC_N; i++)
        {
            int x = rand() % L;
            int y = rand() % L;
            int s0 = lat[x][y];

            float Old_E = Potts_Energy(lat);

            int s1 = 1 + rand() % q;
            lat[x][y] = s1;

            float New_E = Potts_Energy(lat);

            if (rand() / (double)RAND_MAX < std::min(1.0, g[Old_E] / g[New_E]))
            {
                g[New_E] *= f;
                hist[New_E] += 1;
            }
            else
            {
                lat[x][y] = s0;
                g[Old_E] *= f;
                hist[Old_E] += 1;
            }
            if(i % 100 == 0)
            {
                double deviation;
                if (isFlat(hist, f, deviation))
                {
                    f = std::sqrt(f);
                    hist.clear();
                    for (int E = E_min; E <= E_max; E += 1)
                    {
                        hist[E] = 0;
                    }
                }
            }
        }
    }
    return g;
}

int main()
{
    srand(time(NULL));

    int L = 15;
    int q = 2;
    std::vector<std::vector<int>> lat = lattice(L, q);
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            std::cout << lat[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << "\n"
              << Potts_Energy(lat) << std::endl;

    std::map<int, double> g = WangLandauPotts(lat, 2000000, 2);

    save_data(g, "g.txt");

    return 0;
}
