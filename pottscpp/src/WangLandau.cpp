#include "WangLandau.h"

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

std::map<int, double> WangLandauPotts(PottsLattice lat, int MC_N, int q)
{
    int L = lat.lattice.size();
    std::pair<float, float> E_limit = lat.Energy_Limit();
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
    while (f - 1 > 1e-8)
    {   
        #pragma omp parallel
        for (int i = 0; i < MC_N; i++)
        {
            std::cout << "f: " << f << std::endl;
            int x = rand() % L;
            int y = rand() % L;
            int s0 = lat.lattice[x][y];

            float Old_E = lat.Potts_Energy();

            int s1 = 1 + rand() % q;
            lat.lattice[x][y] = s1;

            float New_E = lat.Potts_Energy();

            if (rand() / (double)RAND_MAX < std::min(1.0, g[Old_E] / g[New_E]))
            {
                g[New_E] *= f;
                hist[New_E] += 1;
            }
            else
            {
                lat.lattice[x][y] = s0;
                g[Old_E] *= f;
                hist[Old_E] += 1;
            }
            if(i % 10 == 0)
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

std::map<int, double> WangLandauIsing(IsingLattice lat, int MC_N)
{
    int L = lat.lattice.size();
    std::pair<float, float> E_limit = lat.Energy_Limit();
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
    while (f - 1 > 1e-8)
    {   
        #pragma omp parallel
        for (int i = 0; i < MC_N; i++)
        {
            //print f with 10 decimal places
            std::cout << "f: " << std::fixed << std::setprecision(10) << f << std::endl;
            int x = rand() % L;
            int y = rand() % L;
            int s0 = lat.lattice[x][y];

            float Old_E = lat.Ising_Energy();

            int s1 = -s0;
            lat.lattice[x][y] = s1;

            float New_E = lat.Ising_Energy();

            if (rand() / (double)RAND_MAX < std::min(1.0, g[Old_E] / g[New_E]))
            {
                g[New_E] *= f;
                hist[New_E] += 1;
            }
            else
            {
                lat.lattice[x][y] = s0;
                g[Old_E] *= f;
                hist[Old_E] += 1;
            }
            if(i % 10 == 0)
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