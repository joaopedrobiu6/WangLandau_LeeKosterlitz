#include "WangLandau.h"

bool isFlat(const std::map<int, int> &hist, double h_tol)
{
    // Remove entries with zero counts
    std::vector<int> values;
    for (const auto &entry : hist)
    {
        if (entry.second > 0)
        {
            values.push_back(entry.second);
        }
    }

    if (values.empty())
    {
        return false;
    }

    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double min_value = *std::min_element(values.begin(), values.end());

    if (min_value < h_tol * mean)
    {
        // std::cout << "it's not flat yet" << std::endl;
        return false;
    }
    return true;
}

std::map<int, double> WangLandauPotts(PottsLattice lat, int MC_N, int q, double f_tol, double h_tol, bool NoLog)
{
    srand(time(NULL));
    ;
    int L = lat.lattice.size();
    std::pair<float, float> E_limit = lat.Energy_Limit();
    float E_min = E_limit.first;
    float E_max = E_limit.second;
    float num_bins = -E_min + 1;

    std::map<int, double> lng;
    std::map<int, double> g;
    std::map<int, int> hist; // Histogram H(E)

    for (int E = E_min; E <= E_max; E += 1)
    {
        if (NoLog)
        {
            g[E] = 1.0;
        }
        else
        {
            lng[E] = 0.0;
        }
        hist[E] = 0;
    }

    double f = std::exp(1.0);
    double lnf = std::log(f);

    if (NoLog)
    {
        while (f - 1 > f_tol)
        {
#pragma omp parallel for
            for (int i = 0; i < MC_N; i++)
            {
                int x = rand() % L;
                int y = rand() % L;
                int s0 = lat.lattice[x][y];

                float Old_E = lat.Potts_Energy();

                int s1 = 1 + rand() % q;
                lat.lattice[x][y] = s1;

                float New_E = lat.Potts_Energy();

                if (rand() / (double)RAND_MAX < std::min(1.0, g[Old_E] / g[New_E]))
                {
                    lng[New_E] += lnf;
                    hist[New_E] += 1;
                }
                else
                {
                    lat.lattice[x][y] = s0;
                    g[Old_E] *= f;
                    hist[Old_E] += 1;
                }

                if (i % 1000 == 0)
                {
                    std::cout << "f: " << f << std::endl;
                    if (isFlat(hist, h_tol))
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
    else
    {
        while (lnf > f_tol)
        {
#pragma omp parallel for
            for (int i = 0; i < MC_N; i++)
            {
                int x = rand() % L;
                int y = rand() % L;
                int s0 = lat.lattice[x][y];

                float Old_E = lat.Potts_Energy();

                int s1 = 1 + rand() % q;
                lat.lattice[x][y] = s1;

                float New_E = lat.Potts_Energy();

                if (rand() / (double)RAND_MAX < std::min(1.0, std::exp(lng[Old_E] - lng[New_E])))
                {
                    lng[New_E] += lnf;
                    hist[New_E] += 1;
                }
                else
                {
                    lat.lattice[x][y] = s0;
                    lng[Old_E] += lnf;
                    hist[Old_E] += 1;
                }

                if (i % 1000 == 0)
                {
                    std::cout << "lnf: " << lnf << std::endl;
                    if (isFlat(hist, h_tol))
                    {
                        lnf = lnf / 8;
                        hist.clear();
                        for (int E = E_min; E <= E_max; E += 1)
                        {
                            hist[E] = 0;
                        }
                    }
                }
            }
        }
        return lng;
    }
}