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

void info_dump(std::string folder, int MC_N, int L, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor)
{
    std::ofstream info_file(folder + "/info.txt");
    info_file << "MC samples: " << MC_N << std::endl;
    info_file << "Lattice size: " << L << std::endl;
    info_file << "q: " << q << std::endl;
    info_file << "f tolerance: " << f_tol << std::endl;
    info_file << "Histogram tolerance: " << h_tol << std::endl;
    if (NoLog)
    {
        info_file << "Exp form" << std::endl;
    }
    else
    {
        info_file << "Log form" << std::endl;
    }
    info_file << "Sample Check Interval: " << sampleInterval << std::endl;
    info_file << "f update factor: " << f_factor << std::endl;
    info_file.close();
}

void info_print(int MC_N, int L, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor)
{
    std::cout << "\n*******************************************" << std::endl;
    std::cout << "\nStarting Wang-Landau simulation" << std::endl;
    std::cout << "MC samples: " << MC_N << std::endl;
    std::cout << "Lattice size: " << L << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "f tolerance: " << f_tol << std::endl;
    std::cout << "Histogram tolerance: " << h_tol << std::endl;
    if (NoLog)
    {
        std::cout << "Exp form" << std::endl;
    }
    else
    {
        std::cout << "Log form" << std::endl;
    }
    std::cout << "Sample Check Interval: " << sampleInterval << std::endl;
    std::cout << "f update factor: " << f_factor << std::endl;
    std::cout << "\n*******************************************\n"
              << std::endl;
}

std::map<int, double> WangLandauPotts(PottsLattice lat, int MC_N, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor, std::string filename)
{
    int L = lat.lattice.size();

    info_print(MC_N, L, q, f_tol, h_tol, NoLog, sampleInterval, f_factor);

    // create output folder if it doesn't exist with name equal to filename without extension and in result folder
    std::string folder = "results/" + filename.substr(0, filename.find_last_of("."));
    std::string command = "mkdir -p " + folder;
    system(command.c_str());

    // create info file in output folder with parameters
    info_dump(folder, MC_N, L, q, f_tol, h_tol, NoLog, sampleInterval, f_factor);

    srand(time(NULL));
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

                if (i % sampleInterval == 0)
                {
                    std::cout << "f: " << f << std::endl;
                    if (isFlat(hist, h_tol))
                    {
                        f = pow(f, 1.0 / f_factor);
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

                if (i % sampleInterval == 0)
                {
                    // std::cout << "lnf: " << lnf << std::endl;
                    if (isFlat(hist, h_tol))
                    {
                        lnf = lnf / f_factor;
                        hist.clear();
                        for (int E = E_min; E <= E_max; E += 1)
                        {
                            hist[E] = 0;
                        }
                    }
                }
            }
        }
        filename = folder + "/" + filename;
        save_data(lng, filename);
        return lng;
    }
}