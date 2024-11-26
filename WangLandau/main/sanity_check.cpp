#include <iostream>
#include <vector>
#include <set>

int En (std::vector<std::vector<int>> lat, int L) {
  int i,j;
  int energy=0;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      energy-=((int)(lat[i][j]==lat[i][(j+1)%L])+ (int)(lat[i][j]==lat[(i+1)%L][j]));
    }
  }
  return energy;
}

int energy(std::vector<std::vector<int>> lattice, int L, float J)
{
    double E = 0.0;

    // maybe we can iterate over only half of the lattice to avoid double counting
    // for (int i = 0; i < L; ++i)
    // {
    //     for (int j = 0; j < L; ++j)
    //     {
    //         int s0 = lattice[i][j];
    //         for (auto [dx, dy] : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1}, std::pair{0, 1}})
    //         {
    //             int x = (i + dx + L) % L;
    //             int y = (j + dy + L) % L;
    //             int s1 = lattice[x][y];
    //             E -= (s1 == s0) ? 1 : 0;
    //         }
    //     }
    // }

    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            int s0 = lattice[i][j];

            // Check only forward neighbors (right and down)
            for (auto [dx, dy] : {std::pair{1, 0}, std::pair{0, 1}})
            {
                int x = (i + dx) % L;
                int y = (j + dy) % L;
                int s1 = lattice[x][y];

                // Increment energy if neighbors are the same
                E -= (s1 == s0) ? 1 : 0;
            }
        }
    }

    return J * E;
}

std::vector<std::vector<std::vector<int>>> generateMatrices()
{
    std::vector<std::vector<std::vector<int>>> matrices;

    // Iterate over all possible values for each matrix element
    #pragma omp parallel for
    for (int a = 1; a <= 8; ++a)
    {
        for (int b = 1; b <= 8; ++b)
        {
            for (int c = 1; c <= 8; ++c)
            {
                for (int d = 1; d <= 8; ++d)
                {
                    for (int e = 1; e <= 8; ++e)
                    {
                        for (int f = 1; f <= 8; ++f)
                        {
                            for (int g = 1; g <= 8; ++g)
                            {
                                for (int h = 1; h <= 8; ++h)
                                {
                                    for (int i = 1; i <= 8; ++i)
                                    {
                                        std::vector<std::vector<int>> matrix = {
                                            {a, b, c},
                                            {e, d, f},
                                            {g, h, i}};
                                        matrices.push_back(matrix);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return matrices;
}

int main()
{
    std::vector<std::vector<std::vector<int>>> matrices = generateMatrices();
    std::cout << "Number of matrices: " << matrices.size() << std::endl;
    std::set<double> energyLevels;
    std::set<double> energyLevels2;

    std::vector<std::vector<int>> lattice;
    std::vector<std::vector<int>> case7;
    std::vector<std::vector<int>> case72;


    // set lattice to be ((1, 2), (3, 2))
    std::cout << "finished creating matrices" << std::endl;
    std::cout << "starting energy calculation" << std::endl;

    int E = 0;
    int e = 0;

    #pragma omp parallel for
    for (int i = 0; i < matrices.size(); ++i)
    {
        if (i % 1000000 == 0)
        {
            std::cout << "i = " << i << std::endl;
        }
        lattice = matrices[i];
        E = energy(lattice, 3, 1);
        e = En(lattice, 3);
        
        if (E == -7)
        {
            case7 = lattice;
        }
        if (e == -9)
        {
            case72 = lattice;
        }

        energyLevels.insert(E);
        energyLevels2.insert(e);
    }

    std::cout << "Number of distinct energy levels: " << energyLevels.size() << std::endl;
    // print the energy levels
    for (double E : energyLevels)
    {
        std::cout << E << std::endl;
    }

    std::cout << "Number of distinct energy levels2: " << energyLevels2.size() << std::endl;
    // print the energy levels
    for (double e : energyLevels2)
    {
        std::cout << e << std::endl;
    }

    // print the lattice with energy -7
    std::cout << "Lattice with energy -7:" << std::endl;
    for (int i = 0; i < case7.size(); ++i)
    {
        for (int j = 0; j < case7[i].size(); ++j)
        {
            std::cout << case7[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print the lattice with energy -7
    std::cout << "Lattice with energy2 -7:" << std::endl;
    for (int i = 0; i < case72.size(); ++i)
    {
        for (int j = 0; j < case72[i].size(); ++j)
        {
            std::cout << case72[i][j] << " ";
        }
        std::cout << std::endl;
    }



    return 0;
}