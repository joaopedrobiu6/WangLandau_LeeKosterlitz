#include "utils.h"

void save_data(std::map<int, double> &data, std::string filename)
{
    std::ofstream file;
    file.open(filename);
    for (const auto &entry : data)
    {
        file << entry.first << "\t" << entry.second << std::endl;
    }
    file.close();
};

std::string GetFilename(int L, float J, int q)
{
    std::string J_str = std::to_string(J);
    J_str = J_str.substr(0, J_str.find(".") + 3);
    std::string filename = "L" + std::to_string(L) + "_J" + J_str + "_q" + std::to_string(q) + ".txt";
    return filename;
};