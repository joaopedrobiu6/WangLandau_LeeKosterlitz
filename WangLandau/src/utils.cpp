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
}