#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <map>

void save_data(std::map<int, double> &data, std::string filename);

std::string GetFilename(int L, float J, int q);