#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

#include <json/json.h>

#include "Grid.h"
#include "Axis.h"

using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;
using std::cout;
using std::endl;
using std::cerr;

int main() {
    vector<Axis> a1{Axis(-180, 180, 360), Axis(-180, 180, 360)};
    HistogramFiles hf(a1, "eABF_WHAM");
    ifstream ifs_traj("../alad_300ns.colvars.traj");
    if (!ifs_traj.is_open()) {
        cerr << "Not opened!" << endl;
        cerr << strerror(errno) << endl;
    }
    unsigned long int step;
    double dihed_a, r_dihed_a, dihed_b, r_dihed_b;
    string line;
    while (std::getline(ifs_traj, line)) {
        if (std::sscanf(line.c_str(), "%lu %lf %lf %lf %lf", &step, &dihed_a, &r_dihed_a, &dihed_b, &r_dihed_b)==5) {
            hf.store(vector<double>{r_dihed_a, r_dihed_b}, line);
        }
    }
    hf.saveInfo("windows_info.dat");
    return 0;
}
