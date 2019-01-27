#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <json/json.h>

#include "Grid.h"
#include "Axis.h"
#include "Tools.h"

using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::cout;
using std::endl;
using std::cerr;

ReweightHistogram getPibias1D(const vector<Axis>& ax, const string& filename, const size_t dimension, const bool normalize = true) {
    ReweightHistogram Pibias(ax);
    vector<double> cvs(dimension, 0);
    std::cout << "Reading " << filename << std::endl;
    ifstream ifs_traj(filename.c_str());
    string line;
    string token{" "};
    vector<string> fields;
    while (std::getline(ifs_traj, line)) {
        splitString(line, token, fields);
        if (fields.empty()) continue;
        if (fields[0].compare("#") != 0) {
            for (size_t i = 0; i < dimension; ++i) {
                cvs[i] = std::stod(fields[2*i+1]);
            }
            Pibias.store(cvs, 1.0, 0);
        }
        fields.clear();
    }
    if (normalize) {
        Pibias.normalize();
    }
    return Pibias;
}

HistogramValue genHarmonicPotential(const vector<Axis>& ax, const vector<double>& centre, double K) {
    HistogramValue potential(ax);
    vector<vector<double>> table = potential.getTable();
    double grid_size = potential.getGridSize();
    size_t dimension = potential.getDimension();
    vector<double> pos(dimension, 0);
    for (size_t i = 0; i < grid_size; ++i) {
        double energy = 0;
        for (size_t j = 0; j < dimension; ++j) {
            pos[j] = table[j][i];
            double dist = ax[j].distance(centre[j], pos[j]);
            energy += 0.5 * K * dist * dist;
        }
        potential.set(pos, energy);
    }
    return potential;
}

template <typename T>
T calcError(const vector<T>& v1, const vector<T>& v2) {
    T err = T();
    for (size_t i = 0; i < v1.size(); ++i) {
        err += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    err = std::sqrt(err / v1.size());
    return err;
}

void wham(const Json::Value& obj) {
    // Read temperature and boltzmann constant
    double temperature = obj["temperature"].asDouble();
    const double boltzmann = obj["boltzmann"].asDouble();
    const double kbt = temperature * boltzmann;
    const double beta = 1.0 / kbt;

    // Other settings
    const string output_prefix = obj["output"].asString();
    int maxIteration = 1000000;
    double tolerance = 1e-7;
    if (obj.isMember("max_iteration")) maxIteration = obj["max_iteration"].asInt();
    if (obj.isMember("tolerance")) tolerance = obj["tolerance"].asDouble();

    // Grid settings
    const Json::Value& grid_config = obj["grid"];
    const unsigned dimension = grid_config.size();
    vector<Axis> a1(dimension);
    for (unsigned i = 0; i < dimension; ++i) {
        double lower = grid_config[i]["lower"].asDouble();
        double upper = grid_config[i]["upper"].asDouble();
        unsigned bins = grid_config[i]["bins"].asUInt();
        bool periodic = false;
        if (obj.isMember("periodic")) periodic = obj["periodic"].asBool();
        a1[i] = Axis(lower, upper, bins, periodic);
    }

    // Number of windows
    const Json::Value& windows = obj["windows"];
    const unsigned N = windows.size();
    // Construct biased histograms and potentials
    vector<ReweightHistogram> PiList(N);
    vector<HistogramValue> potentials(N);
    // Determine line format of the colvars traj file
    for (unsigned i = 0; i < N; ++i) {
        // Get biased distribution from colvars traj file
        PiList[i] = getPibias1D(a1, windows[i]["filename"].asString(), dimension, true);
        // Determine the centers of bias potentials
        const Json::Value& centers = windows[i]["centers"];
        if (centers.size() != dimension) {
            cerr << "Dimension mismatches!" << endl;
            std::abort();
        }
        vector<double> c_i(dimension);
        for (unsigned j = 0; j < dimension; ++j) {
            c_i[j] = centers[j].asDouble();
        }
        const double k_i = windows[i]["spring constant"].asDouble();
        potentials[i] = genHarmonicPotential(a1, c_i, k_i);
    }

    // Total unbiased histogram
    HistogramValue Pu(a1);
    // F_i
    FILE* outfile_F;
    outfile_F = fopen(string(output_prefix + "_F.dat").c_str(), "w");
    fprintf(outfile_F, "# iteration error");
    for (unsigned i = 0; i < N; ++i) {
        fprintf(outfile_F, " F[%d]", i);
    }
    fprintf(outfile_F, "\n");
    vector<double> F(N, 0);
    // F_i in last iteration step, used to compute the error
    vector<double> F_old(N, 0);

    // WHAM iteration
    for (int iter = 0; iter < maxIteration; ++iter) {
        // Backup old F value
        F_old = F;
        HistogramValue num(a1);
        HistogramValue denom(a1);
        for (unsigned j = 0; j < N; ++j) {
            double f_j = F[j];
            auto n_j = PiList[j].getTotalCount();
            num = num + n_j * PiList[j];
            // Compute n_j * exp(-beta (w_j_xi - f_j))
            // w_j_xi is the j-th bias potential at reaction coordinate xi
            HistogramValue tmp = potentials[j];
            tmp.applyFunction([beta, f_j, n_j](double x){return n_j * std::exp(-beta * (x - f_j));});
            denom = denom + tmp;
        }
        Pu = num / denom;
        for (unsigned j = 0; j < N; ++j) {
            // Compute exp(-beta (w_j_xi))
            HistogramValue tmp = potentials[j];
            tmp.applyFunction([beta](double x){return std::exp(-beta * x);});
            const vector<double>& raw_p = Pu.getRawData();
            const vector<double>& raw_weight = tmp.getRawData();
            // Compute integral ∫Pu(xi)exp(-beta (w_j_xi))dxi
            // which is exp(-beta * f_j)
            double integral = 0;
            for (size_t i = 0; i < raw_p.size(); ++i) {
                integral += raw_p[i] * raw_weight[i];
            }
            F[j] = -kbt * std::log(integral);
        }
        // Measure convergence of WHAM
        double error = calcError(F, F_old);
        fprintf(outfile_F, " %d %10.7lf", iter, error);
        for (size_t j = 0; j < N; ++j) {
            fprintf(outfile_F, " %10.7lf", F[j]);
        }
        fprintf(outfile_F, "\n");
        if (error < tolerance) {
            break;
        } 
        if (iter == maxIteration - 1) {
            cout << "Warning: Maximum steps of iteration reach!" << '\n';
            if (error >= tolerance) {
                cout << "Warning: Numerical solution of the WHAM equations is still unstable!" << '\n';
            }
        }
    }
    Pu.normalize();
    Pu.writeToFile(string(output_prefix + "_Pu.dat"));
    HistogramValue FES = convertToFreeEnergy(Pu, kbt);
    FES.writeToFile(string(output_prefix + "_FES.dat"));

    // Output PiList
    HistogramView hv_pi(PiList);
    hv_pi.writeToFile(string(output_prefix + "_Pibiased.dat"));

    fclose(outfile_F);
}

int main(int argc, char* argv[]) {
    ifstream ifs_config(argv[1]);
    if (!ifs_config.is_open()) {
        cerr << "Error reading config file!" << endl;
        return 1;
    }
    Json::CharReaderBuilder reader;
    Json::Value wham_config;
    string json_error;
    Json::parseFromStream(reader, ifs_config, &wham_config, &json_error);
    wham(wham_config);
    return 0;
}
