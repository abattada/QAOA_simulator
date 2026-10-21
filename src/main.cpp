#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "qaoa.hpp"
#include "state.hpp"

// p=1 QAOA for MaxCut with an exhaustive (gamma, beta) grid search.
//
// Per grid point: build e^{-i gamma C} H^n |0> in one pass from the cut table,
// apply Rx(2 beta) to every qubit (cache-blocked), then compute <C> and the
// probability of an optimal cut exactly, and optionally draw `shots` exact
// samples. One state vector is reused for the whole grid; all parallelism is
// inside the kernels.

static bool parseInt(const char* s, long long& out) {
    try {
        size_t pos = 0;
        out = std::stoll(s, &pos);
        return pos == std::string(s).size();
    } catch (...) {
        return false;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <graph_name> <shots> <steps> [seed]\n"
                  << "  <graph_name>  file in ./graph/\n"
                  << "  <shots>       exact samples per grid point (0 = no sampling)\n"
                  << "  <steps>       grid resolution: steps x steps points over (gamma, beta) in [0, pi]^2\n"
                  << "  [seed]        base RNG seed for sampling (default 1)\n";
        return 1;
    }

    const std::string inputName = argv[1];
    long long shotsArg = 0, stepsArg = 0, seedArg = 1;
    if (!parseInt(argv[2], shotsArg) || shotsArg < 0 || shotsArg > 100000000) {
        std::cerr << "Error: <shots> must be an integer in 0..1e8\n";
        return 1;
    }
    if (!parseInt(argv[3], stepsArg) || stepsArg < 1 || stepsArg > 10000) {
        std::cerr << "Error: <steps> must be an integer in 1..10000\n";
        return 1;
    }
    if (argc == 5 && !parseInt(argv[4], seedArg)) {
        std::cerr << "Error: [seed] must be an integer\n";
        return 1;
    }
    const int shots = static_cast<int>(shotsArg);
    const int steps = static_cast<int>(stepsArg);
    const uint64_t baseSeed = static_cast<uint64_t>(seedArg);

    // Read the graph before creating any output file
    Graph graph;
    std::string err;
    if (!loadGraph("./graph/" + inputName, graph, err)) {
        std::cerr << "Error: " << err << "\n";
        return 1;
    }
    const int n = graph.n;
    const int m = static_cast<int>(graph.edges.size());

    // Output file results/<graph>/<shots>_<steps>_<run>, first unused run index
    namespace fs = std::filesystem;
    fs::path outputDir = fs::path("results") / inputName;
    fs::create_directories(outputDir);
    fs::path candidate;
    for (int run = 0;; run++) {
        candidate = outputDir / (std::to_string(shots) + "_" + std::to_string(steps) + "_" + std::to_string(run));
        if (!fs::exists(candidate)) break;
    }
    if (!freopen(candidate.string().c_str(), "w", stdout)) {
        std::cerr << "Error: could not open output file " << candidate.string() << "\n";
        return 1;
    }

    const double t0 = omp_get_wtime();
    int maxCut = 0;
    const std::vector<uint16_t> cut = buildCutTable(graph, maxCut);
    const double tCut = omp_get_wtime() - t0;

    State state;
    if (!allocState(state, n)) {
        std::cerr << "Error: could not allocate the state vector for " << n << " qubits\n";
        return 1;
    }

    double bestExpected = -1.0, bestGamma = 0.0, bestBeta = 0.0, pOptAtBest = 0.0;
    int bestSampled = -1;

    for (int ig = 0; ig < steps; ig++) {
        for (int ib = 0; ib < steps; ib++) {
            const double gamma = (steps == 1) ? 0.0 : M_PI * double(ig) / (steps - 1);
            const double beta  = (steps == 1) ? 0.0 : M_PI * double(ib) / (steps - 1);

            applyCostLayer(state, cut, m, gamma);
            applyMixingLayer(state, n, beta);

            // Deterministic per-point seed, independent of the thread count
            std::mt19937_64 rng(splitmix64(baseSeed) ^ splitmix64(uint64_t(ig) * steps + ib + 1));
            const Measurement r = measure(state, cut, maxCut, shots, rng);

            if (r.expectedCut > bestExpected) {
                bestExpected = r.expectedCut;
                bestGamma = gamma;
                bestBeta = beta;
                pOptAtBest = r.pOpt;
            }
            bestSampled = std::max(bestSampled, r.sampledMaxCut);

            std::cout << "Gamma: " << gamma << ", Beta: " << beta
                      << std::setprecision(12)
                      << ", Expected Cut: " << r.expectedCut
                      << ", P(opt): " << r.pOpt;
            if (shots > 0) std::cout << ", Sampled Max Cut: " << r.sampledMaxCut;
            std::cout << std::setprecision(6) << "\n";
        }
    }
    destroyState(state);
    const double wall = omp_get_wtime() - t0;

    // Summary: once into the results file (stdout), once to the terminal (stderr)
    std::stringstream summary;
    summary << std::setprecision(6);
    summary << "=========================================\n";
    summary << " QAOA (p=1) with " << n << " qubits.\n";
    summary << " Input graph: " << inputName << " (edges=" << m << ")\n";
    summary << " Grid search steps = " << steps << ", shots = " << shots
            << ", seed = " << baseSeed << "\n";
    summary << "-----------------------------------------\n";
    summary << " True max cut      = " << maxCut << "\n";
    summary << " Best expected cut = " << bestExpected << "\n";
    if (maxCut > 0) summary << " Approx ratio      = " << bestExpected / maxCut << "\n";
    else            summary << " Approx ratio      = n/a\n";
    summary << " P(opt) at best    = " << pOptAtBest << "\n";
    if (shots > 0)  summary << " Best sampled cut  = " << bestSampled << "\n";
    summary << " Best gamma        = " << bestGamma << "\n";
    summary << " Best beta         = " << bestBeta << "\n";
    summary << "=========================================\n";

    std::cout << summary.str() << std::flush;
    // Timing goes to the terminal only, so result files are reproducible
    std::cerr << summary.str()
              << " OpenMP threads = " << omp_get_max_threads()
              << ", wall time = " << std::fixed << std::setprecision(2) << wall
              << " s (cut table " << tCut << " s)\n";
    return 0;
}
