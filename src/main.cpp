#include <cmath>
#include <string>
#include <fstream>
#include <random> 
#include <filesystem> 
#include <omp.h> 
#include "state.hpp"
#include "circuit.hpp"

using namespace std;

// (1) Perform H gate on every node
void initHGates(State& state, const int& numQubits) {
    for(int q = 0; q < numQubits; q++){
        hGate(state, q);
    }
}

// (2) Apply cost Hamiltonian for all edges
void applyCostHamiltonian(State& state, const std::vector<std::pair<int,int>>& edges, double gamma) {
    for (auto & e : edges) {
        int i = e.first;
        int j = e.second;
        cnotGate(state, i, j);
        rzGate(state, j, 2.0 * gamma);
        cnotGate(state, i, j);
    }
}

// (3) Apply mixing Hamiltonian: H_B = Σ X_i -> each qubit rotateX(2 * beta)
void applyMixingHamiltonian(State& state, const int& numQubits, double beta) {
    for (int q=0; q<numQubits; q++) {
        rxGate(state, q, 2.0 * beta);
    }
}

// (4) Calculate cut value from measurement result
int calcCutValue(const std::vector<int>& measureRes,
                 const std::vector<std::pair<int,int>>& edges) {
    int cutCount = 0;
    for (auto & e : edges) {
        if (measureRes[e.first] != measureRes[e.second]) {
            cutCount++;
        }
    }
    return cutCount;
}

// (5) Evaluate circuit via multiple shots, returning average cut value
double evaluateQaoaCircuit(State& state,
                           const std::vector<std::pair<int, int>>& edges,
                           int shots, 
                           int numQubits) {
    int maxCut = 0;

    // 1) Precompute probability of measuring |0> for each qubit, in parallel
    std::vector<double> prob0(numQubits, 0.0);

    #pragma omp parallel for
    for (int q = 0; q < numQubits; ++q) {
        double localProb0 = 0.0;

        // Optionally vectorize the loop over amplitudes
        #pragma omp simd
        for (ull i = 0; i < state.numAmpTotal; ++i) {
            // If the q-th qubit of index i is 0, add amplitude²
            if (((i >> q) & 1) == 0) {
                localProb0 += state.stateVec[i].real * state.stateVec[i].real +
                              state.stateVec[i].imag * state.stateVec[i].imag;
            }
        }
        prob0[q] = localProb0;
    }

    // 2) Parallelize over shots to measure each qubit and compute cut
    //    We'll do a manual reduction on maxCut
    #pragma omp parallel
    {
        // Each thread has its own RNG
        std::mt19937_64 rng(std::time(nullptr) + 13 * omp_get_thread_num());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        int localMaxCut = 0;

        #pragma omp for
        for (int shot = 0; shot < shots; ++shot) {

            std::vector<int> measurementResults(numQubits);

            // 2a) For each qubit, sample from precomputed probability prob0[q]
            for (int q = 0; q < numQubits; ++q) {
                double randVal = dist(rng);
                measurementResults[q] = (randVal <= prob0[q]) ? 0 : 1;
            }

            // 2b) Compute the cut value for this measurement
            int tmp = calcCutValue(measurementResults, edges);
            if (tmp > localMaxCut) {
                localMaxCut = tmp;
            }
        }

        // 2c) Reduce into a global maxCut
        #pragma omp critical
        {
            if (localMaxCut > maxCut) {
                maxCut = localMaxCut;
            }
        }
    }

    return maxCut;
}



int main(int argc, char* argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] 
                  << " <input_file> <shots> <steps>" 
                  << std::endl;
        return 1;
    }

    // 用 stoi() 轉成整數
    std::string input_name = argv[1];
    std::string file_path = "./graph/" + input_name;
    int shots = std::stoi(argv[2]);   // Number of measurement repetitions
    int steps = std::stoi(argv[3]);   // Grid search resolution for gamma, beta

    // 建立輸出目錄, e.g. ./results/<input_file_name>/<shots>/
    namespace fs = std::filesystem;
    fs::path output_dir = fs::path("results") 
                          / input_name ;
    fs::create_directories(output_dir);
    fs::path candidate;
    // 找到尚未存在的檔名 steps_<increment>
    int increment = 0;
    while (true) {
        // 您可以自行決定是否要加檔名副檔名，如 .txt 或 .log
        candidate = output_dir 
                             / (std::to_string(shots) + "_" + std::to_string(steps) + "_" 
                                + std::to_string(increment) );
        if (!fs::exists(candidate)) {
            // 找到可用的檔名後，將 stdout 導向該檔案
            freopen(candidate.string().c_str(), "w", stdout);
            break;
        }
        increment++;
    }

    // 接下來開始讀取圖形資料
    FILE *f = fopen(file_path.c_str(), "r");
    if (!f) {
        std::cerr << "Error: could not open file " << file_path << std::endl;
        return 1;
    }

    int numQubits = 0;
    int numEdges = 0;
    fscanf(f, "%d", &numQubits);
    fscanf(f, "%d", &numEdges);

    std::vector<std::pair<int,int>> edges;
    edges.reserve(numEdges);
    for (int i = 0; i < numEdges; i++) {
        int a, b;
        fscanf(f, "%d %d", &a, &b);
        edges.push_back({a, b});
    }
    fclose(f);

    double bestCut = -1;   
    double bestGamma = 0;
    double bestBeta = 0;

    std::vector<std::string> lines(steps * steps);

    // Parallel over (ig, ib)
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int ig = 0; ig < steps; ig++) {
        for (int ib = 0; ib < steps; ib++) {

            double gamma = M_PI * (double(ig) / (steps - 1));
            if (steps == 1) {
                gamma = 0.0;
            }
            double beta  = M_PI * (double(ib) / (steps - 1));
            if (steps == 1) {
                beta = 0.0;
            }

            // (0) Initialize to |000...0>
            State state;
            createState(state, numQubits);

            // (1) Apply Hadamard on all qubits -> get |++++...+>
            initHGates(state, numQubits);

            // (2) Cost Hamiltonian
            applyCostHamiltonian(state, edges, gamma);

            // (3) Mixing Hamiltonian
            applyMixingHamiltonian(state, numQubits, beta);

            // (4) Evaluate with multiple shots
            double maxCut = evaluateQaoaCircuit(state, edges, shots, numQubits);

            // Thread-safe update of the global best
            #pragma omp critical
            {
                if (maxCut > bestCut) {
                    bestCut   = maxCut;
                    bestGamma = gamma;
                    bestBeta  = beta;
                }
            }

            // Build the line you want to output
            std::stringstream ss;
            ss << "Gamma: " << gamma
            << ", Beta: "  << beta
            << ", Max Cut: " << maxCut;

            // Store in our buffer at index (ig * steps + ib)
            lines[ig * steps + ib] = ss.str();

            // Clean up
            destroyState(state);
        }
    }

    // Now, after the parallel region finishes, we output lines in ascending (ig, ib) order
    // This ensures the same iteration order as a serial loop would have printed.
    for (int ig = 0; ig < steps; ig++) {
        for (int ib = 0; ib < steps; ib++) {
            std::cout << lines[ig * steps + ib] << "\n";
        }
    }

    // Output summary
    std::cout << "=========================================\n";
    std::cout << " QAOA (p=1) with " << numQubits << " qubits.\n";
    std::cout << " Input graph: " << input_name 
              << " (edges=" << edges.size() << ")\n";
    std::cout << " Grid search steps = " << steps 
              << ", shots = " << shots << std::endl;
    std::cout << "-----------------------------------------\n";
    std::cout << " Best found cut   = " << bestCut << "\n";
    std::cout << " Best gamma       = " << bestGamma << "\n";
    std::cout << " Best beta        = " << bestBeta  << "\n";
    std::cout << "=========================================\n";

    freopen("/dev/tty", "w", stdout);

    // Output summary
    std::cout << "=========================================\n";
    std::cout << " QAOA (p=1) with " << numQubits << " qubits.\n";
    std::cout << " Input graph: " << input_name 
              << " (edges=" << edges.size() << ")\n";
    std::cout << " Grid search steps = " << steps 
              << ", shots = " << shots << std::endl;
    std::cout << "-----------------------------------------\n";
    std::cout << " Best found cut   = " << bestCut << "\n";
    std::cout << " Best gamma       = " << bestGamma << "\n";
    std::cout << " Best beta        = " << bestBeta  << "\n";
    std::cout << "=========================================\n";

    return 0;
}
