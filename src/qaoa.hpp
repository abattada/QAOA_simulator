#ifndef _QAOA_HPP_
#define _QAOA_HPP_

#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "state.hpp"

// Fused p=1 QAOA kernels for MaxCut.
//
// Conventions (the gate-by-gate equivalent is checked in tests/test_qaoa.cpp):
//   cost layer   : |x> -> e^{-i gamma cut(x)} |x>,  applied to H^n |0...0>
//                  == per edge CNOT(i,j) . P_j(-gamma) . CNOT(i,j)
//   mixing layer : Rx(2 beta) on every qubit,
//                  Rx(2b) = [[cos b, -i sin b], [-i sin b, cos b]]

struct Graph {
    int n = 0;                                   // number of nodes = qubits
    std::vector<std::pair<int, int>> edges;
    std::vector<ull> adj;                        // adj[i] = neighbour bitmask of node i
};

// Read "<n> <m> then m lines 'u v'". The graph must be simple (no self loops,
// no repeated edges) and have 1..32 nodes. On failure returns false and sets err.
bool loadGraph(const std::string& path, Graph& g, std::string& err);

// cut(x) for every bitstring x in [0, 2^n). Also returns the exact MaxCut.
std::vector<uint16_t> buildCutTable(const Graph& g, int& maxCut);

// Overwrite `state` with e^{-i gamma C} H^n |0...0>, i.e. psi[x] = e^{-i gamma cut(x)} / sqrt(2^n).
// One pass over the state; m is the number of edges (largest possible cut value).
void applyCostLayer(State& state, const std::vector<uint16_t>& cut, int m, double gamma);

// Rx(2 beta) on qubits 0..n-1. Qubits below `blockBits` are applied together
// inside L2-sized blocks (one pass); higher qubits go in groups of up to three
// (one pass per group).
void applyMixingLayer(State& state, int n, double beta, int blockBits = 14);

struct Measurement {
    double expectedCut = 0.0;   // <C> = sum_x |psi_x|^2 cut(x)
    double pOpt = 0.0;          // probability of measuring a maximum cut
    int sampledMaxCut = -1;     // best cut among `shots` exact samples (-1 if shots == 0)
};

// Exact <C> and P(opt) in one pass; if shots > 0, draws `shots` bitstrings
// from the exact joint distribution |psi_x|^2 (inverse CDF, one more pass).
// The sampled basis indices are written to *samples when it is non-null.
Measurement measure(const State& state, const std::vector<uint16_t>& cut, int maxCut,
                    int shots, std::mt19937_64& rng, std::vector<ull>* samples = nullptr);

inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

#endif
