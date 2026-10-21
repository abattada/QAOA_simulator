// Correctness tests for the gate kernels and the fused QAOA kernels.
// Run from the repository root (`make test`); exits non-zero on failure.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "circuit.hpp"
#include "qaoa.hpp"
#include "state.hpp"

static int failures = 0;

static void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "pass" : "FAIL", what.c_str());
    if (!ok) failures++;
}

static double maxDiff(const State& a, const State& b) {
    double d = 0.0;
    for (ull i = 0; i < a.numAmpTotal; i++) {
        d = std::max(d, std::hypot(a.stateVec[i].real - b.stateVec[i].real,
                                   a.stateVec[i].imag - b.stateVec[i].imag));
    }
    return d;
}

static void randomState(State& s, int n, uint64_t seed) {
    createState(s, n);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> d(0.0, 1.0);
    double norm = 0.0;
    for (ull i = 0; i < s.numAmpTotal; i++) {
        s.stateVec[i] = {d(rng), d(rng)};
        norm += s.stateVec[i].real * s.stateVec[i].real + s.stateVec[i].imag * s.stateVec[i].imag;
    }
    norm = std::sqrt(norm);
    for (ull i = 0; i < s.numAmpTotal; i++) {
        s.stateVec[i].real /= norm;
        s.stateVec[i].imag /= norm;
    }
}

static void copyState(State& dst, const State& src, int n) {
    createState(dst, n);
    for (ull i = 0; i < src.numAmpTotal; i++) dst.stateVec[i] = src.stateVec[i];
}

static Graph randomGraph(int n, int m, uint64_t seed) {
    std::vector<std::pair<int, int>> all;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) all.push_back({i, j});
    std::mt19937_64 rng(seed);
    std::shuffle(all.begin(), all.end(), rng);
    Graph g;
    g.n = n;
    g.adj.assign(n, 0);
    for (int k = 0; k < m; k++) {
        // Random orientation, so both CNOT directions are exercised
        auto [a, b] = all[k];
        if (rng() & 1) std::swap(a, b);
        g.edges.push_back({a, b});
        g.adj[a] |= 1ULL << b;
        g.adj[b] |= 1ULL << a;
    }
    return g;
}

// Gate-by-gate reference: H^n, then per edge CNOT(i,j) P_j(-gamma) CNOT(i,j),
// then Rx(2 beta) on every qubit.
static void referenceCircuit(State& s, const Graph& g, double gamma, double beta, bool mixing) {
    createState(s, g.n);
    for (int q = 0; q < g.n; q++) hGate(s, q);
    for (auto [i, j] : g.edges) {
        cnotGate(s, i, j);
        phaseGate(s, j, -gamma);
        cnotGate(s, i, j);
    }
    if (mixing)
        for (int q = 0; q < g.n; q++) rxGate(s, q, 2.0 * beta);
}

static void fusedCircuit(State& s, const Graph& g, const std::vector<uint16_t>& cut,
                         double gamma, double beta, bool mixing) {
    allocState(s, g.n);
    applyCostLayer(s, cut, static_cast<int>(g.edges.size()), gamma);
    if (mixing) applyMixingLayer(s, g.n, beta);
}

static void bruteForce(const State& s, const std::vector<uint16_t>& cut, int maxCut,
                       double& expected, double& pOpt) {
    expected = 0.0;
    pOpt = 0.0;
    for (ull x = 0; x < s.numAmpTotal; x++) {
        double p = s.stateVec[x].real * s.stateVec[x].real + s.stateVec[x].imag * s.stateVec[x].imag;
        expected += p * cut[x];
        if (cut[x] == maxCut) pOpt += p;
    }
}

static void testCnotBasics() {
    std::printf("1. CNOT on basis states\n");
    State s;
    createState(s, 2);
    s.stateVec[0].real = 0.0;
    s.stateVec[2].real = 1.0;                   // |q1=1, q0=0>
    cnotGate(s, 1, 0);
    check(s.stateVec[3].real == 1.0 && s.stateVec[2].real == 0.0, "CNOT(c=1,t=0)|10> = |11>");
    destroyState(s);

    State b;
    createState(b, 2);
    hGate(b, 0);
    cnotGate(b, 0, 1);
    const double r = 1.0 / std::sqrt(2.0);
    check(std::fabs(b.stateVec[0].real - r) < 1e-15 && std::fabs(b.stateVec[3].real - r) < 1e-15 &&
          b.stateVec[1].real == 0.0 && b.stateVec[2].real == 0.0,
          "H(0); CNOT(0,1)|00> = Bell state");
    destroyState(b);
}

static void testCnotPermutation() {
    std::printf("2. CNOT equals the index permutation (n = 10)\n");
    const int n = 10;
    const std::pair<int, int> pairs[] = {{0, 1}, {1, 0}, {0, 9}, {9, 0}, {3, 7}, {7, 3}, {4, 5}};
    for (auto [c, t] : pairs) {
        State in, out;
        randomState(in, n, 100 + c * 16 + t);
        copyState(out, in, n);
        cnotGate(out, c, t);
        double d = 0.0;
        for (ull x = 0; x < in.numAmpTotal; x++) {
            ull src = ((x >> c) & 1) ? x ^ (1ULL << t) : x;
            d = std::max(d, std::hypot(out.stateVec[x].real - in.stateVec[src].real,
                                       out.stateVec[x].imag - in.stateVec[src].imag));
        }
        check(d == 0.0, "control " + std::to_string(c) + ", target " + std::to_string(t));
        destroyState(in);
        destroyState(out);
    }
}

static void testCostLayer(const Graph& g, const std::string& name) {
    int maxCut = 0;
    auto cut = buildCutTable(g, maxCut);
    for (double gamma : {0.0, 0.37, 1.1, 2.9}) {
        State ref, fused;
        referenceCircuit(ref, g, gamma, 0.0, false);
        fusedCircuit(fused, g, cut, gamma, 0.0, false);
        double d = maxDiff(ref, fused);
        char buf[128];
        std::snprintf(buf, sizeof buf, "%s, gamma = %.2f: max |diff| = %.1e", name.c_str(), gamma, d);
        check(d < 1e-12, buf);
        destroyState(ref);
        destroyState(fused);
    }
}

static void testMixingLayer(int n, int blockBits) {
    for (double beta : {0.0, 0.41, 1.3, 2.7}) {
        State ref, blk;
        randomState(ref, n, 7 + n);
        copyState(blk, ref, n);
        for (int q = 0; q < n; q++) rxGate(ref, q, 2.0 * beta);
        applyMixingLayer(blk, n, beta, blockBits);
        double d = maxDiff(ref, blk);
        char buf[128];
        std::snprintf(buf, sizeof buf, "n = %d, blockBits = %d, beta = %.2f: max |diff| = %.1e",
                      n, blockBits, beta, d);
        check(d < 1e-13, buf);
        destroyState(ref);
        destroyState(blk);
    }
}

static void testEndToEnd(const Graph& g) {
    std::printf("5. <C> and P(opt) against the gate-by-gate circuit (qml_bi_4)\n");
    int maxCut = 0;
    auto cut = buildCutTable(g, maxCut);
    const int m = static_cast<int>(g.edges.size());
    check(maxCut == m, "brute-force MaxCut of K_{4,4} = 16");
    std::mt19937_64 rng(1);
    for (auto [gamma, beta] : std::vector<std::pair<double, double>>{
             {0.3, 0.7}, {1.1, 0.5}, {0.5, 0.4}, {2.4, 2.2}}) {
        State ref, fused;
        referenceCircuit(ref, g, gamma, beta, true);
        fusedCircuit(fused, g, cut, gamma, beta, true);
        double e, p;
        bruteForce(ref, cut, maxCut, e, p);
        Measurement r = measure(fused, cut, maxCut, 0, rng);
        char buf[160];
        std::snprintf(buf, sizeof buf, "gamma = %.2f, beta = %.2f: <C> = %.6f (ref %.6f), P(opt) = %.6f (ref %.6f)",
                      gamma, beta, r.expectedCut, e, r.pOpt, p);
        check(std::fabs(r.expectedCut - e) < 1e-10 && std::fabs(r.pOpt - p) < 1e-12, buf);
        destroyState(ref);
        destroyState(fused);
    }
    for (auto [gamma, beta] : std::vector<std::pair<double, double>>{{0.0, 1.2}, {1.7, 0.0}}) {
        State s;
        fusedCircuit(s, g, cut, gamma, beta, true);
        Measurement r = measure(s, cut, maxCut, 0, rng);
        char buf[128];
        std::snprintf(buf, sizeof buf, "gamma = %.1f, beta = %.1f: <C> = %.12f = m/2", gamma, beta, r.expectedCut);
        check(std::fabs(r.expectedCut - m / 2.0) < 1e-12, buf);
        destroyState(s);
    }
}

static void testMarginals(const Graph& g, const std::string& name) {
    int maxCut = 0;
    auto cut = buildCutTable(g, maxCut);
    double worst = 0.0;
    for (auto [gamma, beta] : std::vector<std::pair<double, double>>{{0.3, 0.7}, {1.1, 0.5}, {1.8188, 2.3149}}) {
        State s;
        fusedCircuit(s, g, cut, gamma, beta, true);
        for (int q = 0; q < g.n; q++) {
            double p0 = 0.0;
            for (ull x = 0; x < s.numAmpTotal; x++)
                if (!((x >> q) & 1))
                    p0 += s.stateVec[x].real * s.stateVec[x].real + s.stateVec[x].imag * s.stateVec[x].imag;
            worst = std::max(worst, std::fabs(p0 - 0.5));
        }
        destroyState(s);
    }
    char buf[128];
    std::snprintf(buf, sizeof buf, "%s: max |P(q=0) - 1/2| = %.1e", name.c_str(), worst);
    check(worst < 1e-12, buf);
}

static void testSampling(const Graph& g) {
    std::printf("7. exact sampling (qml_bi_4, 20000 shots)\n");
    int maxCut = 0;
    auto cut = buildCutTable(g, maxCut);
    State s;
    fusedCircuit(s, g, cut, 1.1, 0.5, true);

    std::map<int, double> exact;
    for (ull x = 0; x < s.numAmpTotal; x++)
        exact[cut[x]] += s.stateVec[x].real * s.stateVec[x].real + s.stateVec[x].imag * s.stateVec[x].imag;

    const int shots = 20000;
    std::vector<ull> a, b;
    std::mt19937_64 rng1(splitmix64(42)), rng2(splitmix64(42));
    Measurement ra = measure(s, cut, maxCut, shots, rng1, &a);
    measure(s, cut, maxCut, shots, rng2, &b);
    check(a == b, "same seed gives identical samples");

    std::map<int, int> counts;
    for (ull x : a) counts[cut[x]]++;
    bool ok = true;
    double worstZ = 0.0;
    for (auto [c, p] : exact) {
        const double mean = shots * p;
        const double sigma = std::sqrt(shots * p * (1 - p));
        const double z = std::fabs(counts[c] - mean) / std::max(sigma, 1.0);
        worstZ = std::max(worstZ, z);
        if (z > 4.0) ok = false;
    }
    for (auto [c, k] : counts)
        if (exact[c] == 0.0 && k > 0) ok = false;
    char buf[128];
    std::snprintf(buf, sizeof buf, "cut-value histogram matches |psi|^2 (worst deviation %.2f sigma)", worstZ);
    check(ok, buf);
    int best = 0;
    for (ull x : a) best = std::max<int>(best, cut[x]);
    check(ra.sampledMaxCut == best, "sampledMaxCut = max cut over the samples");
    destroyState(s);
}

int main() {
    Graph bi4;
    std::string err;
    if (!loadGraph("graph/qml_bi_4", bi4, err)) {
        std::printf("cannot load graph/qml_bi_4 (run from the repo root): %s\n", err.c_str());
        return 1;
    }
    const Graph rnd = randomGraph(10, 22, 2026);

    testCnotBasics();
    testCnotPermutation();

    std::printf("3. fused cost layer equals H^n + CNOT.P(-gamma).CNOT per edge\n");
    testCostLayer(bi4, "qml_bi_4");
    testCostLayer(rnd, "random G(10, 22)");

    std::printf("4. blocked mixing layer equals sequential Rx\n");
    testMixingLayer(10, 3);    // high-qubit groups of 3, 3, 1
    testMixingLayer(12, 4);    // groups of 3, 3, 2
    testMixingLayer(16, 14);   // block shrinks so every thread gets a block
    testMixingLayer(8, 14);    // single block

    testEndToEnd(bi4);

    std::printf("6. every single-qubit marginal of the QAOA state is 1/2\n");
    testMarginals(bi4, "qml_bi_4");
    testMarginals(rnd, "random G(10, 22)");

    testSampling(bi4);

    std::printf(failures ? "\n%d check(s) FAILED\n" : "\nall checks passed\n", failures);
    return failures ? 1 : 0;
}
