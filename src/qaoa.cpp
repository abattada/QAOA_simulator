#include "qaoa.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>

#include <omp.h>

// ---------------------------------------------------------------- graph I/O

bool loadGraph(const std::string& path, Graph& g, std::string& err) {
    std::ifstream in(path);
    if (!in) { err = "could not open file " + path; return false; }

    long long n = 0, m = 0;
    if (!(in >> n >> m)) { err = "malformed graph file header in " + path; return false; }
    if (n < 1 || n > 32) { err = "number of nodes must be in 1..32, got " + std::to_string(n); return false; }
    if (m < 0 || m > 65535) { err = "number of edges must be in 0..65535, got " + std::to_string(m); return false; }

    g = Graph();
    g.n = static_cast<int>(n);
    g.adj.assign(n, 0);
    std::set<std::pair<int, int>> seen;
    for (long long k = 0; k < m; k++) {
        long long a, b;
        if (!(in >> a >> b)) { err = "malformed edge list in " + path; return false; }
        if (a < 0 || a >= n || b < 0 || b >= n) {
            err = "edge (" + std::to_string(a) + ", " + std::to_string(b) + ") out of range in " + path;
            return false;
        }
        if (a == b) { err = "self loop on node " + std::to_string(a) + " in " + path; return false; }
        if (!seen.insert({std::min(a, b), std::max(a, b)}).second) {
            err = "repeated edge (" + std::to_string(a) + ", " + std::to_string(b) + ") in " + path;
            return false;
        }
        g.edges.push_back({static_cast<int>(a), static_cast<int>(b)});
        g.adj[a] |= 1ULL << b;
        g.adj[b] |= 1ULL << a;
    }
    return true;
}

// ---------------------------------------------------------------- cut table

std::vector<uint16_t> buildCutTable(const Graph& g, int& maxCut) {
    const ull N = 1ULL << g.n;
    const int n = g.n;
    const ull* adj = g.adj.data();
    std::vector<uint16_t> cut(N);
    int best = 0;
    #pragma omp parallel for reduction(max:best) schedule(static)
    for (ull x = 0; x < N; x++) {
        // Each crossing edge is counted once, from its endpoint on the 1 side.
        int c = 0;
        for (int i = 0; i < n; i++) {
            if ((x >> i) & 1) c += __builtin_popcountll(adj[i] & ~x);
        }
        cut[x] = static_cast<uint16_t>(c);
        best = std::max(best, c);
    }
    maxCut = best;
    return cut;
}

// ---------------------------------------------------------------- cost layer

void applyCostLayer(State& state, const std::vector<uint16_t>& cut, int m, double gamma) {
    const ull N = state.numAmpTotal;
    const double norm = 1.0 / std::sqrt(static_cast<double>(N));
    std::vector<Complex> tab(m + 1);
    for (int c = 0; c <= m; c++) {
        tab[c].real =  norm * std::cos(gamma * c);
        tab[c].imag = -norm * std::sin(gamma * c);
    }
    const Complex* t = tab.data();
    const uint16_t* cx = cut.data();
    Complex* psi = state.stateVec;
    #pragma omp parallel for schedule(static)
    for (ull x = 0; x < N; x++) psi[x] = t[cx[x]];
}

// ---------------------------------------------------------------- mixing layer

static inline void rxPair(Complex& up, Complex& lo, double c, double s) {
    const Complex u = up, l = lo;
    // up' = c*up - i*s*lo,  lo' = -i*s*up + c*lo
    up.real =  u.real * c + l.imag * s;
    up.imag =  u.imag * c - l.real * s;
    lo.real =  u.imag * s + l.real * c;
    lo.imag = -u.real * s + l.imag * c;
}

// Apply Rx on the G consecutive qubits t0 .. t0+G-1 in one pass: each
// iteration loads the 2^G amplitudes that differ only in those bits, applies
// the G butterfly stages in registers and stores them back.
template <int G>
static void mixHighGroup(State& state, int t0, double c, double s) {
    const ull N = state.numAmpTotal;
    const ull lowMask = (1ULL << t0) - 1;
    Complex* psi = state.stateVec;
    #pragma omp parallel for schedule(static)
    for (ull i = 0; i < N >> G; i++) {
        const ull base = ((i >> t0) << (t0 + G)) | (i & lowMask);
        Complex a[1 << G];
        for (int j = 0; j < (1 << G); j++) a[j] = psi[base + (ull(j) << t0)];
        for (int q = 0; q < G; q++)
            for (int j = 0; j < (1 << G); j++)
                if (!((j >> q) & 1)) rxPair(a[j], a[j | (1 << q)], c, s);
        for (int j = 0; j < (1 << G); j++) psi[base + (ull(j) << t0)] = a[j];
    }
}

void applyMixingLayer(State& state, int n, double beta, int blockBits) {
    const double c = std::cos(beta);
    const double s = std::sin(beta);
    const ull N = state.numAmpTotal;

    // Low qubits: blocks of 2^k amplitudes (2^14 * 16 B = 256 KB, fits in L2).
    // Shrink the block if there would be fewer blocks than threads.
    int k = std::min(blockBits, n);
    const ull threads = static_cast<ull>(omp_get_max_threads());
    while (k > 10 && (N >> k) < threads) k--;

    const ull blockSize = 1ULL << k;
    Complex* psi = state.stateVec;
    #pragma omp parallel for schedule(static)
    for (ull b = 0; b < (N >> k); b++) {
        Complex* blk = psi + b * blockSize;
        for (int t = 0; t < k; t++) {
            const ull stride = 1ULL << t;
            for (ull base = 0; base < blockSize; base += 2 * stride)
                for (ull o = 0; o < stride; o++)
                    rxPair(blk[base + o], blk[base + o + stride], c, s);
        }
    }

    // High qubits: groups of up to three, one pass per group. The per-amplitude
    // operation order (qubit 0, 1, ..., n-1) is the same as gate-by-gate Rx.
    for (int t = k; t < n; t += 3) {
        switch (std::min(3, n - t)) {
            case 1: mixHighGroup<1>(state, t, c, s); break;
            case 2: mixHighGroup<2>(state, t, c, s); break;
            default: mixHighGroup<3>(state, t, c, s); break;
        }
    }
}

// ---------------------------------------------------------------- measurement

Measurement measure(const State& state, const std::vector<uint16_t>& cut, int maxCut,
                    int shots, std::mt19937_64& rng, std::vector<ull>* samples) {
    const ull N = state.numAmpTotal;
    const Complex* psi = state.stateVec;
    const uint16_t* cx = cut.data();

    // A fixed number of chunks (independent of the thread count) keeps the
    // floating-point summation order, and therefore the output, identical for
    // any OMP_NUM_THREADS. 96 divides evenly among 6 or 12 threads.
    const ull C = std::min<ull>(N, 96);
    std::vector<double> chunkP(C), chunkE(C), chunkOpt(C);
    auto chunkBegin = [N, C](ull ch) { return N * ch / C; };

    #pragma omp parallel for schedule(static)
    for (ull ch = 0; ch < C; ch++) {
        double p = 0.0, e = 0.0, o = 0.0;
        const ull b = chunkBegin(ch), end = chunkBegin(ch + 1);
        #pragma omp simd reduction(+:p,e,o)
        for (ull x = b; x < end; x++) {
            const double px = psi[x].real * psi[x].real + psi[x].imag * psi[x].imag;
            p += px;
            e += px * cx[x];
            o += (cx[x] == maxCut) ? px : 0.0;
        }
        chunkP[ch] = p; chunkE[ch] = e; chunkOpt[ch] = o;
    }

    Measurement r;
    double total = 0.0;
    for (ull ch = 0; ch < C; ch++) {
        total += chunkP[ch];
        r.expectedCut += chunkE[ch];
        r.pOpt += chunkOpt[ch];
    }
    // Normalise away the O(1e-15) drift of the state norm.
    r.expectedCut /= total;
    r.pOpt /= total;

    if (shots <= 0) {
        if (samples) samples->clear();
        return r;
    }

    // Exact sampling by inverse CDF: sort the uniform draws, find which chunk
    // each falls into from the chunk prefix sums, then walk each chunk once.
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<double> u(shots);
    for (int k = 0; k < shots; k++) u[k] = dist(rng) * total;
    std::sort(u.begin(), u.end());

    std::vector<double> prefix(C + 1, 0.0);
    for (ull ch = 0; ch < C; ch++) prefix[ch + 1] = prefix[ch] + chunkP[ch];

    std::vector<ull> picked(shots);
    #pragma omp parallel for schedule(dynamic)
    for (ull ch = 0; ch < C; ch++) {
        auto lo = std::lower_bound(u.begin(), u.end(), prefix[ch]) - u.begin();
        auto hi = (ch + 1 == C) ? static_cast<long>(shots)
                                : std::lower_bound(u.begin(), u.end(), prefix[ch + 1]) - u.begin();
        if (lo >= hi) continue;
        ull x = chunkBegin(ch);
        const ull last = chunkBegin(ch + 1) - 1;
        double cum = prefix[ch];
        for (long k = lo; k < hi; k++) {
            while (x < last) {
                const double px = psi[x].real * psi[x].real + psi[x].imag * psi[x].imag;
                // px > 0 guard: never land on a zero-probability state when
                // rounding leaves cum slightly above u[k]
                if (px > 0.0 && cum + px > u[k]) break;
                cum += px;
                x++;
            }
            picked[k] = x;
        }
    }

    int best = 0;
    for (ull x : picked) best = std::max<int>(best, cx[x]);
    r.sampledMaxCut = best;
    if (samples) samples->swap(picked);
    return r;
}
