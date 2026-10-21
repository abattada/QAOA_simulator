#include "state.hpp"

#include <cstdio>
#include <cstdlib>

bool allocState(State& state, int numQubits) {
    state.numAmpTotal = 1ULL << numQubits;
    // aligned_alloc needs a size that is a multiple of the alignment
    size_t bytes = state.numAmpTotal * sizeof(Complex);
    bytes = (bytes + 63) & ~size_t(63);
    state.stateVec = static_cast<Complex*>(std::aligned_alloc(64, bytes));
    return state.stateVec != nullptr;
}

bool createState(State& state, int numQubits) {
    if (!allocState(state, numQubits)) return false;
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal; i++) {
        state.stateVec[i].real = 0.0;
        state.stateVec[i].imag = 0.0;
    }
    state.stateVec[0].real = 1.0;
    return true;
}

void printState(const State& state, int num) {
    // Debug helper: print to stderr so it bypasses any stdout redirection
    for (int i = 0; i < num; i++) {
        fprintf(stderr, "Statevec idx: %d, real: %f, imag: %f\n", i,
                state.stateVec[i].real, state.stateVec[i].imag);
    }
}

void destroyState(State& state) {
    std::free(state.stateVec);
    state.stateVec = nullptr;
    state.numAmpTotal = 0;
}
