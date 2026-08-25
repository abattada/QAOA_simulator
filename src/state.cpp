#include "state.hpp"

void createState(State& state, const int& numQubits) {
    state.numAmpTotal = 1ull << numQubits;
    state.stateVec = (Complex* )calloc(state.numAmpTotal, sizeof(Complex));
    memset(state.stateVec, 0, state.numAmpTotal * sizeof(Complex));
    state.stateVec[0].real = 1.0f;
    return;
}

void printState(const State& state, int num) {
    // Debug helper: print to stderr so it bypasses any stdout redirection
    for (int i = 0; i < num; i++) {
        fprintf(stderr, "Statevec idx: %d, real: %f, imag: %f\n", i, state.stateVec[i].real, state.stateVec[i].imag);
    }
    return;
}

void destroyState(const State& state) {
    free(state.stateVec);
    return;
}
