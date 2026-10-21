#include "circuit.hpp"

#include <cmath>
#include <utility>

ull insert_bit_0(ull task, int targ) {
    ull mask = (1ULL << targ) - 1;
    return ((task >> targ) << (targ + 1)) | (task & mask);
}

ull insert_bit_1(ull task, int targ) {
    return insert_bit_0(task, targ) | (1ULL << targ);
}

void hGate(State& state, int targQubit) {
    const double coeff = 1.0 / std::sqrt(2.0);
    const ull stride = 1ULL << targQubit;
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull up_off = insert_bit_0(i, targQubit);
        ull lo_off = up_off + stride;

        Complex up = state.stateVec[up_off];
        Complex lo = state.stateVec[lo_off];

        state.stateVec[up_off].real = coeff * (up.real + lo.real);
        state.stateVec[up_off].imag = coeff * (up.imag + lo.imag);
        state.stateVec[lo_off].real = coeff * (up.real - lo.real);
        state.stateVec[lo_off].imag = coeff * (up.imag - lo.imag);
    }
}

void rxGate(State& state, int targQubit, double angle) {
    const double c = std::cos(angle / 2);
    const double s = std::sin(angle / 2);
    const ull stride = 1ULL << targQubit;
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull up_off = insert_bit_0(i, targQubit);
        ull lo_off = up_off + stride;

        Complex up = state.stateVec[up_off];
        Complex lo = state.stateVec[lo_off];

        // up' = c*up - i*s*lo,  lo' = -i*s*up + c*lo
        state.stateVec[up_off].real =  up.real * c + lo.imag * s;
        state.stateVec[up_off].imag =  up.imag * c - lo.real * s;
        state.stateVec[lo_off].real =  up.imag * s + lo.real * c;
        state.stateVec[lo_off].imag = -up.real * s + lo.imag * c;
    }
}

void phaseGate(State& state, int targQubit, double theta) {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    // Diagonal gate: only the half of the state with the target bit set changes.
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull off = insert_bit_1(i, targQubit);
        double re = state.stateVec[off].real;
        double im = state.stateVec[off].imag;
        state.stateVec[off].real = re * c - im * s;
        state.stateVec[off].imag = re * s + im * c;
    }
}

void cnotGate(State& state, int controlQubit, int targetQubit) {
    // Enumerate indices with control = 1 and target = 0; each one names a
    // distinct (target=0, target=1) pair, so every pair is swapped exactly once.
    const int lo = controlQubit < targetQubit ? controlQubit : targetQubit;
    const int hi = controlQubit < targetQubit ? targetQubit : controlQubit;
    const ull controlBit = 1ULL << controlQubit;
    const ull targetBit = 1ULL << targetQubit;
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 2; i++) {
        ull off0 = insert_bit_0(insert_bit_0(i, lo), hi) | controlBit;
        std::swap(state.stateVec[off0], state.stateVec[off0 | targetBit]);
    }
}
