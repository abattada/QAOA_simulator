#include "circuit.hpp"

ull insert_bit_0(ull task, const char targ) {
    ull mask = (1ULL << targ) - 1;
    return ((task >> targ) << (targ+1)) | (task & mask);
}

ull insert_bit_1(ull task, const char targ) {
    ull mask = (1ULL << targ) - 1;
    return ((task >> targ) << (targ+1)) | (task & mask) | (1ULL << targ);
}

void rxGate(State& state, const int targQubit, const double& angle) {
    double sinAngle, cosAngle;
    sincosf64(angle/2, &sinAngle, &cosAngle);
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull up_off = insert_bit_0(i, targQubit);
        ull lo_off = up_off + (1ull << targQubit);

        Complex up = state.stateVec[up_off];
        Complex lo = state.stateVec[lo_off];

        state.stateVec[up_off].real =  up.real * cosAngle + lo.imag * sinAngle;
        state.stateVec[up_off].imag =  up.imag * cosAngle - lo.real * sinAngle;
        state.stateVec[lo_off].real =  up.imag * sinAngle + lo.real * cosAngle;
        state.stateVec[lo_off].imag = -up.real * sinAngle + lo.imag * cosAngle;
    }
}

void hGate(State& state, const int targQubit) {
    const double coeff = 1 / sqrt(2.f);
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull up_off = insert_bit_0(i, targQubit);
        ull lo_off = up_off + (1ull << targQubit);

        Complex up = state.stateVec[up_off];
        Complex lo = state.stateVec[lo_off];

        state.stateVec[up_off].real = coeff * (up.real + lo.real);
        state.stateVec[up_off].imag = coeff * (up.imag + lo.imag);
        state.stateVec[lo_off].real = coeff * (up.real - lo.real);
        state.stateVec[lo_off].imag = coeff * (up.imag - lo.imag);
    }
}


void rzGate(State& state, const int targQubit, const double& angle) {
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal; i++) {
        bool isTargetOn = (i >> targQubit) & 1;

        // Apply the phase shift only to |1>
        if (isTargetOn) {
            double realPart = state.stateVec[i].real;
            double imagPart = state.stateVec[i].imag;

            state.stateVec[i].real = realPart * cos(angle / 2) - imagPart * sin(angle / 2);
            state.stateVec[i].imag = realPart * sin(angle / 2) + imagPart * cos(angle / 2);
        }
    }
}

void cpGate(State& state, const double angle, const std::vector<int>& targQubits) {

    int q1 = targQubits[0], q2 = targQubits[1];
    for (ull i = 0; i < state.numAmpTotal; i++) {
        bool isQ1On = (i >> q1) & 1; // Check the q1 bit
        bool isQ2On = (i >> q2) & 1; // Check the q2 bit

        if (isQ1On && isQ2On) { // Both qubits in |1>
            double realPart = state.stateVec[i].real;
            double imagPart = state.stateVec[i].imag;

            state.stateVec[i].real = realPart * cos(angle) - imagPart * sin(angle);
            state.stateVec[i].imag = realPart * sin(angle) + imagPart * cos(angle);
        }
    }
}

void cnotGate(State& state, const int controlQubit, const int targetQubit) {
    #pragma omp parallel for
    for (ull i = 0; i < state.numAmpTotal >> 1; i++) {
        ull controlOnOff = insert_bit_1(i, controlQubit);
        ull targetOff = controlOnOff & ~(1ULL << targetQubit);
        ull targetOn = targetOff | (1ULL << targetQubit);

        std::swap(state.stateVec[targetOff], state.stateVec[targetOn]);
    }
}

