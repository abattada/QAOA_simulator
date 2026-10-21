#ifndef _STATE_HPP_
#define _STATE_HPP_

typedef unsigned long long ull;

typedef struct {
    double real;
    double imag;
} Complex;

typedef struct {
    Complex* stateVec;
    ull numAmpTotal;
} State;

// Allocate 2^numQubits amplitudes (64-byte aligned) without initialising them.
// Returns false if the allocation fails.
bool allocState(State& state, int numQubits);

// Allocate and initialise to |0...0>. Returns false if the allocation fails.
bool createState(State& state, int numQubits);

void printState(const State& state, int num);
void destroyState(State& state);

#endif
