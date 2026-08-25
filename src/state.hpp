#ifndef _STATE_HPP_
#define _STATE_HPP_

#include <cstring>
#include <cstdio>
#include <vector>
#include <iostream>

typedef unsigned long long ull;

typedef struct {
    double real;
    double imag;
} Complex;

typedef struct {
    Complex* stateVec;
    ull numAmpTotal;
} State;

void createState(State& state, const int& numQubits);
void printState(const State& state, int num);
void destroyState(const State& state);
#endif
