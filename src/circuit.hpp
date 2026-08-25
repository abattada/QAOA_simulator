#ifndef _CIRCUIT_HPP_
#define _CIRCUIT_HPP_

#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "state.hpp"


ull insert_bit_0(ull task, const char targ);
ull insert_bit_1(ull task, const char targ);
void rxGate(State& state, const int targQubit, const double& angle);
void hGate(State& state, const int targQubit);
void rzGate(State& state, const int targQubit, const double& angle);
void cnotGate(State& state, const int controlQubit, const int targetQubit);
void cpGate(State& state, const double angle, const std::vector<int>& targQubits);

#endif
