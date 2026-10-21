#ifndef _CIRCUIT_HPP_
#define _CIRCUIT_HPP_

#include "state.hpp"

// Generic gate kernels. The QAOA driver does not call these (it uses the fused
// kernels in qaoa.hpp); they form the gate-by-gate reference path that the
// tests compare the fused kernels against.

// Insert a 0 (resp. 1) bit at position `targ` of `task`, shifting the higher
// bits up by one. Enumerating task = 0 .. 2^(n-1)-1 visits every index whose
// bit `targ` is 0 (resp. 1) exactly once.
ull insert_bit_0(ull task, int targ);
ull insert_bit_1(ull task, int targ);

void hGate(State& state, int targQubit);

// Rx(angle) = [[cos(angle/2), -i sin(angle/2)], [-i sin(angle/2), cos(angle/2)]]
void rxGate(State& state, int targQubit, double angle);

// Phase gate P(theta) = diag(1, e^{i theta})
void phaseGate(State& state, int targQubit, double theta);

void cnotGate(State& state, int controlQubit, int targetQubit);

#endif
