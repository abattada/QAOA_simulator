"""Brute-force p=1 QAOA reference in pure Python (no numpy).

Builds the state amplitude by amplitude with the same conventions as the C++
simulator (cost layer e^{-i gamma cut(x)} on |+>^n, then Rx(2 beta) on every
qubit) and prints <C> and the probability of measuring a maximum cut.
Intended for small graphs (n <= 12) to cross-check the simulator output.

Usage:
    python3 tests/ref_qaoa.py graph/qml_bi_4 <gamma> <beta>
    python3 tests/ref_qaoa.py graph/qml_bi_4 --grid <ig> <ib> <steps>
"""
import argparse
import cmath
import math


def load_graph(path):
    tokens = open(path).read().split()
    n, m = int(tokens[0]), int(tokens[1])
    edges = [(int(tokens[2 + 2 * k]), int(tokens[3 + 2 * k])) for k in range(m)]
    return n, edges


def qaoa(n, edges, gamma, beta):
    N = 1 << n
    cut = [sum(((x >> i) ^ (x >> j)) & 1 for i, j in edges) for x in range(N)]
    psi = [cmath.exp(-1j * gamma * cut[x]) / math.sqrt(N) for x in range(N)]
    c, s = math.cos(beta), math.sin(beta)
    for q in range(n):
        for x in range(N):
            if not (x >> q) & 1:
                y = x | (1 << q)
                a, b = psi[x], psi[y]
                psi[x], psi[y] = c * a - 1j * s * b, -1j * s * a + c * b
    probs = [abs(a) ** 2 for a in psi]
    max_cut = max(cut)
    expected = sum(p * k for p, k in zip(probs, cut))
    p_opt = sum(p for p, k in zip(probs, cut) if k == max_cut)
    return expected, p_opt, max_cut


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("graph")
    parser.add_argument("values", nargs="+", type=float,
                        help="gamma beta, or ig ib steps with --grid")
    parser.add_argument("--grid", action="store_true",
                        help="interpret values as grid indices ig ib steps (same formula as the C++ driver)")
    args = parser.parse_args()

    if args.grid:
        ig, ib, steps = (int(v) for v in args.values)
        gamma = 0.0 if steps == 1 else math.pi * ig / (steps - 1)
        beta = 0.0 if steps == 1 else math.pi * ib / (steps - 1)
    else:
        gamma, beta = args.values

    n, edges = load_graph(args.graph)
    expected, p_opt, max_cut = qaoa(n, edges, gamma, beta)
    print(f"gamma={gamma:.12g} beta={beta:.12g} maxcut={max_cut} "
          f"expected_cut={expected:.12g} p_opt={p_opt:.12g}")


if __name__ == "__main__":
    main()
