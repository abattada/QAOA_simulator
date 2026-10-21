# QAOA_simulator

**A from-scratch, OpenMP-parallel state-vector simulator for the Quantum Approximate Optimization Algorithm (QAOA) on MaxCut, with fused kernels that run a full p=1 grid search at 24 qubits in about a minute.**

![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)
![OpenMP](https://img.shields.io/badge/OpenMP-kernel--level-2a78d6)
![Python 3](https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

No quantum SDK and no external libraries: every kernel is hand-written C++. Instead of applying the textbook circuit gate by gate, the simulator exploits the structure of p=1 QAOA. The cost layer is diagonal, so it is written in a single pass from a precomputed cut table, and the mixing layer is cache-blocked. Each `(γ, β)` grid point passes over the state vector 7 to 8 times, where a gate-by-gate simulation of the same circuit needs about 500 passes at 24 qubits. The expected cut, the probability of the optimal cut, and bitstring samples are all computed exactly.

![Runtime benchmark vs QuEST](docs/figures/runtime_benchmark.png)

## Highlights

- **One-pass cost layer.** The cut value of every bitstring is precomputed once per graph with `popcount` over neighbour bitmasks (2 bytes per bitstring). Each grid point then writes `ψ[x] = e^{-iγ·cut(x)} / √2ⁿ` directly, which replaces the Hadamard layer and all `3m` CNOT/phase gates.
- **Cache-blocked mixing layer.** `Rx(2β)` on the 14 lowest qubits is applied inside 256 KB blocks that stay in L2, in one pass over memory. The remaining qubits are processed three at a time, one pass per group. At 24 qubits the whole mixing layer costs 5 passes instead of 24.
- **Exact measurement.** `⟨C⟩ = Σ |ψₓ|² cut(x)` and `P(opt)`, the probability of measuring a maximum cut, come out of one pass. Optional shots are drawn from the exact joint distribution by inverse-CDF sampling, with a deterministic seed per grid point.
- **Kernel-level OpenMP.** One state vector is reused for the whole grid, and the parallelism lives inside each kernel. Output is bit-identical for any thread count.
- **Tested against the gate-by-gate circuit.** `make test` checks the fused kernels against `H`, `CNOT`, phase and `Rx` gate kernels, and `tests/ref_qaoa.py` is an independent pure-Python brute force.
- **Complete pipeline.** Dataset download (PennyLane HamLib), C++ simulation, CSV export, heatmaps and benchmark figures, each one command.

## The algorithm

QAOA (Farhi, Goldstone & Gutmann, 2014) is a variational quantum algorithm for combinatorial optimization. For MaxCut on a graph `G = (V, E)` with `m` edges, the p=1 circuit prepares `|+⟩ⁿ`, applies the **cost layer** `e^{-iγ H_C}` and the **mixing layer** `e^{-iβ H_B}`, and measures in the computational basis. The measured bitstring is a bipartition of the graph, and its cut value is the number of edges crossing it. The classical outer loop tunes `(γ, β)` to maximize the expected cut `⟨C⟩`. Here that loop is an exhaustive grid search, which doubles as a map of the optimization landscape:

<p align="center">
  <img src="docs/figures/heatmaps/qml_bi_10.png" width="560" alt="Grid-search landscape for qml_bi_10 (20 qubits): expected cut over (gamma, beta), best point ringed">
</p>

Landscapes for all nine benchmark instances (8 to 24 qubits) are in [docs/figures/heatmaps/](docs/figures/heatmaps/):

| 8 qubits (`qml_bi_4`) | 16 qubits (`qml_bi_8`) | 24 qubits (`qml_bi_12`) |
|---|---|---|
| ![heatmap 8 qubits](docs/figures/heatmaps/qml_bi_4.png) | ![heatmap 16 qubits](docs/figures/heatmaps/qml_bi_8.png) | ![heatmap 24 qubits](docs/figures/heatmaps/qml_bi_12.png) |

**Conventions.** The cost layer multiplies each basis state by `e^{-iγ·cut(x)}`. In gates, that is `H` on every qubit followed by `CNOT(i,j) · P_j(-γ) · CNOT(i,j)` for each edge, where `P(θ) = diag(1, e^{iθ})`. The mixing layer is `Rx(2β) = [[cos β, -i sin β], [-i sin β, cos β]]` on every qubit.

**Why the simulator never runs the H layer or the CNOTs.** The state after `H^⊗n` is the constant `1/√2ⁿ`, and every gate of the cost layer is diagonal in the computational basis. Their product is the single diagonal phase `e^{-iγ·cut(x)}`, which `applyCostLayer` writes in one pass. The test suite confirms the fused state matches the gate-by-gate state to about 1e-16.

**Grid range.** Cut values are integers, so the state is periodic in `γ` with period `2π`. It is also unchanged up to a global phase when `β` moves by `π/2`, because the global bit flip `X^⊗n` preserves every cut. Complex conjugation maps `(γ, β)` to `(2π-γ, π-β)` with the same probabilities. Together these mean the grid `(γ, β) ∈ [0, π]²` contains the global optimum, and the `β` axis shows the landscape twice. On the complete bipartite benchmarks there is one more mirror symmetry, so the best `⟨C⟩` appears at two grid points with identical values, and the log reports the first.

## Results

Full 20×20 grid with 1000 exact shots per point on the HamLib complete bipartite graphs `K_{N,N}`, whose maximum cut is `N²`.

| qubits | graph | max cut | best `⟨C⟩` | ratio `⟨C⟩ / max` | `P(opt)` at best point | best sampled cut |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 8 | `qml_bi_4` | 16 | 10.58 | 0.661 | 0.231 | 16 |
| 10 | `qml_bi_5` | 25 | 16.05 | 0.642 | 0.175 | 25 |
| 12 | `qml_bi_6` | 36 | 22.49 | 0.625 | 0.124 | 36 |
| 14 | `qml_bi_7` | 49 | 30.18 | 0.616 | 0.082 | 49 |
| 16 | `qml_bi_8` | 64 | 39.01 | 0.610 | 0.065 | 64 |
| 18 | `qml_bi_9` | 81 | 48.89 | 0.604 | 0.050 | 81 |
| 20 | `qml_bi_10` | 100 | 59.80 | 0.598 | 0.036 | 100 |
| 22 | `qml_bi_11` | 121 | 71.72 | 0.593 | 0.025 | 121 |
| 24 | `qml_bi_12` | 144 | 84.62 | 0.588 | 0.017 | 144 |

The approximation ratio of p=1 QAOA decreases slowly with size on these graphs. The probability of an optimal cut at the best point stays far above the uniform baseline of `2/2ⁿ`, so 1000 shots per point recover the optimum on every instance. The raw logs and CSVs are in [results/](results/).

### What is exact and what is approximate

- **Approximate by design.** p=1 QAOA is a heuristic: it yields a good cut on average, not a guaranteed optimum. The grid search samples the landscape at a fixed resolution instead of running an optimizer.
- **Exact.** The state evolution, `⟨C⟩`, `P(opt)` and the sampled bitstrings are exact up to floating-point rounding. The true maximum cut is computed by brute force from the cut table.

## Benchmarks vs QuEST

**Setup.** Both simulators ran a full 20×20 `(γ, β)` grid search of p=1 QAOA on the same graphs, timed with `time(1)` on a Ryzen 5 3600X (6 cores, 12 threads, 16 MB L3). This simulator ran with 6 threads pinned one per core. The QuEST numbers are an earlier record from 12 threads.

| qubits | edges | **this work** (wall) | QuEST (wall) | **this work** (CPU) | QuEST (CPU) | **this work** (peak RSS) |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 20 | 100 | **2.4 s** | 1m 09s | **14 s** | 13m 36s | 22 MB |
| 22 | 121 | **16 s** | 11m 30s | **1m 36s** | 135m 21s | 76 MB |
| 24 | 144 | **1m 14s** | 66m 39s | **7m 22s** | 709m 24s | 292 MB |

**Read this as "same task, different algorithm".** QuEST applied the circuit gate by gate: `2m` CNOTs, `m` Z rotations and `2n` single-qubit gates per grid point, each a full pass over the state. This simulator uses the structure of p=1 QAOA to do the same work in 7 to 8 passes. The QuEST driver is not in this repository. Its grid outputs in [benchmarks/quest_results/](benchmarks/quest_results/) show a maximum-cut distribution very close to the previous, flawed version of this simulator, so it likely used the same independent-marginal sampler described below. Its timing is still a fair measure of gate-by-gate simulation cost.

**Memory.** Peak memory is one state vector (16 bytes per amplitude) plus the cut table (2 bytes per bitstring). The QuEST readings come from [benchmarks/memory_screenshots/](benchmarks/memory_screenshots/). The `mysim_*` screenshots in that folder were taken with the previous version, which kept about 12 state vectors alive, and no longer apply.

![Memory benchmark](docs/figures/memory_benchmark.png)

**Threads.** The kernels are limited by memory bandwidth, not arithmetic, so SMT siblings add nothing. At 22 qubits, 12 threads on this 6-core CPU took 16.8 s against 16.1 s for 6 threads. The benchmark script therefore uses one thread per physical core.

## Correctness note

Versions of this repository before October 2026 did not simulate QAOA. The CNOT kernel swapped every amplitude pair twice, which made it an identity when run serially and a data race in parallel. The state therefore never became entangled, and the cost layer reduced to single-qubit phases. The measurement step then sampled each qubit independently from its marginal probability. For a correct MaxCut QAOA state that sampler draws uniform random bitstrings, because global bit-flip symmetry makes every single-qubit marginal exactly 1/2. The earlier claim that the search "found the optimal cut on every instance" came from these two effects together on bipartite graphs, not from QAOA.

The current code fixes the CNOT kernel and replaces the sampler with exact computation and exact sampling. `make test` checks that the fused kernels reproduce the gate-by-gate circuit, that every single-qubit marginal is 1/2, and that sampled cut values follow `|ψ|²`. The C++ output for `qml_bi_4` also matches the independent Python brute force in [tests/ref_qaoa.py](tests/ref_qaoa.py) to 12 significant digits.

## Quick start

```bash
# Build and test (requires g++ with OpenMP)
make
make test

# Run: ./qaoa_main <graph_name> <shots> <steps> [seed]
#   <graph_name>  file in ./graph/
#   <shots>       exact samples per grid point (0 = expectation values only)
#   <steps>       grid resolution: steps × steps points over (γ, β) ∈ [0, π]²
#   [seed]        base RNG seed (default 1); each grid point derives its own
OMP_NUM_THREADS=6 OMP_PROC_BIND=true OMP_PLACES=cores ./qaoa_main qml_bi_10 1000 20
```

Set `OMP_NUM_THREADS` to the number of physical cores. Each grid point writes one line:

```
Gamma: 0.330694, Beta: 1.98416, Expected Cut: 59.800105442, P(opt): 0.0363475404588, Sampled Max Cut: 100
```

The full grid log is written to `results/qml_bi_10/<shots>_<steps>_<run>`, where the run index auto-increments so existing logs are never overwritten. A summary goes to the end of the log and to the terminal, which also shows the thread count and wall time:

```
=========================================
 QAOA (p=1) with 20 qubits.
 Input graph: qml_bi_10 (edges=100)
 Grid search steps = 20, shots = 1000, seed = 1
-----------------------------------------
 True max cut      = 100
 Best expected cut = 59.8001
 Approx ratio      = 0.598001
 P(opt) at best    = 0.0363475
 Best sampled cut  = 100
 Best gamma        = 2.8109
 Best beta         = 1.15743
=========================================
```

The same seed reproduces a log byte for byte, whatever the thread count.

### Benchmarks, plots and CSV

```bash
scripts/run_benchmarks.sh            # 4–24 qubit sweep: results/, benchmarks/my_*_record.txt

python3 -m venv .venv                # needs the python3-venv package
.venv/bin/pip install matplotlib numpy
scripts/make_figures.sh              # all README figures from results/ and benchmarks/

python3 scripts/transform_to_csv.py results/qml_bi_10/1000_20_0       # grid log to CSV
python3 tests/ref_qaoa.py graph/qml_bi_4 --grid 5 3 20                # brute-force check of one grid point
```

### Graphs

Benchmark graphs `K_{N,N}` for 2 to 50 qubits ship in [graph/](graph/); the simulator accepts up to 32 nodes, memory permitting. The input format is plain text: number of nodes, number of edges, then one `u v` pair per line. Graphs must be simple, with no self loops or repeated edges. Two generators are included:

```bash
python3 scripts/generate_graph.py       # random n-node, m-edge graph (G(n, m) model)
pip install pennylane
python3 scripts/download_maxcut.py      # fetch HamLib MaxCut graphs via PennyLane
```

## Repository layout

```
QAOA_simulator/
├── src/
│   ├── state.{hpp,cpp}       # 2^n amplitude array (64-byte aligned)
│   ├── circuit.{hpp,cpp}     # gate kernels: H, Rx, phase, CNOT (reference path)
│   ├── qaoa.{hpp,cpp}        # graph loader, cut table, fused cost/mixing layers, exact measurement
│   └── main.cpp              # p=1 QAOA (γ, β) grid search driver
├── tests/
│   ├── test_qaoa.cpp         # fused kernels vs gate-by-gate circuit (make test)
│   └── ref_qaoa.py           # independent pure-Python brute force
├── scripts/                # benchmark sweep, figures, CSV, graph download/generation
├── graph/                  # bundled benchmark graphs (K_{N,N}, 2–50 qubits)
├── results/                # grid-search logs + CSV (8–24 qubits)
├── benchmarks/             # timing and memory records, QuEST outputs and screenshots
└── docs/figures/           # figures used in this README
```

## References

- E. Farhi, J. Goldstone, S. Gutmann, *A Quantum Approximate Optimization Algorithm*, [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) (2014).
- T. Jones, A. Brown, I. Bush, S. C. Benjamin, *QuEST and High Performance Simulation of Quantum Computers*, [Sci. Rep. 9, 10736](https://www.nature.com/articles/s41598-019-47174-9) (2019).
- N. P. D. Sawaya et al., *HamLib: A library of Hamiltonians for benchmarking quantum algorithms and hardware*, [arXiv:2306.13126](https://arxiv.org/abs/2306.13126) (2023). MaxCut instances loaded via [PennyLane Datasets](https://pennylane.ai/datasets).

---

*Originally developed as the final project for a Quantum Design Automation course (Dec 2024).*
