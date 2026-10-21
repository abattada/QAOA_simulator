#!/usr/bin/env bash
# Time the full 20x20 grid search on the bundled K_{N,N} graphs (4..24 qubits)
# and rewrite benchmarks/my_time_record.txt and benchmarks/my_memory_record.txt.
# Each run also writes its grid log to results/qml_bi_<k>/<shots>_20_<run>.
#
# Usage (from the repo root, after `make`):
#     scripts/run_benchmarks.sh [shots]        # default 1000 shots
#
# Threads: one per physical core. The kernels are memory-bandwidth bound, so
# SMT siblings only add contention; override with OMP_NUM_THREADS if needed.
set -euo pipefail

SHOTS=${1:-1000}
STEPS=20
CORES=$(lscpu -p=CORE 2>/dev/null | grep -v '^#' | sort -u | wc -l)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$CORES}
export OMP_PROC_BIND=${OMP_PROC_BIND:-true}
export OMP_PLACES=${OMP_PLACES:-cores}

TIME_OUT=benchmarks/my_time_record.txt
MEM_OUT=benchmarks/my_memory_record.txt
TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

fmt() {  # seconds -> "XmY.YYYs", the format of the QuEST record
    awk -v t="$1" 'BEGIN { m = int(t / 60); printf "%dm%.3fs", m, t - 60 * m }'
}

echo "threads: $OMP_NUM_THREADS, shots: $SHOTS, steps: $STEPS"
time_body=""
mem_body=""
for k in $(seq 2 12); do
    /usr/bin/time -f "%e %U %M" -o "$TMP" ./qaoa_main "qml_bi_$k" "$SHOTS" "$STEPS" 2>/dev/null
    read -r real user rss_kb < "$TMP"
    printf "  %2d qubits: real %8.2f s  user %8.2f s  maxRSS %8.1f MB\n" \
        $((2 * k)) "$real" "$user" "$(awk -v r="$rss_kb" 'BEGIN { print r / 1024 }')"
    sep=$([ "$k" -lt 12 ] && echo "," || echo "")
    time_body+="    $k:\n    {\n        real:\t$(fmt "$real"),\n        user:\t$(fmt "$user")\n    }$sep\n"
    mem_body+="    $k:\n    {\n        maxrss_mb:\t$(awk -v r="$rss_kb" 'BEGIN { printf "%.1f", r / 1024 }')\n    }$sep\n"
done

printf "{\n%b}\n" "$time_body" > "$TIME_OUT"
printf "{\n%b}\n" "$mem_body" > "$MEM_OUT"
echo "wrote $TIME_OUT and $MEM_OUT"
