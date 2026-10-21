"""Convert a raw grid-search log into a CSV file next to it.

Columns: Gamma, Beta, Expected Cut, P(opt) and, when the run drew samples,
Sampled Max Cut. Rows are sorted by (Gamma, Beta).

Usage:
    python3 scripts/transform_to_csv.py results/qml_bi_10/1000_20_0
"""
import argparse
import csv
import re

NUM = r"([0-9.eE+-]+)"
LINE_RE = re.compile(
    rf"Gamma:\s*{NUM},\s*Beta:\s*{NUM},\s*Expected Cut:\s*{NUM},\s*P\(opt\):\s*{NUM}"
    rf"(?:,\s*Sampled Max Cut:\s*{NUM})?")

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("file_name", help="raw simulator output file")
args = parser.parse_args()

with open(args.file_name, "r", encoding="utf-8") as f:
    rows = [m.groups() for m in map(LINE_RE.search, f) if m]

if not rows:
    raise SystemExit(f"no grid lines found in {args.file_name}")

sampled = rows[0][4] is not None
header = ["Gamma", "Beta", "Expected Cut", "P(opt)"] + (["Sampled Max Cut"] if sampled else [])
rows = [r if sampled else r[:4] for r in rows]
rows.sort(key=lambda r: (float(r[0]), float(r[1])))

with open(args.file_name + ".csv", "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    writer.writerows(rows)

print(f"wrote {args.file_name}.csv ({len(rows)} rows)")
