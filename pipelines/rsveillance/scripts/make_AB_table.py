#!/usr/bin/env python3
"""Build one tidy per-sample RSV A/B table from rsveillance summary outputs.

Joins the full sample list with the per-target alignment stats and the final
A/B calls, emitting one row per input sample (including samples that never got
a call) plus a no-call list for re-sequencing.

Usage:
    python3 make_AB_table.py <samples.txt> <results/summary> [nocall_out.txt] > table.tsv

Output columns (TSV, to stdout):
    sample | best_type | call | coverage_pct | mean_depth | verdict

  best_type     genome the sample covered best (RSVA / RSVB / none)
  call          final call from final_calls.txt (RSVA / RSVB / RSVA/B / no-call)
  coverage_pct  best goodcov% (fraction of genome >=10x), as a percentage
  mean_depth    mean depth at that best-covered genome
  verdict       PASS          -> got a confident call
                drop_lowcov   -> mapped to A and/or B but <MINCOV coverage
                no_mash_call  -> too few reads for mash to assign a type
"""
import sys
import csv
import os

MINCOV = 0.8  # matches call_RSVAB.py -c; only used for labelling here

# per-target alignstats columns (no header):
#   0 sample 1 target 2 subfac 3 reads 4 aligned 5 paired
#   6 meandepth 7 goodcov 8 cov 9 gsize
COL_MEANDEPTH, COL_GOODCOV, COL_GSIZE = 6, 7, 9


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make_AB_table.py <samples.txt> <summary_dir> [nocall_out.txt]")
    samples_file = sys.argv[1]
    summary_dir = sys.argv[2]
    nocall_out = sys.argv[3] if len(sys.argv) > 3 else "nocall_samples.txt"

    # 1. all input samples (the full set)
    all_samples = [l.strip() for l in open(samples_file) if l.strip()]

    # 2. per-(sample, target) stats from the per-target alignstats files
    stats = {}  # sample -> [(target, coverage_fraction, mean_depth), ...]
    for tgt in ("RSVA", "RSVB"):
        path = os.path.join(summary_dir, f"{tgt}_alignstats.txt")
        if not os.path.exists(path):
            continue
        for line in open(path):
            c = line.rstrip("\n").split("\t")
            if len(c) <= COL_GSIZE:
                continue
            try:
                meandepth = float(c[COL_MEANDEPTH])
                goodcov = float(c[COL_GOODCOV])
                gsize = float(c[COL_GSIZE])
            except ValueError:
                continue  # skip any stray header/blank line
            covpc = goodcov / gsize if gsize else 0.0
            stats.setdefault(c[0], []).append((c[1], covpc, meandepth))

    # 3. final A/B/mixed calls
    calls = {}
    with open(os.path.join(summary_dir, "final_calls.txt")) as f:
        next(f, None)  # skip header
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2 and c[0]:
                calls[c[0]] = c[1]

    # 4. one tidy row per sample
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(
        ["sample", "best_type", "call", "coverage_pct", "mean_depth", "verdict"]
    )
    nocall = []
    for s in all_samples:
        rows = stats.get(s, [])
        if rows:
            best_type, covpc, depth = max(rows, key=lambda r: r[1])
            covstr, depthstr = f"{covpc * 100:.1f}", f"{depth:.0f}"
        else:
            best_type, covstr, depthstr = "none", "0.0", "0"
        call = calls.get(s, "")
        if call:
            verdict = "PASS"
        else:
            verdict = "drop_lowcov" if rows else "no_mash_call"
            nocall.append(s)
        writer.writerow([s, best_type, call or "no-call", covstr, depthstr, verdict])

    with open(nocall_out, "w") as out:
        out.write("\n".join(nocall) + ("\n" if nocall else ""))

    sys.stderr.write(
        f"{len(all_samples)} samples | {len(all_samples) - len(nocall)} called | "
        f"{len(nocall)} no-call -> {nocall_out}\n"
    )


if __name__ == "__main__":
    main()
