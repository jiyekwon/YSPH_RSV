#!/usr/bin/env python3
"""Build one tidy per-sample hMPV table from the pipeline's summary outputs.

Joins the full sample list with the per-target alignment stats and the final
subgroup calls, emitting one row per input sample (including samples that never
got a call) plus a no-call list for re-sequencing.

Usage:
    python3 make_hMPV_table.py <samples.txt> <results/summary> [nocall_out.txt] > table.tsv

Output columns (TSV, to stdout):
    sample | best_target | call | lineage | coinfection | coverage_pct | mean_depth | verdict

  best_target   reference the sample covered best (a subgroup, or 'none')
  call          final call from final_calls.txt (a subgroup, 'HMPVA/B', or 'no-call')
  lineage       A / B / A+B (from final_calls.txt)
  coinfection   yes / no
  coverage_pct  best good-cov% (fraction of genome >=10x), as a percentage
  mean_depth    mean depth at that best-covered reference
  verdict       PASS          -> got a confident call
                drop_lowcov   -> mapped to a reference but below coverage cutoff
                no_mash_call  -> too few reads for mash to assign a reference
"""
import sys
import csv
import os
import glob

MINCOV = 0.8  # matches call_hMPV.py -c; only used for labelling here

# per-target alignstats columns (no header):
#   0 sample 1 target 2 subfac 3 reads 4 aligned 5 paired
#   6 meandepth 7 goodcov 8 cov 9 gsize
COL_MEANDEPTH, COL_GOODCOV, COL_GSIZE = 6, 7, 9


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make_hMPV_table.py <samples.txt> <summary_dir> [nocall_out.txt]")
    samples_file = sys.argv[1]
    summary_dir = sys.argv[2]
    nocall_out = sys.argv[3] if len(sys.argv) > 3 else "nocall_samples.txt"

    # 1. all input samples (the full set)
    all_samples = [l.strip() for l in open(samples_file) if l.strip()]

    # 2. per-(sample, target) stats from every per-target alignstats file
    #    (skip the merged final_alignstats.txt, which has a different layout)
    stats = {}  # sample -> [(target, coverage_fraction, mean_depth), ...]
    for path in glob.glob(os.path.join(summary_dir, "*_alignstats.txt")):
        if os.path.basename(path) == "final_alignstats.txt":
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

    # 3. final subgroup calls
    #    final_calls.txt cols: sample | call | subgroup | lineage | coinfection
    calls = {}
    calls_path = os.path.join(summary_dir, "final_calls.txt")
    if os.path.exists(calls_path):
        with open(calls_path) as f:
            next(f, None)  # skip header
            for line in f:
                c = line.rstrip("\n").split("\t")
                if not c or not c[0]:
                    continue
                calls[c[0]] = {
                    "call": c[1] if len(c) > 1 else "",
                    "lineage": c[3] if len(c) > 3 else "",
                    "coinfection": c[4] if len(c) > 4 else "",
                }

    # 4. one tidy row per sample
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(
        ["sample", "best_target", "call", "lineage", "coinfection",
         "coverage_pct", "mean_depth", "verdict"]
    )
    nocall = []
    for s in all_samples:
        rows = stats.get(s, [])
        if rows:
            best_target, covpc, depth = max(rows, key=lambda r: r[1])
            covstr, depthstr = f"{covpc * 100:.1f}", f"{depth:.0f}"
        else:
            best_target, covstr, depthstr = "none", "0.0", "0"
        call = calls.get(s)
        if call:
            verdict = "PASS"
            writer.writerow([s, best_target, call["call"], call["lineage"],
                             call["coinfection"], covstr, depthstr, verdict])
        else:
            verdict = "drop_lowcov" if rows else "no_mash_call"
            nocall.append(s)
            writer.writerow([s, best_target, "no-call", "", "",
                             covstr, depthstr, verdict])

    with open(nocall_out, "w") as out:
        out.write("\n".join(nocall) + ("\n" if nocall else ""))

    sys.stderr.write(
        f"{len(all_samples)} samples | {len(all_samples) - len(nocall)} called | "
        f"{len(nocall)} no-call -> {nocall_out}\n"
    )


if __name__ == "__main__":
    main()
