#!/usr/bin/env python
"""Assign hMPV subgroup calls from per-target alignment stats.

Unlike the RSV caller (a fixed A-vs-B pivot on two references), hMPV runs
against a pool of N subgroup references (e.g. A2a, A2b2WT, A2b2-111dup, B1,
B2b). A sample is aligned against every reference mash assigns it to, so it can
pass coverage on several closely-related *within-lineage* references at once.
So a 2-way ratio doesn't apply. Instead we:

  1. keep every (sample, target) alignment passing --coverage (good-cov fraction),
  2. call the SUBGROUP = the passing target with the highest good-coverage
     (mean depth breaks ties) — this is the best-matching reference,
  3. flag an A/B COINFECTION only when the sample passes on BOTH an A-lineage
     and a B-lineage reference and their good-cov ratio is within --ratio
     (neither lineage clearly dominates).

Lineage is read from the target name: HMPVA* -> A, HMPVB* -> B.

Outputs (mirrors the RSV caller's two files, consumed by rule call_hmpv):
  {out}_calls.txt       header:  sample  call  subgroup  lineage  coinfection
  {out}_alignstats.txt  headerless: sample target call <alignstats cols> covpc goodcovpc lineage
                        one row per passing (sample, target)
"""
import sys
import argparse
import pandas as pd

# per-target alignstats columns (no header), as emitted by rules/stats.smk:
ALIGNSTATS_COLS = ["sample", "target", "subfactor", "reads", "aligned",
                   "paired", "meandepth", "goodcov", "cov", "gsize"]


def lineage_of(target):
    t = str(target).upper()
    if t.startswith("HMPVA"):
        return "A"
    if t.startswith("HMPVB"):
        return "B"
    return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignstats", "-i", required=True, help="catted per-target alignstats")
    ap.add_argument("--coverage", "-c", type=float, default=0.8,
                    help="min good-coverage fraction to keep an alignment")
    ap.add_argument("--ratio", "-r", type=float, default=0.95,
                    help="min A/B good-cov ratio to call an A+B coinfection")
    ap.add_argument("--out", "-o", required=True, help="output prefix")
    args = ap.parse_args()

    df = pd.read_csv(args.alignstats, sep="\t", names=ALIGNSTATS_COLS)
    df["goodcovpc"] = df["goodcov"] / df["gsize"]
    df["covpc"] = df["cov"] / df["gsize"]
    df["lineage"] = df["target"].map(lineage_of)

    n_before = len(df)
    passing = df[df["goodcovpc"] >= args.coverage].copy()
    print("keeping {}/{} alignments with good-cov >= {}".format(
        len(passing), n_before, args.coverage), file=sys.stderr)

    calls = []
    for sample, g in passing.groupby("sample"):
        g = g.sort_values(["goodcovpc", "meandepth"], ascending=False)
        best = g.iloc[0]
        subgroup = best["target"]
        lineage = best["lineage"]

        # best good-cov per lineage (NaN if the sample never passes that lineage)
        bestA = g.loc[g["lineage"] == "A", "goodcovpc"].max()
        bestB = g.loc[g["lineage"] == "B", "goodcovpc"].max()

        coinf = False
        if pd.notna(bestA) and pd.notna(bestB):
            ratio = bestA / bestB
            if args.ratio < ratio < (1.0 / args.ratio):
                coinf = True
                lineage = "A+B"

        call = "HMPVA/B" if coinf else subgroup
        calls.append({"sample": sample, "call": call, "subgroup": subgroup,
                      "lineage": lineage, "coinfection": "yes" if coinf else "no"})

    calls_df = pd.DataFrame(
        calls, columns=["sample", "call", "subgroup", "lineage", "coinfection"])
    calls_df.to_csv("{}_calls.txt".format(args.out), sep="\t", index=False)

    # per-passing-alignment stats, with the sample's final call attached
    call_map = {c["sample"]: c["call"] for c in calls}
    out = passing[passing["sample"].isin(call_map)].copy()
    out["call"] = out["sample"].map(call_map)
    out_cols = (["sample", "target", "call"]
                + [c for c in ALIGNSTATS_COLS if c not in ("sample", "target")]
                + ["covpc", "goodcovpc", "lineage"])
    out[out_cols].to_csv("{}_alignstats.txt".format(args.out),
                         sep="\t", index=False, header=False)

    n_coinf = sum(1 for c in calls if c["coinfection"] == "yes")
    print("{} samples called ({} A/B coinfections)".format(len(calls), n_coinf),
          file=sys.stderr)


if __name__ == "__main__":
    main()
