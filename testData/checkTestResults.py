#!/usr/bin/env python3
"""Compare a rohan run on the simulated data against what was simulated.

rohan's MCMC is seeded from the wall clock (see randomProb() in libgab), so the
output is not bit reproducible and this checks the results are within tolerance
of the truth rather than comparing against a stored file.

Exits 0 if every check passes, 1 otherwise.
"""

import argparse
import gzip
import os
import sys

# The reported heterozygosity of a window is allowed to be off by this factor
# from the simulated rate.  Windows are 500kb and carry ~500 het sites, so the
# sampling noise alone is a few percent; the slack is for the model.
H_FACTOR = 2.0
# A window inside the ROH should have essentially no heterozygosity
H_ROH_MAX = 1e-4
# Percentage points of the genome that may be misclassified as in/out of ROH
ROH_PERCENT_TOLERANCE = 10.0

failures = []
checks = 0


def check(condition, message):
    global checks
    checks += 1
    if condition:
        print("  ok       %s" % message)
    else:
        print("  FAILED   %s" % message)
        failures.append(message)


def parse_truth(path):
    truth = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            truth[fields[0]] = fields[1:]
    return truth


def parse_summary(path):
    """Pull the two theta lines and the ROH percentage out of summary.txt."""
    summary = {}
    with open(path) as f:
        for line in f:
            if line.startswith("Genome-wide theta outside ROH:"):
                summary["thetaOutside"] = [float(x) for x in line.split("\t")[1:4]]
            elif line.startswith("Genome-wide theta inc. ROH:"):
                summary["thetaAll"] = [float(x) for x in line.split("\t")[1:4]]
            elif line.startswith("Segments in ROH(%)"):
                summary["rohPercent"] = float(line.split("\t")[1].split(" ")[0])
    return summary


def parse_hest(path):
    windows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            windows.append({
                "chr": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "h": float(fields[4]),
                "hLow": float(fields[6]),
                "hHigh": float(fields[7]),
            })
    return windows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, help="prefix of the rohan output")
    parser.add_argument("--truth", required=True, help="the .truth.txt written by simulateTestData.py")
    parser.add_argument("--roh-start", type=int, required=True, help="0-based start of the simulated ROH")
    parser.add_argument("--roh-end", type=int, required=True, help="0-based end of the simulated ROH")
    parser.add_argument("--theta", type=float, required=True, help="simulated het. rate outside the ROH")
    args = parser.parse_args()

    truth = parse_truth(args.truth)
    trueRohPercent = 100.0 * (args.roh_end - args.roh_start) / int(truth["chromosome"][1].split(" ")[0])
    # the realised rate, which differs a little from the requested --theta
    trueTheta = float(truth["het sites outside ROH"][1].split(" = ")[1])

    print("checking rohan output '%s' against '%s'" % (args.prefix, args.truth))

    print("\noutput files:")
    for suffix in [".hEst.gz", ".summary.txt", ".rginfo.gz",
                   ".mid.hmmp.gz", ".min.hmmp.gz", ".max.hmmp.gz",
                   ".het.pdf", ".hmm.pdf"]:
        path = args.prefix + suffix
        check(os.path.isfile(path) and os.path.getsize(path) > 0,
              "%s exists and is not empty" % path)

    if failures:
        # everything below reads these files, so there is no point going on
        print("\n%d of %d checks FAILED" % (len(failures), checks))
        return 1

    print("\nlocal heterozygosity (%s.hEst.gz):" % args.prefix)
    windows = parse_hest(args.prefix + ".hEst.gz")
    check(len(windows) > 0, "at least one window was estimated")

    for w in windows:
        label = "%s:%d-%d h=%g" % (w["chr"], w["start"], w["end"], w["h"])
        # a window counts as being in the ROH only if it lies entirely inside it
        if w["start"] >= args.roh_start and w["end"] <= args.roh_end:
            check(w["h"] < H_ROH_MAX,
                  "%s is in the ROH, h < %g" % (label, H_ROH_MAX))
        elif w["end"] <= args.roh_start or w["start"] >= args.roh_end:
            check(trueTheta / H_FACTOR < w["h"] < trueTheta * H_FACTOR,
                  "%s is outside the ROH, within %gx of %g" % (label, H_FACTOR, trueTheta))
        else:
            print("  skipped  %s straddles a ROH boundary" % label)

    print("\ngenome-wide results (%s.summary.txt):" % args.prefix)
    summary = parse_summary(args.prefix + ".summary.txt")
    check("thetaOutside" in summary, "summary.txt reports theta outside ROH")
    check("rohPercent" in summary, "summary.txt reports the fraction in ROH")

    if "rohPercent" in summary:
        check(abs(summary["rohPercent"] - trueRohPercent) <= ROH_PERCENT_TOLERANCE,
              "%.2f%% of the genome called ROH, simulated %.2f%% (tolerance %g points)" % (
                  summary["rohPercent"], trueRohPercent, ROH_PERCENT_TOLERANCE))

    if "thetaOutside" in summary:
        mid, low, high = summary["thetaOutside"]
        check(trueTheta / H_FACTOR < mid < trueTheta * H_FACTOR,
              "theta outside ROH %g, within %gx of the simulated %g" % (mid, H_FACTOR, trueTheta))
        check(low <= trueTheta <= high,
              "the simulated theta %g falls inside the reported [%g,%g]" % (trueTheta, low, high))

    print("")
    if failures:
        print("%d of %d checks FAILED:" % (len(failures), checks))
        for f in failures:
            print("  %s" % f)
        return 1

    print("all %d checks passed" % checks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
