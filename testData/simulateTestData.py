#!/usr/bin/env python3
"""Generate a small, self-contained data set for "make test".

Simulates a diploid chromosome that carries one long run of homozygosity in the
middle and a known heterozygosity rate everywhere else, then simulates single
end reads off the two haplotypes.  Everything is driven by a fixed seed so the
data (and therefore the expected result of the test) is reproducible.

Writes:
  <prefix>.fa   the reference (haplotype 1 without the simulated variants)
  <prefix>.sam  the aligned reads, coordinate sorted
  <prefix>.truth.txt  what was simulated, for the human reading the test output
"""

import argparse
import random
import sys

BASES = "ACGT"
# The transition partner of each base; used to keep the simulated Ts/Tv ratio
# in the same ballpark as rohan's default --tstv 2.1
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="simulated", help="output file prefix")
    parser.add_argument("--sam", default=None, help="where to write the SAM (default: <prefix>.sam)")
    parser.add_argument("--chr", default="testchr", help="chromosome name")
    parser.add_argument("--length", type=int, default=3000000, help="chromosome length in bp")
    parser.add_argument("--roh-start", type=int, default=1000000, help="0-based start of the ROH")
    parser.add_argument("--roh-end", type=int, default=2000000, help="0-based end of the ROH")
    parser.add_argument("--theta", type=float, default=1e-3, help="het. rate outside the ROH")
    parser.add_argument("--theta-roh", type=float, default=1e-6, help="het. rate inside the ROH")
    parser.add_argument("--coverage", type=float, default=8.0, help="average depth of coverage")
    parser.add_argument("--readlen", type=int, default=100, help="read length in bp")
    parser.add_argument("--qual", type=int, default=30, help="base quality of every base")
    parser.add_argument("--mapq", type=int, default=60, help="mapping quality of every read")
    parser.add_argument("--seed", type=int, default=42, help="random seed")
    parser.add_argument("--sample", default="testsample", help="sample name for the @RG line")
    return parser.parse_args()


def simulate_reference(rnd, length):
    """A random sequence at 40% GC, which is roughly what DNAprof assumes."""
    weights = (0.3, 0.2, 0.2, 0.3)
    return rnd.choices(BASES, weights=weights, k=length)


def simulate_variants(rnd, ref, args):
    """Return the second haplotype plus the list of heterozygous positions.

    Haplotype 1 is the reference itself, so every variant is a het site: the
    truth we are asking rohan to recover.
    """
    hap2 = list(ref)
    hetPositions = []

    for pos in range(len(ref)):
        inRoh = args.roh_start <= pos < args.roh_end
        theta = args.theta_roh if inRoh else args.theta
        if rnd.random() >= theta:
            continue
        refBase = ref[pos]
        # 2.1 transitions for every transversion, as per rohan's default --tstv
        if rnd.random() < 2.1 / 3.1:
            altBase = TRANSITION[refBase]
        else:
            altBase = rnd.choice([b for b in BASES if b != refBase and b != TRANSITION[refBase]])
        hap2[pos] = altBase
        hetPositions.append(pos)

    return hap2, hetPositions


def simulate_reads(rnd, ref, hap2, args, out):
    """Write coordinate sorted SAM records for reads drawn from both haplotypes."""
    errorRate = 10.0 ** (-args.qual / 10.0)
    qualString = chr(33 + args.qual) * args.readlen
    haplotypes = (ref, hap2)

    nReads = int(round(args.coverage * args.length / args.readlen))
    # Sorting positions up front is what makes the output coordinate sorted;
    # samtools then only has to compress it.
    starts = sorted(rnd.randrange(0, args.length - args.readlen) for _ in range(nReads))

    for i, start in enumerate(starts):
        hap = haplotypes[rnd.getrandbits(1)]
        seq = hap[start:start + args.readlen]

        if errorRate > 0:
            seq = list(seq)
            for j in range(args.readlen):
                if rnd.random() < errorRate:
                    seq[j] = rnd.choice([b for b in BASES if b != seq[j]])

        # Half the reads on the reverse strand.  The stored SEQ of a reverse
        # strand read is already the forward strand sequence, so only the flag
        # changes; that is enough to exercise rohan's strand handling.
        flag = 16 if rnd.getrandbits(1) else 0

        out.write("\t".join([
            "read%d" % i,
            str(flag),
            args.chr,
            str(start + 1),           # SAM is 1-based
            str(args.mapq),
            "%dM" % args.readlen,
            "*", "0", "0",
            "".join(seq),
            qualString,
            "RG:Z:testrg",
        ]) + "\n")

    return nReads


def write_fasta(path, name, seq, lineWidth=60):
    with open(path, "w") as f:
        f.write(">%s\n" % name)
        for i in range(0, len(seq), lineWidth):
            f.write("".join(seq[i:i + lineWidth]) + "\n")


def main():
    args = parse_args()
    rnd = random.Random(args.seed)

    ref = simulate_reference(rnd, args.length)
    hap2, hetPositions = simulate_variants(rnd, ref, args)

    write_fasta(args.prefix + ".fa", args.chr, ref)

    samPath = args.sam if args.sam else args.prefix + ".sam"
    with open(samPath, "w") as out:
        out.write("@HD\tVN:1.6\tSO:coordinate\n")
        out.write("@SQ\tSN:%s\tLN:%d\n" % (args.chr, args.length))
        out.write("@RG\tID:testrg\tSM:%s\tPL:ILLUMINA\tLB:testlib\n" % args.sample)
        nReads = simulate_reads(rnd, ref, hap2, args, out)

    outsideRoh = sum(1 for p in hetPositions if not (args.roh_start <= p < args.roh_end))
    insideRoh = len(hetPositions) - outsideRoh
    lengthOutsideRoh = args.length - (args.roh_end - args.roh_start)

    with open(args.prefix + ".truth.txt", "w") as f:
        f.write("chromosome\t%s\t%d bp\n" % (args.chr, args.length))
        f.write("ROH\t%s:%d-%d\t%d bp (%.1f%% of the chromosome)\n" % (
            args.chr, args.roh_start + 1, args.roh_end,
            args.roh_end - args.roh_start,
            100.0 * (args.roh_end - args.roh_start) / args.length))
        f.write("het sites outside ROH\t%d\ttheta = %g\n" % (
            outsideRoh, float(outsideRoh) / lengthOutsideRoh))
        f.write("het sites inside ROH\t%d\ttheta = %g\n" % (
            insideRoh, float(insideRoh) / (args.roh_end - args.roh_start)))
        f.write("het sites genome-wide\t%d\ttheta = %g\n" % (
            len(hetPositions), float(len(hetPositions)) / args.length))
        f.write("reads\t%d\tlength %d, target coverage %gx\n" % (
            nReads, args.readlen, args.coverage))
        f.write("seed\t%d\n" % args.seed)

    sys.stderr.write("simulated %d reads and %d het sites (%d in the ROH)\n" % (
        nReads, len(hetPositions), insideRoh))


if __name__ == "__main__":
    main()
