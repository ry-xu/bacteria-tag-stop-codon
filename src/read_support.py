"""Raw-read support for in-frame TAG sites of a MAG: stream reads of an SRA run (NCBI AWS mirror via fastq-dump,
nothing stored but hits), map to the MAG with minimap2, and count the codon each read shows at every in-frame TAG
(genes called with TAG read as sense, table 15). A real reassignment shows TAG fixed in reads; assembly errors don't.

usage: python src/read_support.py <SRR> <mag_acc> <total_spots> [streams=4]  -> results/read_support/<mag>_<SRR>.txt
       <SRR> may also be a GSA run path such as gsa2/CRA013436/CRR938779 (CNCB download server): then the first
       <total_spots> read pairs of its _f1/_r2 files are streamed (nothing stored but hits).
Reads with > MAXDIV mismatches per aligned base (NM tag, default 3%) are ignored, so relatives' reads don't count.
<mag_acc> may be a comma-separated list: reads are then mapped competitively to all of them at once and each MAG
gets its own report (e.g. a UAG=Trp? candidate and a UAG=Glu reader from the same sample).
The same tally is made for in-frame TGG (Trp) codons as a control: the codons seen at sites where the MAG has TAG
vs those seen at its TGG sites show which sense codon TAG swaps with between strains.
"""
import bisect
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import ROOT, fetch, read_fasta, revcomp

OUT = ROOT / "results" / "read_support"
MAXDIV = 0.03


def tag_sites(contigs, codon="TAG"):
    finder = pyrodigal.GeneFinder()
    finder.train(*[s.encode() for s in contigs.values()], translation_table=15)
    sites = {}
    for c, s in contigs.items():
        for g in finder.find_genes(s.encode()):
            nt = g.sequence()
            for i in range(0, len(nt) - 3, 3):
                if nt[i:i + 3] == codon:
                    p0 = g.begin - 1 + i if g.strand == 1 else g.end - 3 - i
                    sites[(c, p0)] = g.strand
    return sites


def tally(work, sites):
    """Codon seen in each (high-identity) read at each site."""
    starts = defaultdict(list)
    for (c, p0) in sites:
        starts[c].append(p0)
    for c in starts:
        starts[c].sort()
    obs = defaultdict(Counter)
    mapped = 0
    for p in sorted(work.glob("part_*.sam")):
        for line in open(p):
            if line[0] == "@":
                continue
            f = line.split("\t")
            if len(f) < 11 or int(f[1]) & 0x904 or f[2] not in starts:
                continue
            nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
            c, pos, cigar, seq = f[2], int(f[3]) - 1, f[5], f[9]
            refmap, r, q = {}, pos, 0
            for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                ln = int(ln)
                if op in "M=X":
                    for k in range(ln):
                        refmap[r + k] = q + k
                    r += ln
                    q += ln
                elif op in "IS":
                    q += ln
                elif op in "DN":
                    r += ln
            if nm > MAXDIV * max(1, len(refmap)):
                continue
            mapped += 1
            lo, hi = bisect.bisect_left(starts[c], pos), bisect.bisect_right(starts[c], r - 3)
            for p0 in starts[c][lo:hi]:
                if all(p0 + k in refmap for k in range(3)):
                    cod = "".join(seq[refmap[p0 + k]] for k in range(3))
                    obs[(c, p0)][revcomp(cod) if sites[(c, p0)] == -1 else cod] += 1
    return obs, mapped


def summarise(sites, obs, codon):
    covered = [k for k in sites if sum(obs[k].values()) >= 3]
    frac = {k: obs[k][codon] / sum(obs[k].values()) for k in covered}
    depth = sorted(sum(obs[k].values()) for k in covered)
    alt = Counter()
    for k in covered:
        if frac[k] < 0.9:
            alt.update({c: n for c, n in obs[k].items() if c != codon})
    return [f"in-frame {codon} sites: {len(sites)}; covered by >= 3 reads: {len(covered)}"
            f" (median depth {depth[len(depth) // 2] if depth else 0})",
            f"{codon} in >= 90% of reads: {sum(v >= 0.9 for v in frac.values())};"
            f" mixed: {sum(0.1 < v < 0.9 for v in frac.values())}; {codon} in <= 10%: {sum(v <= 0.1 for v in frac.values())}",
            f"other codons at non-fixed {codon} sites: {alt.most_common(10)}"]


def main(run, acc, total, streams=4):
    OUT.mkdir(parents=True, exist_ok=True)
    accs = acc.split(",")
    genomes = {a: dict(read_fasta(fetch(a))) for a in accs}
    work = OUT / f"{'+'.join(accs)}_{run.rstrip('/').split('/')[-1]}"
    work.mkdir(exist_ok=True)
    ref = work / "mag.fa"
    ref.write_text("".join(f">{k}\n{v}\n" for g in genomes.values() for k, v in g.items()))
    procs = []
    if run.startswith("gsa"):
        name = run.rstrip("/").split("/")[-1]
        base = f"https://download.cncb.ac.cn/{run.rstrip('/')}/{name}"
        for i, end in enumerate(("_f1.fq.gz", "_r2.fq.gz")):
            cmd = (f"{sys.executable} {ROOT / 'src' / 'range_stream.py'} {base}{end} | gzip -dc 2>/dev/null | head -n {4 * total} | "
                   f"minimap2 -x sr -a -t 2 --sam-hit-only {ref} - 2>/dev/null > {work}/part_{i}.sam")
            procs.append(subprocess.Popen(cmd, shell=True))
        run = name
    else:
        per = total // streams
        for i in range(streams):
            a, b = i * per + 1, (i + 1) * per
            cmd = (f"fastq-dump --stdout --split-spot --skip-technical -N {a} -X {b} {run} 2>/dev/null | "
                   f"minimap2 -x sr -a -t 1 --sam-hit-only {ref} - 2>/dev/null > {work}/part_{i}.sam")
            procs.append(subprocess.Popen(cmd, shell=True))
    for p in procs:
        p.wait()
    for a, contigs in genomes.items():
        sites, trp = tag_sites(contigs), tag_sites(contigs, "TGG")
        obs, mapped = tally(work, sites)
        tobs, _ = tally(work, trp)
        others = f"; mapped competitively with {', '.join(x for x in accs if x != a)}" if len(accs) > 1 else ""
        lines = [f"{a} reads from {run} (first {total} spots): {mapped} mapped reads (<= {MAXDIV:.0%} mismatches){others}"]
        lines += summarise(sites, obs, "TAG") + ["control:"] + summarise(trp, tobs, "TGG")
        per_site = [f"{c}\t{p0}\t{sites[(c, p0)]}\t" + " ".join(f"{k}:{n}" for k, n in obs[(c, p0)].most_common())
                    for (c, p0) in sorted(sites) if obs[(c, p0)]]
        tag = a if len(accs) == 1 else f"{a}_comp"
        (OUT / f"{tag}_{run}.txt").write_text("\n".join(lines) + "\n")
        (OUT / f"{tag}_{run}.sites.tsv").write_text("contig\tpos0\tstrand\tcodons\n" + "\n".join(per_site) + "\n")
        print("\n".join(lines))
    for p in work.glob("part_*.sam"):
        p.unlink()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 4)
