"""Contamination/chimera checks. (1) copy number of single-copy marker loci: distinct genomic loci (merged HSP
spans, e<1e-10, qcov>=0.5) per gene name for ribosomal proteins + other single-copy genes; (2) for each test genome,
per-contig identity to the other test genome (blastn of 1020-nt fragments, out/ani/*.frag.fna) for TAG-bearing vs other contigs.
Usage: python3 -I contam_check.py <dir> acc..."""
import sys, collections, subprocess, re
R = sys.argv[1]
single = re.compile(r'^(rp[lsm][A-Z]|rpo[ABC]|gyrB|recA|infB|fusA|tsf|frr|pheS|pheT|alaS|valS|leuS|ileS|secY|ffh|ftsY|dnaG|smpB|era|obgE|obg|rbfA|tig|pyrG|nusA)$')
for acc in sys.argv[2:]:
    loci = collections.defaultdict(list)
    for line in open(f'{R}/out/{acc}.tblastn.tsv'):
        f = line.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
        gn = q.split('|')[1]
        if not single.match(gn): continue
        qs, qe, ss, se, ev, qlen = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10])
        if ev >= 1e-10 or (qe - qs + 1) / qlen < 0.5: continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        loci[gn].append((s, min(ss, se), max(ss, se)))
    multi = {}
    for gn, L in loci.items():
        L.sort(); m = []
        for s, a, b in L:
            if m and m[-1][0] == s and a <= m[-1][2]: m[-1][2] = max(m[-1][2], b)
            else: m.append([s, a, b])
        if len(m) > 1: multi[gn] = len(m)
    print(f'{acc}: single-copy markers found={len(loci)}, with >1 locus={len(multi)} {multi}')
