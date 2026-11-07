"""Frameshift signature: for each conserved query (Ec|/Bs|), pairs of HSPs (e<1e-5) on the same contig and strand,
in different frames, adjacent in both query (gap <= 30 aa, overlap allowed) and subject (gap <= 150 nt).
Reported per genome together with the number of query-locus hits (best HSP e<1e-10, qcov>=0.5).
Usage: python3 -I frameshift_check.py <dir> acc..."""
import sys, collections
R = sys.argv[1]
for acc in sys.argv[2:]:
    by = collections.defaultdict(list); loci = set()
    for line in open(f'{R}/out/{acc}.tblastn.tsv'):
        f = line.rstrip('\n').split('\t'); q = f[0]
        if not (q.startswith('Ec|') or q.startswith('Bs|')): continue
        s = f[1].split('|')[1] if '|' in f[1] else f[1]
        qs, qe, ss, se, ev, qlen, fr = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8]), int(f[10]), int(f[12])
        by[(q, s, fr > 0)].append((qs, qe, min(ss, se), max(ss, se), fr))
        if ev < 1e-10 and (qe - qs + 1) / qlen >= 0.5: loci.add((s, min(ss, se) // 3000))
    fs = set()
    for (q, s, st), L in by.items():
        L.sort()
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                a, b = L[i], L[j]
                if a[4] == b[4]: continue
                qgap = b[0] - a[1]; sgap = max(a[2], b[2]) - min(a[3], b[3])
                if -20 <= qgap <= 30 and sgap <= 150: fs.add((s, min(a[3], b[3]) // 300))
    print(f'{acc}: frameshift-like adjacent HSP pairs (distinct sites) = {len(fs)}; conserved loci ~{len(loci)}')
