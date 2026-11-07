"""Monte-Carlo dispersion test for in-alignment TAGs across contigs: place the observed number of TAGs
on aligned codons uniformly at random (multinomial by aligned codons per contig), 20000 times; compare
number of TAG-bearing contigs and max TAG per contig with observed. Usage: python3 -I contig_perm.py <dir> acc..."""
import sys, random
R = sys.argv[1]; random.seed(1)
for acc in sys.argv[2:]:
    rows = [l.rstrip('\n').split('\t') for l in open(f'{R}/out/{acc}.contig_aligned.tsv')][1:]
    w = [int(r[2]) for r in rows]; T = sum(int(r[4]) for r in rows)
    obs_n = sum(1 for r in rows if int(r[4]) > 0); obs_max = max(int(r[4]) for r in rows)
    cum = []; s = 0
    for x in w: s += x; cum.append(s)
    import bisect
    le_n = ge_max = 0; sims_n = []
    for _ in range(20000):
        c = {}
        for _ in range(T):
            k = bisect.bisect_right(cum, random.randrange(s)); c[k] = c.get(k, 0) + 1
        n = len(c); m = max(c.values())
        sims_n.append(n); le_n += n <= obs_n; ge_max += m >= obs_max
    sims_n.sort()
    print(f'{acc}: TAG={T} on {obs_n} contigs (null median {sims_n[10000]}, 2.5-97.5% {sims_n[500]}-{sims_n[19500]}; P(n<=obs)={le_n/20000:.3f}); '
          f'max per contig obs={obs_max}, P(max>=obs)={ge_max/20000:.3f}')
