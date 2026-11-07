"""Pool read_support.py per-site tables of one MAG across several runs (e.g. different animals carrying the species).
Per TAG site: in how many runs it is covered (>= 3 reads) and fixed (TAG in >= 90% of reads); and, counting sites
rather than reads, which codon replaces TAG at polymorphic sites (alternative codon seen in >= 2 reads of one run).
usage: python src/read_support_pool.py <mag_tag> <run> [<run> ...]   (reads results/read_support/<mag_tag>_<run>.sites.tsv)"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main(tag, runs):
    obs = defaultdict(dict)
    for run in runs:
        for line in list(open(ROOT / "results" / "read_support" / f"{tag}_{run}.sites.tsv"))[1:]:
            c, p, st, cods = line.rstrip("\n").split("\t")
            obs[(c, int(p))][run] = Counter({k: int(v) for k, v in (x.rsplit(":", 1) for x in cods.split())})
    ncov, nfix_all, nfix_any, nabsent = Counter(), 0, 0, 0
    alt_sites = Counter()
    for k, per in obs.items():
        cov = {r: o for r, o in per.items() if sum(o.values()) >= 3}
        if not cov:
            continue
        ncov[len(cov)] += 1
        fr = [o["TAG"] / sum(o.values()) for o in cov.values()]
        nfix_all += all(f >= 0.9 for f in fr)
        nfix_any += any(f >= 0.9 for f in fr)
        nabsent += any(f <= 0.1 for f in fr)
        alts = {c for o in per.values() for c, n in o.items() if c != "TAG" and n >= 2}
        alt_sites.update(alts)
    n = sum(ncov.values())
    print(f"{tag}: {len(obs)} TAG sites with reads in {len(runs)} runs; {n} covered (>= 3 reads) in >= 1 run "
          f"(by number of runs: {dict(sorted(ncov.items()))})")
    print(f"  TAG fixed (>= 90% of reads) in every covering run: {nfix_all}; in at least one: {nfix_any}; "
          f"TAG <= 10% in any run: {nabsent}")
    print(f"  sites with an alternative codon in >= 2 reads of one run: {sum(alt_sites.values())} -> {alt_sites.most_common(8)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
