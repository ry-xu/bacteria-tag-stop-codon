"""Distances in the GTDB R232 bacterial tree among a family's genomes: pairwise between chosen genomes, the genera under
their MRCA, and each chosen genome's nearest family members. The newick is parsed iteratively (parent pointers +
branch lengths); quoted internal labels may contain ':'.
usage: python src/family_tree_dist.py data/bac120_r232.tree.gz <family list: GB_acc<TAB>g__genus> <GB_acc,GB_acc,...>"""
import gzip, sys, re
s = gzip.open(sys.argv[1], "rt").read().strip()
want = {l.split()[0]: l.split()[1] for l in open(sys.argv[2]) if l.strip()}
readers = sys.argv[3].split(",")
parent, blen, name = [-1], [0.0], [None]
cur, i, n = 0, 0, len(s)
tok = re.compile(r"[^(),:;]+")
while i < n:
    c = s[i]
    if c == "(":
        parent.append(cur); blen.append(0.0); name.append(None); cur = len(parent) - 1; i += 1
    elif c == ",":
        cur = parent[cur]; parent.append(cur); blen.append(0.0); name.append(None); cur = len(parent) - 1; i += 1
    elif c == ")":
        cur = parent[cur]; i += 1
    elif c == ":":
        m = re.match(r":([0-9.eE+-]+)", s[i:]); blen[cur] = float(m.group(1)); i += len(m.group(0))
    elif c == ";":
        break
    else:
        m = tok.match(s, i); lab = m.group(0); i = m.end()
        if lab.startswith("'"):  # quoted labels with support:taxon
            j = s.index("'", i) if not lab.endswith("'") or len(lab) == 1 else i
            if j != i: lab = s[m.start():j + 1]; i = j + 1
        name[cur] = lab
# first child of root bookkeeping: the root node is 0; the first '(' created node 1 under 0. fine.
leaf = {name[k]: k for k in range(len(name)) if name[k] in want}
def path(k):
    p, d = [], 0.0
    while k != -1:
        p.append((k, d)); d += blen[k]; k = parent[k]
    return p
paths = {a: dict(path(k)) for a, k in leaf.items()}
def dist(a, b):
    pa = paths[a]
    k, d = leaf[b], 0.0
    while k not in pa:
        d += blen[k]; k = parent[k]
    return d + pa[k], k
print(len(leaf), "of", len(want), "in tree")
for x in readers:
    if x not in leaf: print("missing", x)
rs = [x for x in readers if x in leaf]
for i, a in enumerate(rs):
    for b in rs[i + 1:]:
        d, k = dist(a, b)
        # how many family leaves under the MRCA, and how many are readers
        under = [x for x in leaf if k in paths[x]]
        print(f"{a} {want[a]} | {b} {want[b]} | dist {d:.3f} | mrca has {len(under)} family genomes: "
              + ",".join(sorted(set(want[x] for x in under))[:12]))
for r in rs:
    nn = sorted((dist(r, x)[0], x) for x in leaf if x != r)[:8]
    print(r, want[r], "nearest:", " ".join(f"{want[x][3:]}:{x[3:]}:{d:.3f}" for d, x in nn))
