"""Outbox figure, version 2 (after the full screen): (a) every bacterial genome tested, TAG read-through vs read-through at
TAA/TGA, sampled runs scaled up by their sampling rates; the TAG readers of lineages 1-6 labelled (lineage 6, whose RF1
has decayed rather than been lost, hollow), known UGA-recoding lineages marked; (b) GTDB R232 tree of family
Aristaeellaceae with the three RF1-losing genera marked.
usage: uv run python src/rank_readers.py && uv run python src/figure2.py  -> results/fig_readthrough_v2.png
Needs GTDB R232 metadata and tree (scripts/setup_reference_data.sh). Reads results/rt/*.stats.json if present, otherwise
the shipped data/screen/screened_genomes.tsv.gz (via rank_readers.load_stats)."""
import csv
import gzip
import json
import random
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rank_readers import load_stats

ROOT = Path(__file__).resolve().parent.parent
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#a3a29b", "#e6e5e0"
VIC, ARI, LAC = "#eb6834", "#2a78d6", "#1baf7a"  # categorical slots 2, 1, 3 of the reference palette
UGA = INK2  # known UGA lineages: a dark neutral, not a fourth hue
READERS = {  # acc: (lineage colour, label)
    "GCA_949290385.1": (VIC, "Victivallales · chicken, Ethiopia"),
    "GCA_965332555.1": (VIC, "Victivallales · red kangaroo, Australia"),
    "GCA_045599485.1": (ARI, "SFMI01 · Tibetan wild ass"),
    "GCA_024404275.1": (ARI, "SFMI01 · horse, France"),
    "GCA_030526225.1": (ARI, "JAUNDF01 · warthog, German zoo"),
    "GCA_045628485.1": (ARI, "ZC3RG03 · Tibetan horse"),
    "GCA_937907055.1": (ARI, "ZC3RG03 · horse, Inner Mongolia"),
    "GCA_056232285.1": (LAC, "Lachnospiraceae · Tibetan horse"),
    "GCA_056261325.1": (LAC, "Lachnospiraceae · Tibetan horse"),
    "GCA_056209055.1": (ARI, "HGM11604 · Tibetan wild ass (UAG = Trp?)"),  # lineage 6, drawn hollow: decayed RF1
    "QTP_TA20190011_bin.237": (ARI, "HGM11604 · a second wild ass"),  # authors' CNGB bins, not in NCBI
    "QTP_TA20190013_bin.296": (ARI, "HGM11604 · a third wild ass"),
}
HOLLOW = {"GCA_056209055.1", "QTP_TA20190011_bin.237", "QTP_TA20190013_bin.296"}


def taxonomy():
    t = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(ROOT / "data" / f"{dom}_metadata.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                t[r["accession"][3:]] = r["gtdb_taxonomy"]
    return t


def archaeal_mags():
    """MAGs newer than GTDB (equid/Tibetan projects) whose NCBI name is archaeal, e.g. a Code-34 Methanimicrococcus."""
    kw = ("methan", "archaeon", "thermoplasm", "halobact", "nitrososphaer", "euryarch", "thermoproteot")
    out = set()
    for p in ("PRJNA1030165_genomes.tsv", "PRJNA438436_genomes.tsv"):
        for line in open(ROOT / "data" / "screen" / "lists" / p):
            x = line.rstrip("\n").split("\t")
            if len(x) > 1 and any(k in x[1].lower() for k in kw):
                out.add(x[0])
    return out


def points(tax):
    uga = {l.split("\t")[1] for l in list(open(ROOT / "results" / "readers_TGA.tsv"))[1:]
           if float(l.split("\t")[4]) >= 0.1 and float(l.split("\t")[6]) >= 8 and int(l.split("\t")[2]) >= 10
           and float(l.split("\t")[7]) >= 4}
    arch = archaeal_mags()
    pts, seen = [], set()
    for d in sorted(load_stats(), key=lambda d: d["acc"]):  # GCA_ sorts before GCF_
        if tax.get(d["acc"], "").startswith("d__Archaea") or d["acc"] in arch:
            continue
        num = d["acc"].split("_", 1)[-1].split(".")[0]  # GenBank and RefSeq copies of one assembly count once
        if d["acc"][:4] in ("GCA_", "GCF_") and num in seen:
            continue
        seen.add(num)
        rates = d.get("rates") or {}
        est = {s: d["readthrough"][s] / rates.get(s, 1.0) for s in ("TAA", "TAG", "TGA")}
        kind = "reader" if d["acc"] in READERS else "uga" if d["acc"] in uga else "other"
        pts.append((d["acc"], max(est["TAA"], est["TGA"]), est["TAG"], kind))
    return pts


def parse_tree(fn):
    s = gzip.open(fn, "rt").read().strip()
    parent, blen, name, cur, i = [-1], [0.0], [None], 0, 0
    tok, num = re.compile(r"[^(),:;]+"), re.compile(r":([0-9.eE+-]+)")
    while i < len(s):
        c = s[i]
        if c == "(":
            parent.append(cur); blen.append(0.0); name.append(None); cur = len(parent) - 1; i += 1
        elif c == ",":
            cur = parent[cur]; parent.append(cur); blen.append(0.0); name.append(None); cur = len(parent) - 1; i += 1
        elif c == ")":
            cur = parent[cur]; i += 1
        elif c == ":":
            m = num.match(s, i); blen[cur] = float(m.group(1)); i = m.end()
        elif c == ";":
            break
        elif c == "'":  # quoted internal label (support:taxon) may contain ':'
            j = s.index("'", i + 1); name[cur] = s[i + 1:j]; i = j + 1
        else:
            m = tok.match(s, i); name[cur] = m.group(0); i = m.end()
    return parent, blen, name


def family_subtree():
    fam = {l.split()[0]: l.split()[1] for l in open(ROOT / "data" / "screen" / "lists" / "aristaeellaceae.txt") if l.strip()}
    parent, blen, name = parse_tree(ROOT / "data" / "bac120_r232.tree.gz")
    leaves = [k for k, n in enumerate(name) if n in fam]
    anc = None  # MRCA = deepest node on every leaf's root path
    for k in leaves:
        p, x = [], k
        while x != -1:
            p.append(x); x = parent[x]
        anc = p[::-1] if anc is None else [a for a, b in zip(anc, p[::-1]) if a == b]
    root = anc[-1]
    children = {}
    keep = set()
    for k in leaves:
        x = k
        while x != root and x not in keep:
            keep.add(x); x = parent[x]
    for k in keep:
        children.setdefault(parent[k], []).append(k)
    return root, children, blen, name, fam


def draw_tree(ax, root, children, blen, name, fam):
    order, depth = [], {root: 0.0}
    stack = [root]
    while stack:  # depth-first, children in input order
        k = stack.pop()
        for c in children.get(k, []):
            depth[c] = depth[k] + blen[c]
        if k not in children:
            order.append(k)
        stack.extend(reversed(children.get(k, [])))
    xpos = {k: i for i, k in enumerate(order)}

    def x_of(k):
        if k in xpos:
            return xpos[k]
        xs = [x_of(c) for c in children[k]]
        xpos[k] = (min(xs) + max(xs)) / 2
        return xpos[k]
    x_of(root)
    for k, cs in children.items():
        xs = [xpos[c] for c in cs]
        ax.plot([min(xs), max(xs)], [depth[k]] * 2, color=MUTED, lw=0.5, solid_capstyle="butt")
        for c in cs:
            ax.plot([xpos[c]] * 2, [depth[k], depth[c]], color=MUTED, lw=0.5, solid_capstyle="butt")
    labels = {"GB_GCA_024404275.1": ("SFMI01: horse, wild ass", "left"), "GB_GCA_045599485.1": None,
              "GB_GCA_937907055.1": ("ZC3RG03: horses", "right"), "GB_GCA_030526225.1": ("JAUNDF01: warthog", "center")}
    base = max(depth[k] for k in order) + 0.03  # labels on one baseline under the deepest tip, with leader lines
    for k in order:
        if name[k] in labels:
            ax.scatter([xpos[k]], [depth[k] + 0.012], s=46, color=ARI, edgecolor=SURFACE, linewidth=1.5, zorder=5)
            if labels[name[k]]:
                lab, ha = labels[name[k]]
                ax.plot([xpos[k]] * 2, [depth[k] + 0.02, base], color=ARI, lw=0.8)
                ax.annotate(lab, (xpos[k], base), ha=ha, va="top", fontsize=12, color=INK,
                            xytext=({"left": -4, "right": 4, "center": 0}[ha], -3), textcoords="offset points")
    ax.invert_yaxis()
    ax.set_xlim(-5, len(order) + 5)
    ax.axis("off")
    return len(order)


def main():
    tax = taxonomy()
    pts = points(tax)
    random.seed(1)
    fig = plt.figure(figsize=(12, 12.2), facecolor=SURFACE)
    ax = fig.add_axes([0.08, 0.40, 0.88, 0.37], facecolor=SURFACE)
    # sampled runs scaled up by their sampling rate land on a coarse grid; spread them a little (log-space jitter)
    jit = lambda v: v * 10 ** random.uniform(-0.07, 0.07) if v > 0 else random.uniform(-0.06, 0.06)
    for kind, col, size, z in (("other", MUTED, 9, 1), ("uga", UGA, 20, 2)):
        xs = [jit(p[1]) for p in pts if p[3] == kind]
        ys = [jit(p[2]) for p in pts if p[3] == kind]
        ax.scatter(xs, ys, s=size, color=col, alpha=0.55 if kind == "other" else 0.9, linewidth=0, zorder=z)
    # labels in one column right of the points, evenly spaced in log y, with thin leader lines
    rd = sorted(((y, x, acc) for acc, x, y, kind in pts if kind == "reader"), reverse=True)
    slots = [10 ** (2.72 - i * (2.72 - 1.38) / (len(rd) - 1)) for i in range(len(rd))]
    for (y, x, acc), ys in zip(rd, slots):
        col, lab = READERS[acc]
        if acc in HOLLOW:
            ax.scatter([x], [y], s=52, facecolor=SURFACE, edgecolor=col, linewidth=1.8, zorder=5)
        else:
            ax.scatter([x], [y], s=60, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=5)
        ax.plot([x * 1.12, 11, 13], [y, ys, ys], color=col, lw=0.8, zorder=4)
        ax.annotate(lab, (14, ys), fontsize=11.5, color=INK, va="center", ha="left", zorder=6,
                    bbox=dict(boxstyle="square,pad=0.1", fc=SURFACE, ec="none", alpha=0.85))
    ax.set_xscale("symlog", linthresh=1, linscale=0.4)
    ax.set_yscale("symlog", linthresh=1, linscale=0.4)
    ax.set_xlim(-0.2, 3000)
    ax.set_ylim(-0.2, 700)
    ticks = [0, 1, 10, 100, 1000]
    ax.set_xticks(ticks, [f"{t:,}" for t in ticks])
    ax.set_yticks(ticks[:4], [f"{t:,}" for t in ticks[:4]])
    ax.minorticks_off()
    ax.tick_params(colors=INK2, labelsize=12, length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlabel("Genes read through at TAA or TGA (whichever is higher)", color=INK2, fontsize=13)
    ax.set_ylabel("Genes read through at TAG", color=INK2, fontsize=13)
    n_other = sum(1 for p in pts if p[3] == "other")
    n_uga = sum(1 for p in pts if p[3] == "uga")
    fig.text(0.04, 0.965, "Bacteria have stopped stopping at TAG at least six separate times",
             fontsize=21, weight="bold", color=INK)
    fig.text(0.04, 0.935, "Each dot is a genome; a gene is read through when its protein continues past the stop codon.",
             fontsize=13, color=INK2)
    fig.text(0.04, 0.913, "The outliers have lost release factor 1, which ends translation at TAG, or (hollow) carry one with a broken stop motif.",
             fontsize=13, color=INK2)
    legend = [("●  Victivallales", VIC, " (Verrucomicrobiota): two genomes, one lineage; what TAG encodes is unresolved"),
              ("●  Aristaeellaceae", ARI, " (Bacillota): five genomes, three genera; TAG appears to code for glutamate"),
              ("●  Lachnospiraceae", LAC, " (Bacillota): two genomes, one lineage; TAG appears to code for glutamine"),
              ("○  Aristaeellaceae HGM11604", ARI, ": three genomes, one species; TAG probably codes for tryptophan"),
              (f"●  {n_uga} genomes", UGA, " of lineages already known to read TGA as tryptophan or glycine"),
              (f"●  {n_other:,}", MUTED, " other bacterial genomes (most tested on a random subsample of genes, scaled up)")]
    for i, (head, col, rest) in enumerate(legend):
        y = 0.885 - 0.021 * i
        t = fig.text(0.04, y, head, fontsize=12.5, color=col, weight="bold")
        fig.canvas.draw()
        bb = t.get_window_extent().transformed(fig.transFigure.inverted())
        fig.text(bb.x1 + 0.004, y, rest, fontsize=12.5, color=INK2)

    tx = fig.add_axes([0.06, 0.03, 0.90, 0.245], facecolor=SURFACE)
    n = draw_tree(tx, *family_subtree())
    fig.text(0.04, 0.325, f"Family Aristaeellaceae ({n} species in the GTDB tree): three genera, three separate losses",
             fontsize=15, weight="bold", color=INK)
    fig.text(0.04, 0.305, "All 1,705 genomes of the family were tested: the 5 TAG readers each sit among close relatives that keep release factor 1.",
             fontsize=13, color=INK2)
    fig.text(0.04, 0.285, "The hollow lineage above is a fourth genus of the family, HGM11604; its three genomes are newer than this tree.",
             fontsize=13, color=INK2)
    out = ROOT / "results" / "fig_readthrough_v2.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(out, n_other, n_uga)


if __name__ == "__main__":
    main()
