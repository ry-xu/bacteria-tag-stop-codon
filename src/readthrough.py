"""Read-through test: does protein homology continue past a gene's stop codon?

For one genome: call genes with pyrodigal (table 11), and for every gene read on past its stop codon, in frame, to
the next stop (max 300 codons). The joined protein (upstream + 'X' + downstream) is aligned to Swiss-Prot with
DIAMOND. A gene "reads through" when one alignment covers the stop with >= FLANK residues on each side. For each
such gene we keep the stop codon, the Swiss-Prot subject and the subject position the stop aligns to, so that
recoded sites can later be matched across genomes (the same site in the same family = conserved in-frame stop, the
Sec/Pyl signature; scattered sites = pseudogenes or sequencing errors).

usage: uv run python src/readthrough.py <acc> [<acc> ...]
       -> results/rt/<acc>.tsv (one row per read-through gene) and results/rt/<acc>.stats.json
"""
import gzip
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import zipfile
import zlib
from pathlib import Path

import pyrodigal
import requests

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "rt"
# Genome cache and Swiss-Prot DIAMOND database; override with env vars to reuse an existing cache.
GEN = Path(os.environ.get("GENOME_DIR", ROOT / "data" / "genomes"))
DB = Path(os.environ.get("SPROT_DB", ROOT / "data" / "sprot.dmnd"))
FLANK = 15
MAXEXT = 300
MINEXT = 15
STOPS = {"TAA", "TAG", "TGA"}
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODE[a + b + c] = _aa[16 * i + 4 * j + k]
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def fetch_ftp(acc, fn):
    """Fallback when the Datasets API is down: the assembly's *_genomic.fna.gz from the NCBI FTP site (over https)."""
    d = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/"
    for attempt in range(4):  # the server answers 503 "Service unavailable" pages when busy: back off and retry
        if attempt:
            time.sleep(30 * attempt)
        try:
            r = requests.get(d, timeout=60)
            sub = sorted(set(re.findall(rf'href="({re.escape(acc)}_[^"/]+)/?"', r.text)))
            if r.status_code != 200 or not sub:
                continue
            r = requests.get(f"{d}{sub[0]}/{sub[0]}_genomic.fna.gz", timeout=300)
            if r.status_code == 200 and r.content[:2] == b"\x1f\x8b":
                fn.write_bytes(r.content)
                return fn
        except requests.RequestException:
            pass
    raise RuntimeError(f"download failed: {acc}")


_api_down_until = 0


def fetch(acc):
    """Genome FASTA from the NCBI Datasets API (FTP site if the API fails), cached in data/genomes/ (gzipped)."""
    GEN.mkdir(parents=True, exist_ok=True)
    fn = GEN / f"{acc}.fna.gz"
    if fn.exists():
        return fn
    global _api_down_until
    if time.time() < _api_down_until:
        return fetch_ftp(acc, fn)
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
    for attempt in range(2):
        try:
            r = requests.get(url, timeout=300)
            if r.status_code == 200:
                break
        except requests.RequestException:
            pass
        time.sleep(5 * (attempt + 1))
    else:
        _api_down_until = time.time() + 600  # don't retry the API for every genome during an outage
        return fetch_ftp(acc, fn)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    names = [n for n in z.namelist() if n.endswith(".fna")]
    if not names:
        raise RuntimeError(f"no fna for {acc}")
    fn.write_bytes(gzip.compress(z.read(names[0])))
    return fn


def read_fasta(fn):
    name, buf = None, []
    with gzip.open(fn, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def joined_proteins(contigs):
    """Yield (gene_id, stop_codon, upstream_aa_len, joined_protein, ext_hit_contig_end, next_gene_codons).
    next_gene_codons: codons from the stop to the start of the next called gene on the same strand and in the same
    frame inside the extension (-1 if none). A short spacer + a new gene is the signature of two adjacent genes whose
    Swiss-Prot homolog is a natural fusion, not of read-through."""
    seqs = dict(contigs)
    total = sum(len(s) for s in seqs.values())
    finder = pyrodigal.GeneFinder(meta=total < 100000)
    if total >= 100000:
        finder.train(*[s.encode() for s in seqs.values()], translation_table=11)
    for cname, seq in seqs.items():
        genes = finder.find_genes(seq.encode())
        rc = None
        L = len(seq)
        starts = {1: sorted(g.begin - 1 for g in genes if g.strand == 1),
                  -1: sorted(L - g.end for g in genes if g.strand == -1)}  # 0-based starts in each strand's own coords
        for i, g in enumerate(genes):
            if g.partial_end if g.strand == 1 else g.partial_begin:
                continue  # no stop codon
            if g.strand == 1:
                s = seq
                start, end = g.begin - 1, g.end  # 0-based, end exclusive, includes stop
            else:
                if rc is None:
                    rc = revcomp(seq)
                s = rc
                L = len(seq)
                start, end = L - g.end, L - g.begin + 1
            cds = s[start:end]
            stop = cds[-3:]
            if stop not in STOPS or len(cds) % 3:
                continue
            ext, j, hit_end = [], end, True
            while j + 3 <= len(s) and len(ext) < MAXEXT:
                c = s[j:j + 3]
                if c in STOPS:
                    hit_end = False
                    break
                ext.append(CODE.get(c, "X"))
                j += 3
            if len(ext) < MINEXT:
                continue
            up = translate(cds[:-3])
            nxt = [(x - end) // 3 for x in starts[g.strand] if end <= x < end + 3 * len(ext) and (x - end) % 3 == 0]
            yield (f"{cname}|{g.begin}|{g.end}|{g.strand}", stop, len(up), up + "X" + "".join(ext), hit_end,
                   nxt[0] if nxt else -1)


def keep(gid, stop, rates):
    """Deterministic subsampling of genes to test, by stop codon (rates: {stop: fraction}); None = test all."""
    if not rates:
        return True
    r = rates.get(stop, 1.0)
    return r >= 1.0 or (zlib.crc32(gid.encode()) % 1000) < r * 1000


OUTFMT = ["qseqid", "sseqid", "pident", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qseq",
          "sseq", "stitle"]


def prepare(acc, rates=None):
    """Download, call genes, build the joined proteins to test. -> (queries [(gid, prot)], meta, stats skeleton)."""
    fn = fetch(acc)
    contigs = list(read_fasta(fn))
    meta, queries = {}, []
    n_stop = {"TAA": 0, "TAG": 0, "TGA": 0}
    for gid, stop, upl, prot, hit_end, nxt in joined_proteins(contigs):
        if not keep(gid, stop, rates):
            continue
        meta[gid] = (stop, upl, len(prot), hit_end, nxt)
        n_stop[stop] += 1
        queries.append((gid, prot))
    stats = dict(acc=acc, contigs=len(contigs), bp=sum(len(s) for _, s in contigs), tested=n_stop,
                 readthrough=None, rates=rates or {})
    return queries, meta, stats


def finish(acc, meta, stats, hit_lines):
    """Keep alignments that span the stop with >= FLANK residues each side; write <acc>.tsv and <acc>.stats.json."""
    rows = {}
    for line in hit_lines:
        qid, sid, pid, qs, qe, ss, se, ev, bs, ql, sl, qseq, sseq, stitle = line.rstrip("\n").split("\t")
        stop, upl, L, hit_end, nxt = meta[qid]
        qs, qe, ss = int(qs), int(qe), int(ss)
        xpos = upl + 1  # 1-based position of the X
        if not (qs <= xpos - FLANK and qe >= xpos + FLANK):
            continue
        # subject position aligned to the X, and the subject residue there
        qi, si, spos, sres, k = qs - 1, ss - 1, None, None, None
        for k, (a, b) in enumerate(zip(qseq, sseq)):
            if a != "-":
                qi += 1
            if b != "-":
                si += 1
            if qi == xpos and a != "-":
                spos, sres = (si, b) if b != "-" else (None, "-")
                break
        qwin, swin = qseq[max(0, k - 12):k + 13], sseq[max(0, k - 12):k + 13]
        if qid not in rows or float(bs) > float(rows[qid]["bitscore"]):
            rows[qid] = dict(acc=acc, gene=qid, stop=stop, up_len=upl, len=L, ext_to_contig_end=int(hit_end),
                             sseqid=sid, stitle=stitle.split(" OS=")[0], pident=pid, qstart=qs, qend=qe,
                             evalue=ev, bitscore=bs, spos=spos, sres=sres, next_gene_codons=nxt,
                             gaps_win=qwin.count("-") + swin.count("-"), qwin=qwin, swin=swin)
    cols = ["acc", "gene", "stop", "up_len", "len", "ext_to_contig_end", "sseqid", "stitle", "pident", "qstart",
            "qend", "evalue", "bitscore", "spos", "sres", "next_gene_codons", "gaps_win", "qwin", "swin"]
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"{acc}.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows.values():
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    rt = {"TAA": 0, "TAG": 0, "TGA": 0}
    for r in rows.values():
        rt[r["stop"]] += 1
    stats["readthrough"] = rt
    (OUT / f"{acc}.stats.json").write_text(json.dumps(stats))
    return stats


def diamond(fasta, hits, threads):
    subprocess.run(["diamond", "blastp", "-q", str(fasta), "-d", str(DB), "-o", str(hits), "--threads", str(threads),
                    "--evalue", "1e-5", "--max-target-seqs", "5", "--max-hsps", "1", "--quiet",
                    "--outfmt", "6", *OUTFMT], check=True)


def run(acc, threads=4, rates=None):
    """rates: optional {stop: fraction of genes ending in that stop to test}; counts in 'tested' are the tested ones."""
    queries, meta, stats = prepare(acc, rates)
    with tempfile.TemporaryDirectory() as td:
        q, hits = Path(td) / "q.faa", Path(td) / "hits.tsv"
        q.write_text("".join(f">{g}\n{p}\n" for g, p in queries))
        diamond(q, hits, threads)
        return finish(acc, meta, stats, open(hits))


if __name__ == "__main__":
    for a in sys.argv[1:]:
        t = time.time()
        print(json.dumps(run(a)), f"{time.time() - t:.0f}s", flush=True)
