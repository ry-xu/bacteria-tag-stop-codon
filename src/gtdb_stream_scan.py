"""Stream GTDB R232's per-genome gene files (gtdb_proteins_nt_reps_r232.tar.gz, 187 GB, never written to disk) and
record, for every species representative:
  - stop-codon usage of complete genes (GTDB's Prodigal calls; last codon of each gene)
  - stop-decoding machinery: proteins (200-800 aa) whose best DIAMOND hit (e<=1e-20) in data/markers.dmnd is a
    marker: RF1, RF2 (bacterial release factors; RF1 reads UAA/UAG, RF2 reads UAA/UGA, so losing one frees its stop
    for reassignment), aRF1, PylS/B/C/D (pyrrolysine), SelA/SelB/SelD/SepSecS/PSTK (selenocysteine), SepRS.
    The database also holds decoy families (EF-Tu, IF2, PheRS, LysRS, radical-SAM, ...) so look-alikes hit those.
Archaea use aRF1 for all stops, so their RF1/RF2 columns are not meaningful.

usage: nice -n 10 uv run python src/gtdb_stream_scan.py [max_genomes]
       -> results/gtdb_stops_rf.tsv (appends; restart skips genomes already written)
"""
import gzip
import subprocess
import sys
import tarfile
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
URL = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_nt_reps_r232.tar.gz"
OUT = ROOT / "results" / "gtdb_stops_rf.tsv"
RFDB = ROOT / "data" / "markers.dmnd"
MARKERS = ["aRF1", "PylS", "PylB", "PylC", "PylD", "SelA", "SelB", "SelD", "SepSecS", "PSTK", "SepRS"]
BATCH = 100
COLS = ["acc", "domain", "genes", "complete", "TAA", "TAG", "TGA", "other_end", "n_RF1", "n_RF2", "RF1_best", "RF2_best", "markers"]
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODE[a + b + c] = _aa[16 * i + 4 * j + k]


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3)).rstrip("*").replace("*", "X")


def parse(raw):
    """Yield (header, seq) from a gzipped FASTA held in memory."""
    name, buf = None, []
    for line in gzip.decompress(raw).decode("ascii", "replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf).upper()
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()


def genome_stats(raw):
    stops, n, comp, other = Counter(), 0, 0, 0
    prots = []
    for h, s in parse(raw):
        n += 1
        f = h.split(" # ")
        strand = f[3].strip() if len(f) > 3 else "1"
        part = h.split("partial=")[1][:2] if "partial=" in h else "00"
        complete_stop = (strand == "1" and part[1] == "0") or (strand == "-1" and part[0] == "0")
        if complete_stop:
            comp += 1
            c = s[-3:]
            if c in ("TAA", "TAG", "TGA"):
                stops[c] += 1
            else:
                other += 1
        L = len(s) // 3
        if 200 <= L <= 800:  # marker proteins are ~250-700 aa; skip the rest to keep DIAMOND fast
            prots.append((f[0].split()[0], translate(s)))
    return n, comp, stops, other, prots


def run_diamond(batch):
    """batch: list of (acc, prots). Returns {acc: {'RF1': [(pid, len)], 'RF2': [...]}}."""
    res = {acc: {k: [] for k in ["RF1", "RF2"] + MARKERS} for acc, _ in batch}
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        with open(q, "w") as f:
            for acc, prots in batch:
                for pid, seq in prots:
                    f.write(f">{acc}|{pid}\n{seq}\n")
        h = Path(td) / "h.tsv"
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(RFDB), "-o", str(h), "--threads", "4", "--quiet",
                        "--evalue", "1e-20", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "sseqid", "pident",
                        "length"], check=True)
        for line in open(h):
            qid, sid, pid, ln = line.rstrip("\n").split("\t")
            acc = qid.split("|")[0]
            kind = sid.split("|")[0]
            if kind in res[acc]:
                res[acc][kind].append((float(pid), int(ln)))
    return res


def flush(batch, out):
    rf = run_diamond([(acc, st[4]) for acc, dom, st in batch])
    for acc, dom, (n, comp, stops, other, _) in batch:
        r1 = sorted(rf[acc]["RF1"], reverse=True)
        r2 = sorted(rf[acc]["RF2"], reverse=True)
        row = [acc, dom, n, comp, stops["TAA"], stops["TAG"], stops["TGA"], other, len(r1), len(r2),
               f"{r1[0][0]}/{r1[0][1]}" if r1 else "", f"{r2[0][0]}/{r2[0][1]}" if r2 else "",
               ";".join(f"{k}:{len(rf[acc][k])}({max(rf[acc][k])[0]})" for k in MARKERS if rf[acc][k])]
        out.write("\t".join(map(str, row)) + "\n")
    out.flush()


def main(limit=None):
    done = set()
    if OUT.exists():
        done = {l.split("\t")[0] for l in open(OUT)}
    out = open(OUT, "a")
    if not done:
        out.write("\t".join(COLS) + "\n")
    t0, seen, batch = time.time(), 0, []
    # DIAMOND runs in one background thread so the download never stalls while a batch is searched (it did ~40% of
    # the time); at most 2 batches queued, writes stay in order (single worker)
    pool, pending = ThreadPoolExecutor(max_workers=1), []
    with requests.get(URL, stream=True, timeout=600) as r:
        r.raise_for_status()
        tar = tarfile.open(fileobj=r.raw, mode="r|gz")
        for m in tar:
            if not m.isfile() or not m.name.endswith("_protein.fna.gz"):
                continue
            acc = Path(m.name).name.replace("_protein.fna.gz", "")
            dom = "archaea" if "/archaea/" in m.name else "bacteria"
            raw = tar.extractfile(m).read()
            seen += 1
            if acc in done:
                continue
            batch.append((acc, dom, genome_stats(raw)))
            if len(batch) >= BATCH:
                pending = [f for f in pending if not f.done()]
                while len(pending) >= 2:
                    pending[0].result()
                    pending = [f for f in pending if not f.done()]
                pending.append(pool.submit(flush, batch, out))
                batch = []
                print(f"{seen} genomes, {time.time() - t0:.0f}s", flush=True)
            if limit and seen >= limit:
                break
    if batch:
        pending.append(pool.submit(flush, batch, out))
    for f in pending:
        f.result()
    pool.shutdown()
    print(f"done: {seen} genomes, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
