"""Like batch.py, but one DIAMOND search per POOL genomes instead of one per genome: workers download and call genes
(readthrough.prepare), the main process searches their pooled queries against Swiss-Prot and writes each genome's
results (readthrough.finish). Same outputs as batch.py; several times less CPU on small subsampled query sets, where
a per-genome DIAMOND run is mostly the fixed cost of scanning Swiss-Prot.

usage: [KEEP_GENOMES=0] uv run python src/batch_pooled.py <acc_list> \
           [workers=4] [diamond_threads=4] [pool=40]
"""
import os
import sys
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import readthrough as rt
from batch import RATES


def prep(task):
    acc, cat = task
    t = time.time()
    cached = (rt.GEN / f"{acc}.fna.gz").exists()
    try:
        return acc, rt.prepare(acc, RATES.get(cat)), time.time() - t
    except Exception as e:  # noqa: BLE001
        return acc, f"ERROR {e}", time.time() - t
    finally:
        if os.environ.get("KEEP_GENOMES", "1") == "0" and not cached:
            (rt.GEN / f"{acc}.fna.gz").unlink(missing_ok=True)


def search(buf, threads, log):
    t = time.time()
    with tempfile.TemporaryDirectory() as td:
        q, hits = Path(td) / "q.faa", Path(td) / "hits.tsv"
        with open(q, "w") as f:
            for k, (acc, (queries, meta, stats)) in enumerate(buf):
                for g, p in queries:
                    f.write(f">{k}#{g}\n{p}\n")
        rt.diamond(q, hits, threads)
        per = [[] for _ in buf]
        for line in open(hits):
            k, rest = line.split("#", 1)
            per[int(k)].append(rest)
    for (acc, (queries, meta, stats)), lines in zip(buf, per):
        s = rt.finish(acc, meta, stats, lines)
        log(f"{acc}\t{s['readthrough']} / {s['tested']}")
    log(f"# searched {len(buf)} genomes, {sum(len(x[1][0]) for x in buf)} queries in {time.time() - t:.0f}s")


if __name__ == "__main__":
    lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    pool = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    tasks = [(x[0].replace("GB_", "").replace("RS_", ""), x[1] if len(x) > 1 else None) for x in lines]
    tasks = [t for t in tasks if not (rt.OUT / f"{t[0]}.stats.json").exists()]
    n = 0

    def log(msg):
        print(f"{time.strftime('%H:%M:%S')}\t{msg}", flush=True)

    log(f"# {len(tasks)} genomes to do ({len(lines) - len(tasks)} already done)")
    buf = []
    with Pool(workers) as p:
        for acc, res, dt in p.imap_unordered(prep, tasks):
            n += 1
            if isinstance(res, str):
                log(f"{acc}\t{res}")
                continue
            buf.append((acc, res))
            if len(buf) >= pool:
                search(buf, threads, log)
                log(f"# {n}/{len(tasks)}")
                buf = []
        if buf:
            search(buf, threads, log)
    log(f"done: {n}/{len(tasks)}")
