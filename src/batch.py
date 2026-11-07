"""Run readthrough.run over a list of accessions in parallel (skips finished ones).
usage: [KEEP_GENOMES=0] uv run python src/batch.py <acc_list.txt> [workers=4] [diamond_threads=2]"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import readthrough as rt

THREADS = 2


RATES = {"CTA": {"TAG": 1, "TAA": 0.2, "TGA": 0.2}, "TTA": {"TAA": 1, "TAG": 0.2, "TGA": 0.2},
         "TCA": {"TGA": 1, "TAA": 0.2, "TAG": 0.2}, "null": {"TAA": 0.3, "TAG": 0.3, "TGA": 0.3},
         "noRF1": {"TAG": 1, "TAA": 0.3, "TGA": 0.3}, "noRF1fast": {"TAG": 0.3, "TAA": 0.1, "TGA": 0.1}, "noRF2": {"TGA": 1, "TAA": 0.3, "TAG": 0.3},
         "noRF2fast": {"TGA": 0.3, "TAA": 0.1, "TAG": 0.1}}
CAT = {}


def one(task):
    # settings travel with each task: Pool workers are spawned fresh and don't see globals set under __main__
    acc, cat, threads = task
    if (rt.OUT / f"{acc}.stats.json").exists():
        return acc, "done", 0
    t = time.time()
    cached = (rt.GEN / f"{acc}.fna.gz").exists()
    try:
        s = rt.run(acc, threads=threads, rates=RATES.get(cat))
        return acc, json.dumps(s["readthrough"]) + " / " + json.dumps(s["tested"]), time.time() - t
    except Exception as e:  # noqa: BLE001
        return acc, f"ERROR {e}", time.time() - t
    finally:
        # KEEP_GENOMES=0: don't keep genomes this run downloaded (large screens); ones cached before are kept
        if os.environ.get("KEEP_GENOMES", "1") == "0" and not cached:
            (rt.GEN / f"{acc}.fna.gz").unlink(missing_ok=True)


if __name__ == "__main__":
    lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
    accs = [x[0].replace("GB_", "").replace("RS_", "") for x in lines]
    # optional 2nd column = screen category (CTA/TTA/TCA/null): test all genes ending in that genome's suppressor stop,
    # a fraction of the others (RATES); without it every gene is tested
    CAT.update({x[0].replace("GB_", "").replace("RS_", ""): x[1] for x in lines if len(x) > 1 and x[1] in RATES})
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    THREADS = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    with Pool(workers) as p:
        tasks = [(a, CAT.get(a), THREADS) for a in accs]
        for i, (acc, msg, dt) in enumerate(p.imap_unordered(one, tasks), 1):
            print(f"{i}/{len(accs)}\t{acc}\t{dt:.0f}s\t{msg}", flush=True)
