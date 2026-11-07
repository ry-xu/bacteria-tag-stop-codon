"""Bulk-download genome FASTAs from NCBI Datasets (batches of 100 accessions per request) into one gzipped FASTA per
genome under <outdir>/. Sequential requests (NCBI, not EBI); skips files already present.
usage: python src/bulk_fetch.py <acc_list> <outdir>"""
import gzip
import io
import sys
import time
import zipfile
from pathlib import Path

import requests

accs = [l.split()[0] for l in open(sys.argv[1]) if l.strip()]
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
todo = [a for a in accs if not (out / f"{a}.fna.gz").exists()]
for i in range(0, len(todo), 100):
    batch = todo[i:i + 100]
    for attempt in range(5):
        try:
            r = requests.post("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download",
                              json={"accessions": batch, "include_annotation_type": ["GENOME_FASTA"]}, timeout=600)
            if r.status_code == 200:
                break
        except requests.RequestException:
            pass
        time.sleep(10 * (attempt + 1))
    else:
        print("failed batch", i, flush=True)
        continue
    z = zipfile.ZipFile(io.BytesIO(r.content))
    for n in z.namelist():
        if n.endswith(".fna"):
            acc = n.split("/")[-2]
            (out / f"{acc}.fna.gz").write_bytes(gzip.compress(z.read(n)))
    print(f"{min(i + 100, len(todo))}/{len(todo)}", flush=True)
