"""Side check: do Victivallales family UBA1829 genomes (incl. the two TAG readers) carry ArfB (alternative ribosome
rescue)? Genomes must already be in the genome cache. usage: uv run python src/arfb_family.py"""
import sys, subprocess, tempfile, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from readthrough import ROOT, read_fasta, GEN
import pyrodigal
accs=[l.split('\t')[0].replace('GB_','').replace('RS_','') for l in open(ROOT/'data'/'screen'/'lists'/'uba1829.txt')]+['GCA_965332555.1']
res={}
for a in accs:
    fn=GEN/f'{a}.fna.gz'
    if not fn.exists(): continue
    c=dict(read_fasta(fn)); f=pyrodigal.GeneFinder(); f.train(*[s.encode() for s in c.values()], translation_table=11 if a not in ('GCA_949290385.1','GCA_965332555.1') else 15)
    with tempfile.TemporaryDirectory() as td:
        q=Path(td)/'q.faa'; q.write_text(''.join(f'>{k}|{g.begin}\n{g.translate()}\n' for k,s in c.items() for g in f.find_genes(s.encode()) if len(g.translate())<260))
        h=Path(td)/'h.tsv'
        subprocess.run(['diamond','blastp','-q',str(q),'-d',str(ROOT/'data'/'arfb'),'-o',str(h),'--threads','1','--quiet','--evalue','1e-10','--max-target-seqs','1','--outfmt','6','qseqid','sseqid','pident','length'],check=True)
        hits=[l.split('\t') for l in open(h)]
    arfb=[x for x in hits if '|ARFB_' in x[1]]
    res[a]=len(arfb)
    print(a, 'ArfB hits', len(arfb), flush=True)
print('genomes',len(res),'with ArfB',sum(v>0 for v in res.values()))
