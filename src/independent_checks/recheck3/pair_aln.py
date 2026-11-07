"""Print local pairwise alignment (BLOSUM62) of all pairs of sequences in a fasta. Usage: python -I pair_aln.py file.faa"""
import sys
from Bio import Align
from Bio.Align import substitution_matrices
seqs = []; k = None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): seqs.append([l[1:], ''])
    elif l: seqs[-1][1] += l
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62'); al.open_gap_score = -11; al.extend_gap_score = -1; al.mode = 'local'
for i in range(len(seqs)):
    for j in range(i + 1, len(seqs)):
        a = al.align(seqs[i][1], seqs[j][1])[0]
        print('#', seqs[i][0], 'vs', seqs[j][0], 'score', a.score); print(a)
