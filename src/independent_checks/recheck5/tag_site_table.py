"""Per-site table of T's TAG codons that pass the 'coding' filter (strict + control codon in-frame in a control gene) in >= 1
control: control codons/residues, outgroup residues (loose), T gene (pyrodigal call containing the codon in-frame, if any) and its
DIAMOND top Swiss-Prot hit; the T ORF (TAG read as X) is used for annotation when no gene call contains the codon.
Usage: python3 -I tag_site_table.py <dir> T ctrl1,ctrl2,... outg1,outg2,..."""
import sys, subprocess, collections
sys.argv = sys.argv[:4] + [sys.argv[4], 'coding']
src = open(sys.argv[0].replace('tag_site_table.py', 'tag_meaning5.py')).read().split("print(f'##### {T}")[0]
exec(src)
ot = {O: load(O, strict=False) for O in outg}
sites = sorted({k for O in ctrls for k, f in tabs[O].items() if f[3] == 'TAG'})
g = read_contigs(T); COMPL2 = str.maketrans('ACGT', 'TGCA')
genes = collections.defaultdict(list)
for l in open(f'{R}/out/genes/{T}.genes.tsv'):
    c, b_, e_, st = l.rstrip('\n').split('\t'); genes[c].append((int(b_), int(e_), int(st)))
def orf_of(c, g0, st):
    s = g[c]
    if st == 1:
        a = g0
        while a - 3 >= 0 and CODE.get(s[a-3:a]) != '*' or (a - 3 >= 0 and s[a-3:a] == 'TAG'): a -= 3
        z = g0
        while z + 3 <= len(s) and (CODE.get(s[z:z+3]) != '*' or s[z:z+3] == 'TAG'): z += 3
        seq = s[a:z]
    else:
        rcs = s[::-1].translate(COMPL2); o = len(s) - g0 - 3
        a = o
        while a - 3 >= 0 and (CODE.get(rcs[a-3:a]) != '*' or rcs[a-3:a] == 'TAG'): a -= 3
        z = o
        while z + 3 <= len(rcs) and (CODE.get(rcs[z:z+3]) != '*' or rcs[z:z+3] == 'TAG'): z += 3
        seq = rcs[a:z]
    return ''.join('X' if seq[i:i+3] == 'TAG' else CODE.get(seq[i:i+3], 'X') for i in range(0, len(seq) - 2, 3))
fa = f'{R}/out/tag_sites_orfs.faa'
with open(fa, 'w') as o:
    for k in sites: o.write(f'>{k[0]}:{k[1]}:{k[2]}\n{orf_of(k[0], int(k[1]), int(k[2]))}\n')
out = subprocess.run(['diamond', 'blastp', '-q', fa, '-d', f'{R}/data/sprot', '--threads', '2', '--quiet', '-k', '1', '-e', '1e-5',
                      '--outfmt', '6', 'qseqid', 'pident', 'stitle'], capture_output=True, text=True, check=True).stdout
ann = {}
for l in out.splitlines():
    q, pid, st = l.split('\t'); ann.setdefault(q, f'{st.split(" OS=")[0].split(" ", 1)[1][:42]} ({float(pid):.0f}%)')
def in_gene(c, g0, st):
    for b_, e_, s_ in genes[c]:
        if s_ != st or not (b_ - 1 <= g0 and g0 + 3 <= e_): continue
        if (st == 1 and (g0 - (b_ - 1)) % 3 == 0) or (st == -1 and (e_ - (g0 + 3)) % 3 == 0):
            return 'internal?' if False else ('gene-end' if (st == 1 and g0 + 3 == e_) or (st == -1 and g0 == b_ - 1) else 'in-gene')
    return 'no-call'
print('site\tcontrols(codon:res)\tconsensus\toutgroups\tT_gene_call\tannotation(ORF, TAG=X)')
for k in sites:
    cr = [f'{O[4:13]}:{tabs[O][k][6]}({tabs[O][k][7]})' for O in ctrls if k in tabs[O]]
    rs = [tabs[O][k][7] for O in ctrls if k in tabs[O]]; cons_ = rs[0] if len(set(rs)) == 1 else 'mixed'
    og = ' '.join(f'{O[4:13]}={ot[O][k][7] if k in ot[O] else "."}' for O in outg)
    key = f'{k[0]}:{k[1]}:{k[2]}'
    print(f'{key}\t{" ".join(cr)}\t{cons_}{"" if len(rs) > 1 else "(1 ctrl)"}\t{og}\t{in_gene(k[0], int(k[1]), int(k[2]))}\t{ann.get(key, "no Swiss-Prot hit")}')
