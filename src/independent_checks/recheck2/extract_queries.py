"""Build conserved-protein queries from NCBI RefSeq annotation of E. coli K-12 (GCF_000005845.2) and
B. subtilis 168 (GCF_000009045.1), selecting CDS by gene name (same gene list as the earlier recheck).
Also writes rf_queries.faa (prfA/prfB of both + E. coli arfB).
Usage: python3 -I extract_queries.py <recheck_dir>"""
import sys, re
R = sys.argv[1]
genes = set("""
rplA rplB rplC rplD rplE rplF rplI rplJ rplK rplL rplM rplN rplO rplP rplQ rplR rplS rplT rplU rplV rplW rplX
rpsA rpsB rpsC rpsD rpsE rpsF rpsG rpsH rpsI rpsJ rpsK rpsL rpsM rpsN rpsO rpsP rpsQ rpsR rpsS rpsT
rpmA rpmB rpmC rpmE rpmF rpmG rpmH rpmI rpmJ
alaS argS asnS aspS cysS glnS gltX glyQ glyS hisS ileS leuS lysS metG pheS pheT proS serS thrS trpS tyrS valS
tufA tuf fusA fus tsf infA infB infC frr smpB lepA prfC
rpoA rpoB rpoC rpoD rpoZ nusA nusG rho
gyrA gyrB parC parE topA recA recG recN recR dnaA dnaB dnaG dnaE dnaN dnaX holB polA ligA ruvA ruvB ruvC uvrA uvrB uvrC uvrD mutL mutS radA
ffh ftsY secA secY secE secG yidC lepB lspA
era obgE obg der engA engB yihA ychF rbfA rimM rsmA ksgA truA truB rluD mnmG gidA mnmE trmE mnmA trmD tsaD gcp tsaB rnc rnhB rnpA rbgA
groL groS dnaK dnaJ grpE clpX clpP lon ftsH htpG tig
pnp rne rnr rph
def map fmt pth ppa pgk eno tpiA gapA gap pgi fbaA pyk
adk gmk cmk pyrG pyrH tmk ndk prsA prs guaA guaB purA purB
ftsZ ftsA mraY murA murAA murB murC murD murE murF murG ddl ddlB ddlA glmS glmM glmU
atpA atpB atpC atpD atpE atpF atpG atpH
coaE dxr dxs ispE ispH ispG ispD ispF accA accB accC accD fabD fabG fabZ fabH fabI acpP plsX plsY pgsA cdsA
thrRS hemB hemL folD folC dut nrdA nrdE thyA thyX metK ileS
""".split())
rfgenes = {'prfA', 'prfB', 'arfB', 'prmC', 'hemK'}
def read_fasta(p):
    d = {}; k = None; buf = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if k: d[k] = ''.join(buf)
            k = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if k: d[k] = ''.join(buf)
    return d
n = 0; seen = set()
with open(f'{R}/data/queries.faa', 'w') as o, open(f'{R}/data/rf_queries.faa', 'w') as orf:
    for acc, org in (('GCF_000005845.2', 'Ec'), ('GCF_000009045.1', 'Bs')):
        d = f'{R}/dl/{acc}/ncbi_dataset/data/{acc}'
        prot = read_fasta(f'{d}/protein.faa')
        for line in open(f'{d}/genomic.gff'):
            f = line.split('\t')
            if len(f) < 9 or f[2] != 'CDS': continue
            m_g = re.search(r';gene=([^;]+)', f[8]); m_p = re.search(r'protein_id=([^;\n]+)', f[8])
            if not m_g or not m_p: continue
            gn, pid = m_g.group(1), m_p.group(1)
            if pid not in prot or (org, gn, pid) in seen: continue
            seen.add((org, gn, pid))
            if gn in genes:
                o.write(f'>{org}|{gn}|{pid}\n{prot[pid]}\n'); n += 1
            if gn in rfgenes:
                orf.write(f'>{org}_{gn}|{pid}\n{prot[pid]}\n')
print(n, 'queries written')
