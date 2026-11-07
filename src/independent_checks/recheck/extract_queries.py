"""Extract conserved-protein queries (E. coli K-12 OX=83333, B. subtilis 168 OX=224308)
from local Swiss-Prot by gene name. Usage: python -I extract_queries.py sprot.fasta.gz out.faa"""
import gzip, re, sys
src, out = sys.argv[1], sys.argv[2]
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
keep = {}
n = 0
with gzip.open(src, 'rt') as f, open(out, 'w') as o:
    hdr = None; seq = []
    def flush():
        global n
        if hdr is None: return
        m_ox = re.search(r' OX=(\d+) ', hdr); m_gn = re.search(r' GN=(\S+)', hdr)
        if not m_ox or not m_gn: return
        ox, gn = m_ox.group(1), m_gn.group(1)
        if ox not in ('83333', '224308') or gn not in genes: return
        acc = hdr.split('|')[1]
        org = 'Ec' if ox == '83333' else 'Bs'
        o.write(f'>{org}|{gn}|{acc}\n' + ''.join(seq) + '\n'); n += 1
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            flush(); hdr = line; seq = []
        else: seq.append(line)
    flush()
print(n, 'queries written')
