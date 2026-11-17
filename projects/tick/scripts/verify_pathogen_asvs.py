#!/usr/bin/env python3
"""
verify_pathogen_asvs.py  —  run ON RON, from project root:
    cd ~/meedlab/tick_project
    python3 verify_pathogen_asvs.py --prep      # step 1: pull pathogen ASV seqs -> FASTA
    #   then run blastn (command is printed for you), then:
    python3 verify_pathogen_asvs.py --report    # step 2: parse hits, flag SILVA vs BLAST disagreements

WHY: the 16S BLAST audit only covered 633 of 89,971 ASVs (top reads). Four of five bacterial
pathogens (Borreliella, Borrelia, Bartonella, Ehrlichia) had ZERO ASVs BLAST-verified — they ride
on the SILVA classifier label. This verifies ONLY the pathogen-labeled ASVs (~150) so the pathogen
calls rest on BLAST, not just the classifier. (Anaplasma mostly done already, included for completeness.)

Reads:
    results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv   (feature_id -> taxon)
    <rep-seqs FASTA>   (auto-located; the DADA2 representative sequences)
Writes:
    results/16S/all/blast_verify/pathogen_asvs.fasta          (step 1)
    results/16S/all/blast_verify/pathogen_asvs_blast.tsv      (you run blastn -> this)
    results/16S/all/blast_verify/pathogen_verification_report.tsv (step 2)
"""
import sys, os, glob, subprocess, pandas as pd

PATHOGENS = ['Borreliella','Borrelia','Bartonella','Ehrlichia','Anaplasma']
TAX   = 'results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv'
OUTD  = 'results/16S/all/blast_verify'
FASTA = f'{OUTD}/pathogen_asvs.fasta'
HITS  = f'{OUTD}/pathogen_asvs_blast.tsv'
REPORT= f'{OUTD}/pathogen_verification_report.tsv'
os.makedirs(OUTD, exist_ok=True)

def find_repseqs():
    if os.path.exists('qiime2/16S/all/rep-seqs-exported-temp/dna-sequences.fasta'):
        return 'qiime2/16S/all/rep-seqs-exported-temp/dna-sequences.fasta'
    # common locations / names for DADA2 rep-seqs
    pats = ['results/16S/**/*rep*seq*.fasta','results/16S/**/*rep*seq*.fa',
            'results/16S/**/dna-sequences.fasta','**/16S*rep*seqs*/**/*.fasta']
    for p in pats:
        hits=glob.glob(p, recursive=True)
        if hits: return hits[0]
    return None

def prep():
    tax = pd.read_csv(TAX, sep='\t')
    # first col = feature_id, some col holds the taxon string
    idcol = tax.columns[0]
    taxcol = [c for c in tax.columns if 'tax' in c.lower()] or [tax.columns[1]]
    taxcol = taxcol[0]
    mask = tax[taxcol].astype(str).str.contains('|'.join(PATHOGENS), case=False, na=False)
    want = set(tax.loc[mask, idcol].astype(str))
    print(f"pathogen-labeled ASVs to verify: {len(want)}")
    for p in PATHOGENS:
        n = tax[taxcol].astype(str).str.contains(p, case=False, na=False).sum()
        print(f"  {p}: {n}")

    rs = find_repseqs()
    if not rs:
        print("!! could not auto-find rep-seqs FASTA. Set it manually:")
        print("   qiime tools export --input-path <rep-seqs.qza> --output-path", OUTD)
        sys.exit(1)
    print("rep-seqs:", rs)
    # pull the wanted sequences
    keep=False; n=0
    with open(rs) as fh, open(FASTA,'w') as out:
        for line in fh:
            if line.startswith('>'):
                fid=line[1:].split()[0].strip()
                keep = fid in want
                if keep: out.write(line); n+=1
            elif keep:
                out.write(line)
    print(f"wrote {n} sequences -> {FASTA}")
    print("\nNOW RUN BLAST (local nt):")
    print(f"  blastn -query {FASTA} -db nt -out {HITS} -max_target_seqs 5 -evalue 1e-10 \\")
    print(f"     -outfmt '6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle' -num_threads 8")
    print(f"\nthen: python3 {sys.argv[0]} --report")

def report():
    tax = pd.read_csv(TAX, sep='\t')
    idcol=tax.columns[0]; taxcol=[c for c in tax.columns if 'tax' in c.lower()][0]
    silva = dict(zip(tax[idcol].astype(str), tax[taxcol].astype(str)))
    cols=['qseqid','sseqid','pident','length','evalue','bitscore','staxids','sscinames','stitle']
    h = pd.read_csv(HITS, sep='\t', names=cols)
    # top hit per ASV (highest bitscore)
    top = h.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    rows=[]
    for _,r in top.iterrows():
        asv=r['qseqid']; silva_tax=silva.get(asv,'?')
        silva_genus=[p for p in PATHOGENS if p.lower() in silva_tax.lower()]
        silva_genus=silva_genus[0] if silva_genus else silva_tax.split(';')[-1]
        blast_name=str(r['sscinames']) if pd.notna(r['sscinames']) else str(r['stitle'])[:40]
        agree = silva_genus.lower() in blast_name.lower() or any(
            g.lower() in blast_name.lower() for g in [silva_genus])
        rows.append({'ASV':asv,'SILVA_genus':silva_genus,'BLAST_top_hit':blast_name,
                     'pident':r['pident'],'AGREE':'yes' if agree else 'NO — CHECK'})
    rep=pd.DataFrame(rows).sort_values(['AGREE','SILVA_genus'])
    rep.to_csv(REPORT, sep='\t', index=False)
    print(f"wrote {REPORT}\n")
    print("=== DISAGREEMENTS (SILVA label != BLAST top hit) — these need your eye ===")
    dis=rep[rep['AGREE']!='yes']
    print(dis.to_string(index=False) if len(dis) else "  none — every pathogen ASV's BLAST top hit matches its SILVA genus ✓")
    print(f"\nsummary: {len(rep)} ASVs checked, {len(dis)} disagreements")
    print("per pathogen agreement:")
    print(rep.groupby('SILVA_genus')['AGREE'].value_counts().to_string())

if __name__=='__main__':
    if '--prep' in sys.argv: prep()
    elif '--report' in sys.argv: report()
    else: print("usage: --prep  (then run blastn)  then  --report")
