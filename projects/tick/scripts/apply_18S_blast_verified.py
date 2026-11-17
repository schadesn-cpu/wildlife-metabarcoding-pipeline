#!/usr/bin/env python3
"""
apply_18S_blast_verified.py
Build results/18S/all/taxonomy_verified/ from the raw 18S classifier export + the BLAST
corrections verified on 2026-07-09 (Verification_Batch_Results.md). Matches the 16S pattern:
  - taxonomy_18S_blast_verified.tsv   (Feature ID / Taxon / Confidence) — corrected taxonomy
  - blast_corrections_applied.tsv     (feature_id / old_taxon / new_taxon / note / confidence / total_reads)

Two kinds of correction:
  1. TARGETED per-ASV corrections for the 10 ASVs BLAST-verified today (species-level pathogen calls,
     tick/host/ked IDs, and the cross-amplification flag on the 23M-read "Rickettsia").
  2. BLANKET flag: EVERY ASV whose raw taxonomy is Bacteria is bacterial cross-amplification in an
     18S run -> relabel with a CROSSAMP prefix so it can be filtered from the eukaryotic community.
     (18S primers should not amplify bacteria; these reads are the endosymbiont 16S bleeding through.)

Nothing is deleted; the feature table is untouched. This only produces a corrected taxonomy the
downstream 18S analyses should read INSTEAD of the raw export.
"""
import os, glob, zipfile, tempfile
import pandas as pd
import biom

RAW = "qiime2/18S/all/taxonomy-exported-temp/taxonomy.tsv"
OUTDIR = "results/18S/all/taxonomy_verified"
os.makedirs(OUTDIR, exist_ok=True)

# ---- targeted corrections (prefix -> (new_taxon, note)) ----
# new_taxon uses a readable, source-consistent string; note explains the BLAST basis.
TARGETED = {
 'bbf858d8': ("Eukaryota;Opisthokonta;Metazoa;Arthropoda;Arachnida;Ixodida;Ixodidae;Ixodes;Ixodes_scapularis",
              "BLAST 100% I. scapularis SSU rRNA (tick host)"),
 '87a7cc32': ("Eukaryota;Opisthokonta;Metazoa;Arthropoda;Arachnida;Ixodida;Ixodidae;Dermacentor",
              "BLAST 100% Dermacentor SSU rRNA (winter tick; genus-level)"),
 'de15f073': ("CROSSAMP;Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Rickettsia",
              "BLAST 100% Rickettsia 16S rRNA — BACTERIAL CROSS-AMPLIFICATION, exclude from eukaryotic community"),
 '88104f95': ("Eukaryota;Opisthokonta;Metazoa;Arthropoda;Arachnida;Ixodida;Ixodidae;Ixodes;Ixodes_scapularis_genomic",
              "BLAST 89.5% I. scapularis BAC — tick genomic DNA, not pathogen"),
 'b174eed7': ("Eukaryota;Opisthokonta;Metazoa;Chordata;Mammalia;Artiodactyla;Cervidae;Odocoileus;Odocoileus_virginianus",
              "BLAST 100% O. virginianus (deer blood-meal host DNA)"),
 '14bc5ef0': ("Eukaryota;Opisthokonta;Metazoa;Arthropoda;Insecta;Diptera;Hippoboscidae;Lipoptena;Lipoptena_cervi",
              "BLAST 100% Lipoptena cervi (deer ked)"),
 'cc42bc05': ("Eukaryota;TSAR;Alveolata;Apicomplexa;Aconoidasida;Piroplasmida;Babesiidae;Babesia;Babesia_odocoilei",
              "BLAST 100% Babesia odocoilei (deer parasite)"),
 'dbba6955': ("Eukaryota;TSAR;Alveolata;Apicomplexa;Aconoidasida;Piroplasmida;Babesiidae;Babesia;Babesia_microti",
              "BLAST 100% Babesia microti (human pathogen)"),
 '506c19f6': ("Eukaryota;TSAR;Alveolata;Apicomplexa;Aconoidasida;Piroplasmida;Theileriidae;Theileria;Theileria_cervi",
              "BLAST 100% Theileria cervi (deer piroplasm)"),
 '1e50fc95': ("Eukaryota;TSAR;Alveolata;Apicomplexa;Coccidia;Eucoccidiorida;Hepatozoidae;Hepatozoon",
              "BLAST 100% Hepatozoon sp. (apicomplexan)"),
}

def load_table_reads():
    q=[x for x in glob.glob('qiime2/18S/all/**/table*.qza',recursive=True) if 'raref' not in x][0]
    td=tempfile.mkdtemp()
    with zipfile.ZipFile(q) as z: z.extractall(td)
    t=biom.load_table(glob.glob(td+'/*/data/feature-table.biom')[0]).to_dataframe(dense=True)
    return t.sum(axis=1)  # reads per ASV

def main():
    tax=pd.read_csv(RAW, sep='\t')
    tax.columns=['Feature ID','Taxon','Confidence'][:len(tax.columns)]
    reads=load_table_reads()

    corrections=[]
    new_taxon=tax['Taxon'].copy()

    # 1. targeted per-ASV corrections
    for i,row in tax.iterrows():
        fid=row['Feature ID']; pref=fid[:8]
        if pref in TARGETED:
            nt,note=TARGETED[pref]
            corrections.append({'feature_id':fid,'old_taxon':row['Taxon'],'new_taxon':nt,
                                'note':note,'confidence':row['Confidence'],
                                'total_reads':int(reads.get(fid,0))})
            new_taxon.at[i]=nt

    # 2. blanket: flag ALL Bacteria-labeled ASVs as cross-amplification (if not already targeted)
    targeted_ids={c['feature_id'] for c in corrections}
    for i,row in tax.iterrows():
        fid=row['Feature ID']
        if fid in targeted_ids: continue
        t=str(row['Taxon'])
        if t.startswith('Bacteria') or ';Bacteria' in t or t.startswith('d__Bacteria'):
            nt='CROSSAMP;'+t
            corrections.append({'feature_id':fid,'old_taxon':t,'new_taxon':nt,
                                'note':'Bacteria in 18S = cross-amplification; exclude from eukaryotic community',
                                'confidence':row['Confidence'],'total_reads':int(reads.get(fid,0))})
            new_taxon.at[i]=nt

    # write verified taxonomy (same 3-col format as 16S)
    out=tax.copy(); out['Taxon']=new_taxon
    out.to_csv(f'{OUTDIR}/taxonomy_18S_blast_verified.tsv', sep='\t', index=False)

    cdf=pd.DataFrame(corrections).sort_values('total_reads',ascending=False)
    cdf.to_csv(f'{OUTDIR}/blast_corrections_applied.tsv', sep='\t', index=False)

    # summary
    n_targeted=len(TARGETED)
    n_crossamp=sum(1 for c in corrections if 'cross-amplification' in c['note'])
    reads_crossamp=sum(c['total_reads'] for c in corrections if 'cross-amplification' in c['note'] or 'CROSS-AMPL' in c['note'])
    reads_total=int(reads.sum())
    print("wrote", f'{OUTDIR}/taxonomy_18S_blast_verified.tsv')
    print("wrote", f'{OUTDIR}/blast_corrections_applied.tsv')
    print(f"\ntargeted per-ASV corrections applied: {len([c for c in corrections if c['feature_id'][:8] in TARGETED])}/{n_targeted}")
    print(f"bacterial cross-amplification ASVs flagged: {n_crossamp}")
    print(f"reads flagged as cross-amplification: {reads_crossamp:,} of {reads_total:,} total 18S reads ({100*reads_crossamp/reads_total:.1f}%)")
    print("\nTop cross-amplification ASVs by reads:")
    ca=cdf[cdf['note'].str.contains('cross-ampl',case=False)].head(6)
    for _,r in ca.iterrows():
        print(f"  {r['feature_id'][:10]}  {r['total_reads']:>12,}  {r['old_taxon'][:55]}")

if __name__=='__main__':
    main()
