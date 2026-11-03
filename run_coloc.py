#!/usr/bin/env python3
"""
run_coloc.py — reproduce the single-signal colocalisation posteriors reported in
"Colocalization-gated validity audit of drug-target Mendelian randomization in
the PPi-TNAP axis" from the deposited regional summary statistics in data/.
Writes results_coloc.csv. All posteriors are COMPUTED here; no value is hard-coded,
and the figure scripts read this file rather than embedding numbers.

Engine: coloc.abf (Wakefield approximate Bayes factors), single-causal-variant
assumption; priors p1=p2=1e-4, p12=1e-5; sd.prior 0.15 (quant) / 0.20 (cc).
coloc.abf depends on each SNP's chi-square (beta/se)^2, so it is invariant to
allele orientation; SNPs are matched across datasets by position and allele set.
Cross-build harmonisation (GRCh37 exposures vs GRCh38 outcomes) uses the modal
allele-matched positional offset (see Methods).

IMPORTANT — the LPA locus is multi-signal. Single-causal coloc.abf is NOT valid
there and is NOT the basis for the positive-control verdict: Lp(a) -> CAVS is
admitted by the conditional/stepwise colocalisation in run_reanalysis.py (shared
causal variant; lead-pair r2 and conditional signal counts). The single-causal
PP.H4 for Lp(a) -> CAVS is reported below only for transparency, flagged in the
`basis` column, and must not be read as the decision variable. Lp(a) -> CAC is
rejected because the CAC outcome carries no genome-wide signal at the locus
(0 independent outcome signals in run_reanalysis.py), not on a single PP.H4.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc_abf

U = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data') + os.sep
_RESTRICTED=['LPA_Lpa_b37.tsv','LPA_CAVS_b38.tsv','axis_CAVS_b38.tsv']
_missing=[f for f in _RESTRICTED if not os.path.exists(U+f)]
if _missing:
    raise SystemExit('Missing input(s) not redistributed here: '+', '.join(_missing)+
        '\nThese derive from UK Biobank (Sinnott-Armstrong) and FinnGen R12. Download them (free) '
        'and place under data/ per data/RESTRICTED_INPUTS.md, then re-run.')

def aset(a, b):
    return '_'.join(sorted([str(a), str(b)]))

def modal_offset(p37, al37a, al37b, p38, al38a, al38b):
    L = pd.DataFrame({'p37': np.asarray(p37), 'ak': [aset(a, b) for a, b in zip(al37a, al37b)]})
    R = pd.DataFrame({'p38': np.asarray(p38), 'ak': [aset(a, b) for a, b in zip(al38a, al38b)]})
    j = pd.merge(L, R, on='ak'); j['d'] = j.p37 - j.p38
    j = j[j.d.abs() < 2_000_000]
    return int(j['d'].value_counts().idxmax()) if len(j) else 0

def run(exp, out, label, basis, type_exp='quant', type_out='cc', window=100_000):
    """exp/out: DataFrames with columns pos37, ak, beta, se. Window centred on the
    outcome lead if the outcome carries a regional signal, else the exposure lead."""
    out_min_z = (out.beta / out.se).abs().max()
    centre_df = out if out_min_z > 5.45 else exp   # |Z|>5.45 ~ p<5e-8
    lead = int(centre_df.loc[(centre_df.beta / centre_df.se).abs().idxmax(), 'pos37'])
    lo, hi = lead - window, lead + window
    E = exp[(exp.pos37 >= lo) & (exp.pos37 <= hi)][['pos37', 'ak', 'beta', 'se']]
    O = out[(out.pos37 >= lo) & (out.pos37 <= hi)][['pos37', 'ak', 'beta', 'se']]
    m = pd.merge(E, O, on=['pos37', 'ak'], suffixes=('_e', '_o'))
    m = m.replace([np.inf, -np.inf], np.nan).dropna()
    m = m[(m.se_e > 0) & (m.se_o > 0)]
    pp, n = coloc_abf(m.beta_e.values, m.se_e.values, m.beta_o.values, m.se_o.values,
                      type_exp, type_out)
    return dict(pair=label, nSNP=n, **{f'PP.{k}': round(v, 4) for k, v in pp.items()}, basis=basis)

rows = []

# ---- Lp(a) -> CAVS  (positive control; single-causal shown for transparency only) ----
lpa = pd.read_csv(U+'LPA_Lpa_b37.tsv', sep='\t')
mk = lpa['MarkerName'].str.extract(r'^(\d+):(\d+)_([ACGT]+)_([ACGT]+)$')
lpa['pos37'] = pd.to_numeric(mk[1], errors='coerce')
lpa = lpa.dropna(subset=['pos37']); lpa['pos37'] = lpa['pos37'].astype(int)
lpa['ak'] = [aset(a, b) for a, b in zip(lpa.REF, lpa.ALT)]
lpa = lpa.rename(columns={'Effect': 'beta', 'StdErr': 'se'})
cavs = pd.read_csv(U+'LPA_CAVS_b38.tsv', sep='\t')
off = modal_offset(lpa.pos37, lpa.REF, lpa.ALT, cavs.pos, cavs.ref, cavs.alt)
cavs['pos37'] = cavs['pos'] + off
cavs['ak'] = [aset(a, b) for a, b in zip(cavs.ref, cavs.alt)]
cavs = cavs.rename(columns={'sebeta': 'se'})
rows.append(run(lpa[['pos37','ak','beta','se']], cavs[['pos37','ak','beta','se']],
                'Lp(a) -> CAVS', 'single-causal only; NOT verdict basis (multi-signal locus -> conditional, run_reanalysis.py)',
                'quant', 'cc'))

# ---- Lp(a) -> CAC  (outcome carries no genome-wide signal at the locus) ----
cac_l = pd.read_csv(U+'LPA_CAC_region_b38.tsv', sep='\t')
offc = modal_offset(lpa.pos37, lpa.REF, lpa.ALT, cac_l.base_pair_location, cac_l.effect_allele, cac_l.other_allele)
cac_l['pos37'] = cac_l['base_pair_location'] + offc
cac_l['ak'] = [aset(a, b) for a, b in zip(cac_l.effect_allele, cac_l.other_allele)]
cac_l = cac_l.rename(columns={'standard_error': 'se'})
rows.append(run(lpa[['pos37','ak','beta','se']], cac_l[['pos37','ak','beta','se']],
                'Lp(a) -> CAC', 'no genome-wide CAC signal at locus (0 outcome signals, run_reanalysis.py)',
                'quant', 'quant'))

# ---- ANKH blood expression -> CAVS  (eQTLGen Z -> beta/se via MAF) ----
blood = pd.read_csv(U+'axis_eQTLGenBlood_b37.tsv', sep='\t')
ankh = blood[blood.GeneSymbol == 'ANKH'].copy(); ankh['pos37'] = ankh['SNPPos'].astype(int)
af = pd.read_csv(U+'ankh_af_region.tsv', sep='\t')
e = pd.merge(ankh, af[['SNP','AF_AlleleB','AlleleA','AlleleB']], on='SNP')
e['maf'] = np.minimum(e.AF_AlleleB.clip(0,1), 1-e.AF_AlleleB.clip(0,1)); e = e[e.maf > 0.01]
Z = e.Zscore.astype(float); Nn = e.NrSamples.astype(float)
den = 2*e.maf*(1-e.maf)*(Nn + Z**2)
e['beta'] = Z/np.sqrt(den); e['se'] = 1/np.sqrt(den)
e['ak'] = [aset(a, b) for a, b in zip(e.AlleleA, e.AlleleB)]
ac = pd.read_csv(U+'axis_CAVS_b38.tsv', sep='\t'); c5 = ac[ac.chrom == 5].copy()
o5 = modal_offset(e.pos37, e.AlleleA, e.AlleleB, c5.pos, c5.ref, c5.alt)
c5['pos37'] = c5['pos'] + o5; c5['ak'] = [aset(a, b) for a, b in zip(c5.ref, c5.alt)]
c5 = c5.rename(columns={'sebeta':'se'})
rows.append(run(e[['pos37','ak','beta','se']], c5[['pos37','ak','beta','se']],
                'ANKH blood -> CAVS', 'single-signal coloc.abf', 'quant', 'cc'))

# ---- ENPP1 aortic expression -> CAC  (both GRCh38) ----
aorta = pd.read_csv(U+'axis_GTExAorta_b38.tsv', sep='\t')
aorta['gene'] = aorta['gene_id'].str.split('.').str[0]
ae = aorta[aorta.gene == 'ENSG00000197594'].copy()
parts = ae['variant_id'].str.split('_', expand=True)
ae['pos37'] = parts[1].astype(int); ae['ak'] = [aset(a, b) for a, b in zip(parts[2], parts[3])]
ae = ae.rename(columns={'slope':'beta','slope_se':'se'})
cac = pd.read_csv(U+'axis_CAC_b38.tsv', sep='\t'); c6 = cac[cac.chromosome == 6].copy()
c6['pos37'] = c6['base_pair_location']; c6['ak'] = [aset(a, b) for a, b in zip(c6.effect_allele, c6.other_allele)]
c6 = c6.rename(columns={'standard_error':'se'})
rows.append(run(ae[['pos37','ak','beta','se']], c6[['pos37','ak','beta','se']],
                'ENPP1 aorta -> CAC', 'single-signal coloc.abf', 'quant', 'quant'))

# ---- ENPP1 blood expression -> CAC  (eQTLGen Z -> beta/se; MAF from CAC EAF by rsid) ----
en = blood[blood.GeneSymbol == 'ENPP1'].copy(); en['pos37'] = en['SNPPos'].astype(int)
mafmap = dict(zip(c6.rsid, c6.effect_allele_frequency))
en['eaf'] = en.SNP.map(mafmap); en = en.dropna(subset=['eaf'])
en['maf'] = np.minimum(en.eaf.clip(0,1), 1-en.eaf.clip(0,1)); en = en[en.maf > 0.01]
Ze = en.Zscore.astype(float); Ne = en.NrSamples.astype(float)
dene = 2*en.maf*(1-en.maf)*(Ne + Ze**2)
en['beta'] = Ze/np.sqrt(dene); en['se'] = 1/np.sqrt(dene)
en['ak'] = [aset(a, b) for a, b in zip(en.AssessedAllele, en.OtherAllele)]
oe = modal_offset(en.pos37, en.AssessedAllele, en.OtherAllele, c6.base_pair_location, c6.effect_allele, c6.other_allele)
c6b = c6.copy(); c6b['pos37'] = c6b['base_pair_location'] + oe
rows.append(run(en[['pos37','ak','beta','se']], c6b[['pos37','ak','beta','se']],
                'ENPP1 blood -> CAC', 'single-signal coloc.abf', 'quant', 'quant'))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_coloc.csv'), index=False)
pd.set_option('display.width', 200); pd.set_option('display.max_colwidth', 80)
print(res.to_string(index=False))
