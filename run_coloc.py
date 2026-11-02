#!/usr/bin/env python3
"""
run_coloc.py — reproduce the principal colocalisation posteriors reported in
"Colocalization-gated validity audit of drug-target Mendelian randomization in
the PPi-TNAP axis" from the deposited regional summary statistics in data/.

Engine: coloc.abf (Wakefield approximate Bayes factors), single-causal-variant
assumption; priors p1=p2=1e-4, p12=1e-5; sd.prior 0.15 (quant) / 0.20 (cc).
coloc.abf depends on each SNP's chi-square (beta/se)^2, so it is invariant to
allele orientation; SNPs are matched across datasets by position and allele set.

Cross-build harmonisation (GRCh37 exposures vs GRCh38 outcomes) uses a constant
positional offset estimated as the modal base-pair difference among allele-matched
variants shared between the two datasets (as described in Methods).

Running this reproduces, from the deposited data:
  Lp(a) -> CAVS         PP.H4 ~ 0.998  (positive control, PASS strong)
  ANKH blood -> CAVS    PP.H4 ~ 0.27   PP.H1 ~ 0.71  (non-colocalizing nominal)
  ENPP1 aorta -> CAC    PP.H4 ~ 0.05   (fail-coloc)
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc_abf

U = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data') + os.sep

def aset(a, b):
    return '_'.join(sorted([str(a), str(b)]))

def modal_offset(p37, al37a, al37b, p38, al38a, al38b):
    L = pd.DataFrame({'p37': p37, 'ak': [aset(a, b) for a, b in zip(al37a, al37b)]})
    R = pd.DataFrame({'p38': p38, 'ak': [aset(a, b) for a, b in zip(al38a, al38b)]})
    j = pd.merge(L, R, on='ak'); j['d'] = j.p37 - j.p38
    j = j[j.d.abs() < 2_000_000]
    return int(j['d'].value_counts().idxmax())

def run(exp, out, label, type_exp='quant', type_out='cc', window=100_000):
    """exp/out: DataFrames with columns pos37, ak, beta, se. Window centred on
    the outcome lead if the outcome carries a regional signal, else the exposure lead."""
    # choose centre: outcome lead if it is genome-wide-ish, else exposure lead
    out_min_p = (out.beta / out.se).abs().max()
    centre_df = out if out_min_p > 5.45 else exp   # |Z|>5.45 ~ p<5e-8
    lead = int(centre_df.loc[(centre_df.beta / centre_df.se).abs().idxmax(), 'pos37'])
    lo, hi = lead - window, lead + window
    E = exp[(exp.pos37 >= lo) & (exp.pos37 <= hi)][['pos37', 'ak', 'beta', 'se']]
    O = out[(out.pos37 >= lo) & (out.pos37 <= hi)][['pos37', 'ak', 'beta', 'se']]
    m = pd.merge(E, O, on=['pos37', 'ak'], suffixes=('_e', '_o'))
    m = m.replace([np.inf, -np.inf], np.nan).dropna()
    m = m[(m.se_e > 0) & (m.se_o > 0)]
    pp, n = coloc_abf(m.beta_e.values, m.se_e.values, m.beta_o.values, m.se_o.values,
                      type_exp, type_out)
    return dict(pair=label, nSNP=n, **{f'PP.{k}': round(v, 4) for k, v in pp.items()})

rows = []

# ---- 1) Lp(a) -> CAVS  (positive control) ----
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
cavs = cavs.rename(columns={'beta': 'beta', 'sebeta': 'se'})
rows.append(run(lpa[['pos37','ak','beta','se']], cavs[['pos37','ak','beta','se']],
                'Lp(a) -> CAVS', 'quant', 'cc'))

# ---- 2) ANKH blood expression -> CAVS  (eQTLGen Z -> beta/se via MAF) ----
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
c5 = c5.rename(columns={'beta':'beta','sebeta':'se'})
rows.append(run(e[['pos37','ak','beta','se']], c5[['pos37','ak','beta','se']],
                'ANKH blood -> CAVS', 'quant', 'cc'))

# ---- 3) ENPP1 aortic expression -> CAC  (both GRCh38) ----
aorta = pd.read_csv(U+'axis_GTExAorta_b38.tsv', sep='\t')
aorta['gene'] = aorta['gene_id'].str.split('.').str[0]
ae = aorta[aorta.gene == 'ENSG00000197594'].copy()
parts = ae['variant_id'].str.split('_', expand=True)
ae['pos37'] = parts[1].astype(int); ae['ak'] = [aset(a, b) for a, b in zip(parts[2], parts[3])]
ae = ae.rename(columns={'slope':'beta','slope_se':'se'})
cac = pd.read_csv(U+'axis_CAC_b38.tsv', sep='\t'); c6 = cac[cac.chromosome == 6].copy()
c6['pos37'] = c6['base_pair_location']; c6['ak'] = [aset(a, b) for a, b in zip(c6.effect_allele, c6.other_allele)]
c6 = c6.rename(columns={'beta':'beta','standard_error':'se'})
rows.append(run(ae[['pos37','ak','beta','se']], c6[['pos37','ak','beta','se']],
                'ENPP1 aorta -> CAC', 'quant', 'cc'))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_coloc.csv'), index=False)
pd.set_option('display.width', 140)
print(res.to_string(index=False))
