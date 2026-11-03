# PPi–TNAP drug-target MR validity audit — analysis, LD-aware colocalisation & figures

Reproducible code and data for *"Colocalization-gated validity audit of drug-target
Mendelian randomization in the PPi–TNAP axis"* (target journal: **G3: Genes|Genomes|Genetics**).
Every reported statistic is **computed from the deposited summary data and 1000 Genomes
European LD** by `run_reanalysis.py` (MR, conditional/LD-aware analyses) and `run_coloc.py`
(single-causal colocalisation posteriors), which write `results_reanalysis.csv` and
`results_coloc.csv`. The figure scripts **read those computed result files** and do not
embed statistics of their own — no headline value is hard-coded. Note that at the
multi-signal *LPA* locus the single-causal `coloc.abf` posterior is **not** the verdict
basis (it is reported in `results_coloc.csv` only for transparency, flagged in the `basis`
column); the positive control is admitted by the conditional analysis in `run_reanalysis.py`.

## Reproduce

> First obtain six inputs that are **not redistributed here** (UK-Biobank–based
> Sinnott-Armstrong biomarker GWAS and FinnGen R12 CAVS) — both are free; follow
> [`data/RESTRICTED_INPUTS.md`](data/RESTRICTED_INPUTS.md). Every other input is included.

```bash
pip install numpy pandas scipy matplotlib
python run_reanalysis.py     # headline MR + colocalisation statistics -> results_reanalysis.csv
python run_coloc.py          # single-signal coloc.abf cross-check
python fig1.py ; python fig2.py ; python fig3.py ; python fig4.py   # Figures 1–4
```

`run_reanalysis.py` prints and writes `results_reanalysis.csv`:

| pair | recomputed | verdict |
|---|---|---|
| Lp(a) → CAVS | OR 1.29 (1.24–1.34), effect allele harmonised; primary CAVS signal shares the Lp(a) signal (lead-pair r²=0.97; CAVS lead \|Z\| 12.4→1.2 on conditioning) | PASS (shared causal variant) |
| Lp(a) → CAC | β +0.13 (P 3×10⁻⁴); CAC has 0 genome-wide signals at LPA | non-colocalizing nominal |
| GGT activity vs GGT1 blood (eQTLGen) | lead-pair r²=0.99, cross-variant r=−0.89 | reversed-direction (blood misleading) |
| GGT activity vs GGT1 liver | liver eQTL lead |Z|≈4 (GTEx n=208) | under-powered |
| ANKH blood → CAVS | coloc.abf PP.H4=0.27, PP.H1=0.71 | non-colocalizing nominal |
| ENPP1 aorta → CAC | coloc.abf PP.H4=0.04 (quantitative) | fail-coloc |
| ENPP1 blood → CAC | Wald β +0.84 (SE 0.39, P 0.033); eQTLGen cis-lead rs9385589 | nominal only (not colocalizing) |

## Verify integrity (GitHub ↔ Zenodo)

Every file is listed with its SHA-256 in `SHA256SUMS.txt`. The Zenodo archive
(all-versions DOI **10.5281/zenodo.23896548**) is built from this same tree, so a
download can be checked against the repository file-by-file:

```bash
sha256sum -c SHA256SUMS.txt     # run inside the unpacked archive or a clone
```

All drug-target Wald ratios harmonise the outcome effect allele to the exposure
effect allele before taking the β-ratio (see `run_reanalysis.py`), so an estimate
is never produced from a mis-oriented allele.

## Protocol (frozen)

- **Builds/scale.** Exposures GRCh37 (per-SD); eQTL GRCh37 (eQTLGen) / GRCh38 (GTEx). CAVS binary → OR; CAC quantitative ln(Agatston+1) → β (never exponentiated). Cross-build harmonisation by the modal allele-matched positional offset.
- **Instrument.** Exposure-only lead cis-variant by max |Z| in a fixed ±100-kb gene window (not min-P, which underflows at LPA); deterministic tie-break by position.
- **Colocalisation.** `coloc.abf` (single causal variant) at single-signal loci. At multi-signal loci (LPA, UGT1A) the single-causal assumption fails; the number of conditionally independent signals is found by stepwise conditional analysis (|z_cond|>5.45) on the 1000 Genomes EUR LD matrix, and shared causality judged by lead-pair r² and reciprocal conditioning. A hand-rolled SuSiE-RSS was tried and rejected for failing internal validation against the conditional analysis; it is not used.
- **LD.** Per-locus 1000 Genomes phase-3 EUR, 503 individuals, GRCh37 (PLINK `--r square`); matrices in `data/LD_*.ld.gz` with `data/LD_*.bim`.

## Contents
`run_reanalysis.py`, `run_coloc.py`, `coloc_abf.py`, `mr_engine.py`, `fig1–4.py`, `figstyle.py`,
`data/` (region summary statistics, `LD_*` matrices, eQTLGen allele frequencies), `DATA.md`, `LICENSE`, `CITATION.cff`.
