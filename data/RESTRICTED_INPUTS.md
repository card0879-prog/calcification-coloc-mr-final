# Inputs you must download yourself (not redistributed here)

Six regional input files are **not included** in this repository, because they derive
from two sources whose data-use policies govern redistribution: the **UK-Biobank**–based
Sinnott-Armstrong biomarker GWAS and **FinnGen R12**. Both are **publicly downloadable
free of charge** from the links below. Download the full summary statistics, extract the
listed genomic window, and save each as the exact file name under `data/`. The analysis
scripts then run unchanged.

## From Sinnott-Armstrong et al. (meta-analysed UK Biobank biomarker GWAS, *Nat Genet* 2021)

Source: NIH figshare, dataset **12355382**
(https://nih.figshare.com/articles/dataset/12355382). Build **GRCh37**.
Expected columns (tab-separated): `CHROM  POS  MarkerName  REF  ALT  Effect  StdErr  P-value  MAF`
(MarkerName format `CHR:POS_REF_ALT`; Effect/StdErr are per-SD on the ALT allele).

| Save as `data/…` | Biomarker | Region to extract (GRCh37) |
|---|---|---|
| `LPA_Lpa_b37.tsv` | Lipoprotein(a) | chr6:160,400,000–161,600,000 |
| `ENZ_GGT1_activity_b37.tsv` | Gamma-glutamyltransferase | chr22:24,800,000–25,250,000 |
| `ENZ_UGT1A1_activity_b37.tsv` | Total bilirubin | chr2:234,500,000–234,850,000 |

## From FinnGen release R12 — endpoint `I9_CAVS_OPERATED` (CAVS outcome)

Source: FinnGen public results (https://finngen.fi/en/access_results). Build **GRCh38**.
Expected columns (tab-separated): `chrom  pos  ref  alt  rsids  pval  beta  sebeta  af_alt`
(native FinnGen summary-statistic columns).

| Save as `data/…` | Region to extract (GRCh38) |
|---|---|
| `LPA_CAVS_b38.tsv` | chr6:160,000,000–161,200,000 |
| `alpl_CAVS_b38.tsv` | chr1:21,000,000–22,100,000 |
| `axis_CAVS_b38.tsv` | chr5:13,700,000–14,000,000; chr6:131,600,000–132,950,000; chr17:47,500,000–47,700,000 (ANKH, ENPP1, PHOSPHO1 windows; concatenate into one file, keeping the header once) |

Please cite the original studies and acknowledge FinnGen and UK Biobank as those sources require.
All other `data/` files (eQTLGen, GTEx/eQTL Catalogue, GWAS Catalog CAC, 1000 Genomes LD)
are included and openly redistributable.
