# Data provenance

Small per-locus subsets of public summary statistics and LD, re-distributed for reproducibility.

- Exposures (activity/biomarker proxy): Sinnott-Armstrong et al., meta-analyzed UK Biobank biomarker GWAS, Nat Genet 2021 (`LPA_Lpa_b37.tsv`, `ENZ_*_activity_b37.tsv`).
- Blood cis-eQTL: eQTLGen Consortium 2019 (`*_eQTLGen_b37.tsv`, `axis_eQTLGenBlood_b37.tsv`); allele frequencies `ankh_af_region.tsv`.
- Tissue cis-eQTL: eQTL Catalogue / GTEx v8 — aorta/coronary/whole blood/liver (`ENZ_*_{liver,blood}GTEx_b38.tsv`, `axis_GTEx*_b38.tsv`).
- CAVS outcome: FinnGen R12, I9_CAVS_OPERATED (`*CAVS_b38.tsv`).
- CAC outcome: Kavousi et al. 2023, GWAS Catalog GCST90278456, quantitative ln(Agatston) (`*CAC_b38.tsv`).
- LD: 1000 Genomes phase-3 EUR, 503 individuals, GRCh37, PLINK `--r square` (`LD_*.ld.gz`, `LD_*.bim`).

Region files shared under CC-BY 4.0 with attribution; code under MIT (see LICENSE).
