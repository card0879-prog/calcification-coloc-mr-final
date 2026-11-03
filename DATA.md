# Data provenance and terms of use

This repository redistributes **small per-locus subsets** of publicly available
summary statistics and LD, solely to make the analysis reproducible. These data
are **third-party data**: they remain under the terms of their original providers,
and the authors of this repository claim **no ownership and grant no licence** over
them. Only the analysis **code** is licensed by us (MIT; see `LICENSE`).

Anyone reusing the data subsets must cite the original sources and comply with each
provider's data-use policy. The `data/` files are regional slices; full datasets are
available from the sources below.

| Source (files) | Provider & citation | Access | Reuse terms |
|---|---|---|---|
| Activity / biomarker-proxy exposures (`LPA_Lpa_b37.tsv`, `ENZ_*_activity_b37.tsv`) | Sinnott-Armstrong et al., meta-analysed **UK Biobank** biomarker GWAS, *Nat Genet* 2021 | Authors' public release | **Derived from UK Biobank** — reuse under UK Biobank's access terms; confirm redistribution is permitted before further sharing |
| Blood cis-eQTL (`*_eQTLGen_b37.tsv`, `axis_eQTLGenBlood_b37.tsv`, `ankh_af_region.tsv`) | eQTLGen Consortium, *Nat Genet* 2021 | eqtlgen.org | Freely available for research; cite eQTLGen |
| Tissue cis-eQTL (`ENZ_*_{liver,blood}GTEx_b38.tsv`, `axis_GTEx*_b38.tsv`) | eQTL Catalogue / **GTEx v8** | ebi.ac.uk/eqtl | Open licence (eQTL Catalogue); cite eQTL Catalogue and GTEx |
| CAVS outcome (`*CAVS_b38.tsv`) | **FinnGen** release R12, `I9_CAVS_OPERATED` | finngen.fi/en/access_results | Public under FinnGen's data policy — acknowledge FinnGen; confirm redistribution terms before further sharing |
| CAC outcome (`*CAC_b38.tsv`) | Kavousi et al. 2023, GWAS Catalog **GCST90278456**, quantitative ln(Agatston) | ebi.ac.uk/gwas | GWAS Catalog harmonised summary statistics; cite the original study |
| LD reference (`LD_*.ld.gz`, `LD_*.bim`) | **1000 Genomes** phase-3 EUR, 503 individuals, GRCh37, PLINK `--r square` | internationalgenome.org | Open |

**Two sources warrant a direct permissions check before this archive is treated as
freely redistributable**: the UK-Biobank–derived Sinnott-Armstrong biomarker GWAS and
FinnGen R12. Both are publicly downloadable, but their data-use policies govern
redistribution. If either does not permit redistributing summary subsets, replace the
corresponding `data/` files with the access links above and run the pipeline against
locally downloaded inputs — the scripts read the same file names from `data/`.

Code: MIT (`LICENSE`). Data: see the table above; no blanket licence is asserted over
third-party summary statistics.
