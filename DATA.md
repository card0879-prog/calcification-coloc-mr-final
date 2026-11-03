# Data provenance and terms of use

This repository redistributes **small per-locus subsets** of publicly available
summary statistics and LD, solely to make the analysis reproducible. These data
are **third-party data**: they remain under the terms of their original providers,
and the authors of this repository claim **no ownership and grant no licence** over
them. Only the analysis **code** is licensed by us (MIT; see `LICENSE`).

Anyone reusing the data subsets must cite the original sources and comply with each
provider's data-use policy. The `data/` files are regional slices; full datasets are
available from the sources below.

> **Two sources are NOT redistributed here** — the UK-Biobank–based Sinnott-Armstrong
> biomarker GWAS and FinnGen R12. Both are publicly downloadable free of charge; download
> them yourself and place the six regional files under `data/` following
> [`data/RESTRICTED_INPUTS.md`](data/RESTRICTED_INPUTS.md). The scripts then run unchanged.
> Everything else needed to run the pipeline **is** included.

| Source (files) | Provider & citation | Access | Reuse terms |
|---|---|---|---|
| Activity / biomarker-proxy exposures (`LPA_Lpa_b37.tsv`, `ENZ_*_activity_b37.tsv`) — **NOT included, download yourself** | Sinnott-Armstrong et al., meta-analysed **UK Biobank** biomarker GWAS, *Nat Genet* 2021 | NIH figshare 12355382 | Derived from UK Biobank; not redistributed here — see `data/RESTRICTED_INPUTS.md` |
| Blood cis-eQTL (`*_eQTLGen_b37.tsv`, `axis_eQTLGenBlood_b37.tsv`, `ankh_af_region.tsv`) | eQTLGen Consortium, *Nat Genet* 2021 | eqtlgen.org | Freely available for research; cite eQTLGen |
| Tissue cis-eQTL (`ENZ_*_{liver,blood}GTEx_b38.tsv`, `axis_GTEx*_b38.tsv`) | eQTL Catalogue / **GTEx v8** | ebi.ac.uk/eqtl | Open licence (eQTL Catalogue); cite eQTL Catalogue and GTEx |
| CAVS outcome (`*CAVS_b38.tsv`) — **NOT included, download yourself** | **FinnGen** release R12, `I9_CAVS_OPERATED` | finngen.fi/en/access_results | Public under FinnGen's data policy — acknowledge FinnGen; not redistributed here — see `data/RESTRICTED_INPUTS.md` |
| CAC outcome (`*CAC_b38.tsv`) | Kavousi et al. 2023, GWAS Catalog **GCST90278456**, quantitative ln(Agatston) | ebi.ac.uk/gwas | GWAS Catalog harmonised summary statistics; cite the original study |
| LD reference (`LD_*.ld.gz`, `LD_*.bim`) | **1000 Genomes** phase-3 EUR, 503 individuals, GRCh37, PLINK `--r square` | internationalgenome.org | Open |

To avoid any question over third-party redistribution rights, the UK-Biobank–derived
Sinnott-Armstrong biomarker GWAS and the FinnGen R12 files are **not shipped** in this
repository; download them yourself per `data/RESTRICTED_INPUTS.md` (both are free). Every
other input needed to run the pipeline is included and openly redistributable.

Code: MIT (`LICENSE`). Data: see the table above; no blanket licence is asserted over
third-party summary statistics.
