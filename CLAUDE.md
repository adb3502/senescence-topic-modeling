# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains three long-form scientific documents and an active R analysis pipeline for single-cell RNA-seq topic modeling. The scientific goal is to detect senescence commitment as an attractor bifurcation via critical slowing down in topic-model coordinates, using GSE223128 (IMR-90 fibroblasts, H2O2 senescence timecourse).

## Documents

**CSD.md** — A review synthesizing literature on critical slowing down (CSD) as a universal precursor to biological state transitions. Covers mathematical foundations (bifurcation theory, quasi-potential landscapes, Ornstein-Uhlenbeck processes), empirical validations in hematopoiesis and EMT, the computational toolkit (BioTIP, scDiffEq, MuTrans, dynamo, CellRank, Waddington-OT), the aging-cancer connection via dynamical instability, and the limitations catalogued by Boettiger and Hastings (prosecutor's fallacy, noise-induced transitions, statistical power).

**NMF.md** — A review comparing unsupervised matrix factorization methods (cNMF, LDA/fastTopics, ICA, MOFA, archetypal analysis, deep generative models) across biological data modalities, with emphasis on population-scale plasma proteomics. Critiques the supervised organ-clock paradigm (Oh et al., Nature 2023) and maps five research gaps where NMF/topic models applied to Olink/SomaScan/DIA-MS data could reshape the field.

**senescence_attractor_proposal.md** — A pre-analysis research proposal (Amruth Deepak Bhat, IISc Bioengineering, March 2026) for treating senescence commitment as an attractor bifurcation detectable via critical slowing down in topic-model coordinates. Uses GSE223128 (Burnaevskiy et al., GeroScience 2023: 7-timepoint scRNA-seq of human fibroblasts under H2O2). The framework layers a multinomial observation model, a Gaussian state-space model, and a dynamic correlated topic model (DCTM), using Gaussian mixture models for subpopulation tracking and optimal transport for flow estimation.

## Dataset: GSE223128

- **5 SPLIT-seq samples**: sc2, sc78 (controls — Nutlin/OxStress/Quiescence conditions), sc5, sc6, sc9 (timecourse — Control → 4hrs → D1 → D2 → D3 → D4 → D7)
- **61,285 author-filtered cells** total
- **SPLIT-seq barcodes**: composite `bc1_bc2_well_bc`; well barcode = last `_`-separated component; maps to condition via per-sample `bc_samples.csv`
- **Raw data**: not tracked in git (4.4GB). Re-download from GEO or copy from source machine

## Analysis Pipeline (R)

All analysis is in R 4.5.3. Shared library at `C:/Program Files/R/4.1` (intentional — shared across R versions, do not treat as a path error).

### Scripts

| Script | Purpose |
|---|---|
| `data/fetch_geo_metadata.R` | Fetch GSE223128 SOFT metadata via GEOquery |
| `data/load_matrices.R` | Load raw MTX, map SPLIT-seq well barcodes to conditions |
| `analysis/01_qc.R` | QC pipeline: rhdf5 extracts author-filtered barcodes from h5ad, Seurat v5 objects, pct_mt/pct_ribo |
| `analysis/01_qc_plots.R` | Publication-quality QC figures (cowplot, Nature palette, SVG+PNG) |
| `analysis/03_lda_poc.R` | Poisson NMF / LDA POC on sc2 and sc78; K-selection with checkpointing |

### Key Implementation Decisions

- **Raw counts for LDA**: do NOT log-normalize. Multinomial generative model assumes raw integer counts; library size handled implicitly
- **fastTopics**: Poisson NMF (KL divergence, SCD + extrapolation). `fit_topic_model()` with `verbose="none"`. Log-likelihood via `mean(loglik_multinom_topic_model(X, fit))`
- **Corpus-style gene filtering**: keep genes detected in 1–95% of cells, remove MT- and RP[SL] genes. Results in ~10,725 genes from 57,905
- **h5ad loading**: SeuratDisk incompatible with Seurat v5. Use `rhdf5` to extract `/obs/_index` barcodes, then filter raw MTX matrices
- **Checkpointing**: each fitted model saved as `analysis/lda_poc/{prefix}_fit_k{k}.rds` before proceeding to next K
- **GPU acceleration**: target workstation has RTX 5060 (WSL2). Use `torch` R package with LibTorch CUDA backend for Poisson NMF

### Gitignore

Large files excluded from git:
- `data/raw/` — GEO downloads (h5ad, MTX, tar)
- `analysis/qc_output/*.rds` — Seurat objects (~288MB)
- `analysis/qc_output/qc_metrics.csv` — cell metrics (~166MB)
- `analysis/lda_poc/*.rds` — model checkpoints

## Key Concepts

- **Critical slowing down**: rising variance and autocorrelation as a system's dominant eigenvalue approaches zero near bifurcation
- **Dynamic Network Biomarker (DNB)**: gene groups with elevated intra-correlation and elevated SD near tipping points (Chen, Liu, Aihara 2012)
- **cNMF**: consensus NMF (Kotliar et al., eLife 2019) — gold standard for gene program discovery in scRNA-seq
- **MOFA/MOFA+**: multi-omics factor analysis with ARD priors, suited for continuous/negative data like Olink NPX
- **Organ clocks**: supervised plasma proteomic aging models using GTEx tissue labels (Oh et al., Nature 2023)
- **Carbonetto et al. proof**: Poisson NMF = multinomial LDA at MLE — justifies using fastTopics for topic modeling

## Writing Style and Conventions

Documents are written in a direct, literature-dense academic register. Claims are anchored to specific papers with author names, journal, and year inline. Refutation criteria are stated explicitly alongside confirmation criteria. The proposal distinguishes framework claims from the specific pseudohypoxia hypothesis (HIF1α soft eigenvector direction) as a separable bet.
