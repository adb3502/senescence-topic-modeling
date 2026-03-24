# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a research writing repository, not a software project. It contains three long-form scientific documents exploring the intersection of dynamical systems theory, single-cell transcriptomics, and computational biology.

## Documents

**CSD.md** — A review synthesizing literature on critical slowing down (CSD) as a universal precursor to biological state transitions. Covers mathematical foundations (bifurcation theory, quasi-potential landscapes, Ornstein-Uhlenbeck processes), empirical validations in hematopoiesis and EMT, the computational toolkit (BioTIP, scDiffEq, MuTrans, dynamo, CellRank, Waddington-OT), the aging-cancer connection via dynamical instability, and the limitations catalogued by Boettiger and Hastings (prosecutor's fallacy, noise-induced transitions, statistical power).

**NMF.md** — A review comparing unsupervised matrix factorization methods (cNMF, LDA/fastTopics, ICA, MOFA, archetypal analysis, deep generative models) across biological data modalities, with emphasis on population-scale plasma proteomics. Critiques the supervised organ-clock paradigm (Oh et al., Nature 2023) and maps five research gaps where NMF/topic models applied to Olink/SomaScan/DIA-MS data could reshape the field.

**senescence_attractor_proposal.md** — A pre-analysis research proposal (Amruth Deepak Bhat, IISc Bioengineering, March 2026) for treating senescence commitment as an attractor bifurcation detectable via critical slowing down in topic-model coordinates. Uses GSE223128 (Burnaevskiy et al., GeroScience 2023: 7-timepoint scRNA-seq of human fibroblasts under H2O2). The framework layers a multinomial observation model, a Gaussian state-space model, and a dynamic correlated topic model (DCTM), using Gaussian mixture models for subpopulation tracking and optimal transport for flow estimation.

## Key Concepts Across All Three Documents

- **Critical slowing down**: rising variance and autocorrelation as a system's dominant eigenvalue approaches zero near bifurcation
- **Dynamic Network Biomarker (DNB)**: gene groups with elevated intra-correlation and elevated SD near tipping points (Chen, Liu, Aihara 2012)
- **cNMF**: consensus NMF (Kotliar et al., eLife 2019) — gold standard for gene program discovery in scRNA-seq
- **MOFA/MOFA+**: multi-omics factor analysis with ARD priors, suited for continuous/negative data like Olink NPX
- **Organ clocks**: supervised plasma proteomic aging models using GTEx tissue labels (Oh et al., Nature 2023)
- **GSE223128**: the primary dataset for the proposal — publicly available on GEO
- **Proposed implementation stack**: Scanpy, Harmony, Pyro or scholar (CTM), sklearn (GMM), Python Optimal Transport (POT), GSEApy

## Writing Style and Conventions

These documents are written in a direct, literature-dense academic register. Claims are anchored to specific papers with author names, journal, and year inline. Refutation criteria are stated explicitly alongside confirmation criteria. The proposal distinguishes framework claims from the specific pseudohypoxia hypothesis (HIF1α soft eigenvector direction) as a separable bet.
