# Senescence Commitment as Attractor Bifurcation: A Dynamical Topic Model Framework

**Amruth Deepak Bhat**\
Department of Bioengineering, Indian Institute of Science\
Independent research proposal — computational biology track

------------------------------------------------------------------------

## 1. The Problem

Senescence research has a biomarker problem. p16, p21, SA-β-galactosidase, LMNB1 loss — none work reliably as standalone markers. All appear late, after commitment is irreversible. None predict which cells will senesce versus recover. The field responds by adding more markers, building longer panels, training classifiers on endpoint data.

This approach treats senescence as a state to be identified. We argue it should be treated as a transition to be predicted.

The deeper issue is conceptual. Standard single-cell analysis pipelines take a gene count matrix, perform PCA, embed in UMAP, apply clustering, and annotate cell types. Each step discards biological information: PCA imposes orthogonality that biology does not respect, UMAP destroys global geometry, clustering forces continuous variation into discrete bins, and annotation imports prior assumptions as ground truth. The pipeline discovers the structure it assumed.

The result is a field stuck debating whether senescence has two states or ten, whether light senescence is reversible or not, whether hybrid states are real or artifacts. These debates are unresolvable within the discrete state framework because the framework is wrong.

------------------------------------------------------------------------

## 2. Prior Work and the Gap

A systematic survey of over eighty papers across three fields reveals that the key components of this framework have matured independently but never converged. The gap is not incremental — it is structural.

**Attractor theory for senescence exists only as assumed-parameter ODE models.** Peng et al. (PLOS ONE, 2018) constructed a 13-variable ODE model of core senescence regulators (p53, ATM, Rb, E2F, AKT, p21) and identified three attractors — homeostasis, arrest, and senescence — via Monte Carlo quasi-potential computation. But the model used assumed parameters with no experimental calibration and no single-cell data. The landscape coordinates were pre-defined, not discovered. Galvis et al. (Journal of the Royal Society Interface, 2019) modeled population-level transitions between proliferating, senescent, and apoptotic states, but this is a population dynamics model, not an attractor landscape in transcriptional space. Li et al. (Science, 2020) is the only study combining actual single-cell data with formal attractor landscape computation in any aging context — but it is in yeast, with a 2D toggle-switch, not mammalian senescence. No study has reconstructed a senescence attractor landscape from empirical high-dimensional single-cell transcriptomic data.

Mohit Kumar Jolly's group (IISc) has produced critical slowing down analyses for EMT transitions (Hari et al., PNAS, 2020), and Biplab Bose's group has modeled GRN bistability extensively, but neither has directed these methods toward senescence.

**The bistable switch structure is mathematically established but not connected to data.** Heldt et al. (PNAS, 2018, Tyson/Novak group) modeled the proliferation-quiescence decision as a bistable hysteresis switch in Cdk2/p21 ratio with explicit saddle-node bifurcations. Yao et al. (Nature Cell Biology, 2008) established the Rb-E2F pathway as a bistable switch with saddle-node bifurcation controlling cell cycle commitment. Hsu et al. (Cell, 2019) demonstrated p21 bistability governing the proliferation-senescence decision at single-cell resolution via live imaging. Purvis et al. (Science, 2012) showed that pulsed versus sustained p53 dynamics encode the recovery versus senescence fate. These models all operate in low-dimensional pre-specified circuits and have never been calibrated against or validated with scRNA-seq data. The proposed framework bridges the scales: topic models discover the relevant low-dimensional manifold from scRNA-seq, then dynamical systems analysis runs on that manifold.

**Critical slowing down has been detected for differentiation and EMT — never for senescence.** Mojtahedi et al. (PLoS Biology, 2016) provided the landmark demonstration of CSD in cell fate commitment, observing decreased cell-cell correlation and increased gene-gene correlation approaching the erythroid-myeloid bifurcation, with prolonged recovery time as a direct CSD measurement. This used single-cell qPCR of 17 genes. Richard et al. (PLoS Biology, 2016) independently showed Shannon entropy peaks 8-24 hours before irreversible commitment in chicken erythroid progenitors. For scRNA-seq specifically, BioTIP (Yang et al., Nucleic Acids Research, 2022) extended this to genome-wide data and validated across six datasets including mouse cardiogenesis and gastrulation; FateNet (Sadria and Bury, Bioinformatics, 2024) trains deep learning on simulated SDEs near bifurcations and transfers to pseudotime-ordered scRNA-seq; spliceJAC (Bocci et al., Molecular Systems Biology, 2022) constructs cell-state-specific Jacobian matrices from RNA velocity. None have been applied to senescence. The senescence commitment transition is among the most clinically relevant tipping points in human cell biology and has received zero attention from the early warning signals community.

**Topic models for scRNA-seq are mature but disconnected from dynamical systems.** At least eight distinct approaches exist: cNMF (Kotliar et al., eLife, 2019) using consensus NMF for gene expression programs; scHPF (Levitin et al., Molecular Systems Biology, 2019) using hierarchical Poisson factorization for count data; scETM (Zhao et al., Nature Communications, 2021) using variational autoencoders scaling to $>10^6$ cells; fastTopics/GoM (Carbonetto et al., Genome Biology, 2023) with a principled differential expression framework; Gavish et al. (Nature, 2023) applying per-sample NMF across 1,163 tumors to identify 41 consensus meta-programs. No CTM or dynamic topic model has been adapted for scRNA-seq. No paper uses topic-model-derived coordinates for attractor identification or bifurcation analysis. No paper applies any topic model to senescence.

**The two senescent fate structure is real and unexplained.** The Burnaevskiy et al. finding has been substantially supported. Wechter et al. (Aging, 2023) identified two divergent senescence programs in WI-38 fibroblasts under replicative and DNA damage senescence. Saul et al. (EMBO Journal, 2025) showed that p21$^+$ and p16$^+$ senescent cells are largely non-overlapping populations following independent RNA velocity trajectories in vivo, confirmed by MERFISH spatial transcriptomics. Tao et al. (Cell Metabolism, 2024) identified six distinct Senescence Identities across 602 samples. The mechanism determining which fate a cell takes — whether metabolic state, epigenetic priming, or stochastic fluctuation — is unknown. No paper has applied attractor theory, bifurcation analysis, or variance-based CSD detection to explain the discrete senescent subtypes.

**Transcriptional variance predicts fate transitions broadly but never for senescence.** Chang et al. (Nature, 2008) showed that transcriptome-wide noise in hematopoietic progenitors predicts lineage choice with 7-fold differential fate propensity. Shaffer et al. (Nature, 2017) demonstrated pre-existing transcriptional variability in rare melanoma cells predicts drug resistance. VarID2 (Herman et al., Genome Biology, 2023) quantifies biological noise at single-cell resolution and found increased noise in aged hematopoietic stem cells — the closest existing work to noise quantification in an aging context, but not senescence. No paper has measured variance dynamics before senescence entry.

**Optimal transport and trajectory inference have barely touched senescence.** Waddington-OT (Schiebinger et al., Cell, 2019) has never been applied to senescence data. The few trajectory analyses that exist use basic pseudotime. Chan et al. (eLife, 2022, Calico) produced the most comprehensive replicative senescence multi-omics timecourse but used no formal trajectory inference algorithm, only UMAP visualization.

**The pseudohypoxia-senescence link is real but bidirectional and computationally unexplored.** HIF-1α can both promote and suppress senescence depending on context. Wu and Picone (Aging Cell, 2023) showed HIF-1α drives a non-canonical SASP across multiple senescence triggers, independent of cGAS/STING. Wiley et al. (Cell Metabolism, 2016) defined mitochondrial dysfunction-associated senescence (MiDAS) via NAD$^+$/NADH decline without measuring HIF-1α stabilization, leaving the pseudohypoxia-MiDAS connection as an untested hypothesis. No single-cell analysis has distinguished HIF-1α-high from HIF-1α-low senescent subpopulations. Topic model decomposition could determine whether pseudohypoxic signaling defines the cytoprotective versus tissue-remodeling fate split.

**The gap is complete.** No existing work combines any two of the three core components — topic models for cell state decomposition, dynamical systems/attractor analysis, and critical slowing down detection — for senescence or any cell fate transition. The nearest partial approximations are instructive: BioTIP detects CSD in scRNA-seq but uses standard clustering and has never been applied to senescence; Peng et al. models senescence as an attractor but without data; Li et al. combines single-cell data with attractor computation but in yeast with 2D models; Hsu et al. demonstrates senescence bistability at single-cell resolution but in a one-dimensional circuit. The proposed framework fills a void at the intersection of three mature methodological traditions by doing what none has done: using topic models to discover the relevant low-dimensional manifold from high-dimensional scRNA-seq, fitting dynamical systems analysis on that manifold, and detecting the critical slowing down signature that theory predicts should precede the bifurcation.

------------------------------------------------------------------------

## 3. The Reframing

We propose treating cell states as attractors in a continuous dynamical landscape rather than discrete bins in a categorical space.

Formally, the true state of a cell is a point $\pi_c$ on a low-dimensional manifold $\mathcal{M}$ embedded in gene expression space. This manifold is carved out by gene regulatory network dynamics:

$$\frac{d\pi}{dt} = f(\pi)$$

Stable cell types are attractors: fixed points $\pi^*_k$ where $f(\pi^*_k) = 0$ and perturbations decay back. Near each attractor, the dynamics are approximately linear, described by the Jacobian $J_k = \frac{\partial f}{\partial \pi}\big|_{\pi^*_k}$.

The eigenvectors of $J_k$ define the natural coordinate system near attractor $k$. Eigenvectors with small eigenvalue magnitude (soft modes) correspond to directions of weak regulatory stiffness — programs that are primed but labile. Eigenvectors with large eigenvalue magnitude (stiff modes) correspond to tightly regulated programs that resist perturbation.

**Bifurcation** occurs when a soft mode eigenvalue approaches zero: the attractor loses stability in that direction. Perturbations along the soft mode no longer decay. The cell drifts irreversibly toward a new attractor. This is senescence commitment — not a switch, not a gradual decline, but a loss of dynamical stability in a specific program direction.

The critical observable consequence: **variance in the soft mode direction increases as the cell approaches bifurcation**. The attractor is losing its grip. Cells spread out along the instability direction before committing. This is critical slowing down, a universal signature of systems approaching bifurcation, documented in climate tipping points, ecological collapse, and financial crashes. It has not been systematically applied to single-cell transcriptomics as a pre-senescence biomarker.

------------------------------------------------------------------------

## 4. The Mathematical Framework

We operationalize this through three layers.

### Layer 1: Observation Model

Each cell $c$ produces an observed count vector $x_c \in \mathbb{Z}_{\geq 0}^G$ via multinomial sampling from its true relative expression state $\pi_c \in \Delta^{G-1}$:

$$x_c \sim \text{Multinomial}(N_c, \pi_c)$$

where $N_c$ is total transcript count (library size). The true biological state is $\pi_c$, a point on the probability simplex. The counts are noisy integer observations of it.

### Layer 2: State Space Model

The true states $\pi_c$ lie on a low-dimensional manifold $\mathcal{M} \subset \Delta^{G-1}$ defined by GRN dynamics. Cells occupy attractor neighborhoods and transition paths between them. Near attractor $k$, the covariance structure of cells reflects the Jacobian eigenvector geometry:

$$\Sigma_k \approx \sum_i \frac{\sigma^2}{2|\lambda_i|} v_i v_i^T$$

where $v_i$ are Jacobian eigenvectors and $\lambda_i$ their eigenvalues. **The top eigenvectors of the empirical covariance matrix approximate the soft eigenvectors of the GRN Jacobian.** Directions of maximum transcriptional variance are directions of minimum regulatory stiffness. Statistics and dynamics converge on the same directions.

### Layer 3: Dynamic Correlated Topic Model

The manifold is approximated by a dynamic correlated topic model (DCTM). Topic gene distributions $\phi_{k,t}$ evolve as a Gaussian random walk over time, capturing attractor drift under sustained damage:

$$\phi_{k,t} | \phi_{k,t-1} \sim \mathcal{N}(\phi_{k,t-1}, \sigma^2_\phi I)$$

Cell topic weights follow a logistic-normal distribution with full covariance, allowing positive correlations between topics (unlike standard LDA which uses a Dirichlet prior and enforces effective anticorrelation through the simplex constraint):

$$\eta_{c,t} \sim \mathcal{N}(\mu_t, \Sigma_t), \quad \theta_{c,t} = \text{softmax}(\eta_{c,t})$$

The observation model connects to the data:

$$x_{c,t} \sim \text{Multinomial}\left(N_{c,t}, \sum_k \theta_{ck,t} \phi_{k,t}\right)$$

The covariance matrix $\Sigma_t$ evolves over time, encoding changes in program-program coupling as the damage response progresses. Its eigenvalue spectrum is the primary readout.

Population-level flow between timepoints is estimated via **optimal transport** rather than RNA velocity. The transport map between consecutive timepoints minimizes displacement cost in topic space under the Fisher information metric, avoiding the unstable splicing ratio assumptions that make RNA velocity unreliable. This gives a flow field on topic space grounded in the actual population dynamics.

------------------------------------------------------------------------

## 5. The Hypothesis

**HIF1α target genes define the soft eigenvector direction along which the normal cell attractor loses stability under oxidative stress.**

The mechanistic argument: HIF1α is basally expressed but constitutively degraded by PHD enzymes under normoxia via oxygen-dependent prolyl hydroxylation. This degradation is metabolic and post-translational — a soft restoring force, not an epigenetic lock. Oxidative stress perturbs mitochondrial function, destabilizes TCA cycle intermediates including succinate and fumarate that drive PHD activity, and reduces PHD-mediated HIF1α degradation. The pseudohypoxic program, always latent at low probability, becomes progressively less regulated. The restoring force weakens. The eigenvalue magnitude of this mode decreases toward zero.

If this is correct: - The soft eigenvector of $\Sigma_t$ restricted to normal-state cells should show anomalous eigenvalue growth under damage - Genes with highest loadings on this eigenvector should be enriched for HIF1α targets (VEGF, LDHA, PDK1, BNIP3, GLUT1) - Eigenvalue growth should precede classical senescence marker expression - The cytoprotective versus tissue-remodeling fate split identified by Burnaevskiy et al. should be predictable from topic space position at bifurcation

------------------------------------------------------------------------

## 6. The Dataset

**GSE223128**: Burnaevskiy, Oshima & Mendenhall (GeroScience, 2023). Large-scale scRNA-seq timecourse of human fibroblasts under oxidative stress (55 μM H2O2). Timepoints: untreated, 4h, 1 day, 2 days, 3 days, 4 days, 7 days post-stress. Additional conditions: nutlin-3a-induced senescence-like arrest, quiescence (0.2% FBS). The original study identified two transcriptionally distinct senescent fates (cytoprotective CP and tissue-remodeling TR) and noted metabolic divergence at day 1 as the earliest distinguishing feature.

This dataset is ideal because: - 7 timepoints provide genuine temporal resolution for dynamic modeling - The metabolic divergence at day 1 is consistent with the pseudohypoxia hypothesis - The CP/TR fate split provides a binary outcome for bifurcation proximity prediction - Quiescence condition provides a negative control: quiescent cells should not show bifurcation dynamics - Nutlin condition provides a mechanistically distinct senescence path for comparison

------------------------------------------------------------------------

## 7. Analysis Plan

### Step 1: Preprocessing and quality control

Standard scRNA-seq pipeline. Library size normalization, log1p transformation, highly variable gene selection, batch correction across timepoints using Harmony. Reproduce the original UMAP and cluster structure to verify data integrity.

### Step 2: Topic model fitting

Fit a correlated topic model (CTM) with $K$ topics selected by held-out perplexity and biological coherence. Initialize with $K \in \{5, 8, 12, 15\}$ and compare. Use the Python package `scholar` or implement in Pyro. Validate topic coherence by gene set enrichment of top genes per topic against MSigDB Hallmarks.

The CTM produces per-cell topic weight vectors $\theta_{c,t}$ and topic gene distributions $\phi_k$.

### Step 3: Dynamic extension

Fit the dynamic CTM jointly across timepoints, allowing topic gene distributions to drift. Compare static versus dynamic fit by held-out likelihood. Examine how topic gene distributions shift over the 7-day timecourse.

### Step 4: Subpopulation identification and bifurcation proximity measurement

A critical complication: after day 1, the cell population may have already split into multiple coexisting subpopulations, each occupying distinct regions of topic space and evolving along different trajectories. Computing a single covariance matrix over all cells near the normal topic conflates cells from different dynamical regimes, producing eigenvalue growth that reflects population heterogeneity rather than attractor instability. This must be addressed before bifurcation proximity can be measured.

**Step 4a: Soft subpopulation decomposition.** At each timepoint, fit a Gaussian mixture model (GMM) in topic space with $M$ components, where $M$ is selected by BIC. Each mixture component defines a soft subpopulation: cell $c$ belongs to component $m$ with probability $r_{cm}$, the posterior responsibility. This preserves soft assignment and avoids hard clustering.

The number of components $M$ may increase across timepoints as the population diversifies. Tracking $M_t$ over time is itself a signal: attractor proliferation, the emergence of new stable subpopulations, manifests as an increase in the number of mixture components needed to fit the data.

Subpopulation correspondence across timepoints is established using optimal transport: the mixture component at $t+1$ that receives the most probability mass from component $m$ at $t$ is its successor. This gives a lineage graph of subpopulations over time without tracking individual cells.

**Step 4b: Per-subpopulation covariance and eigendecomposition.** For each subpopulation $m$ at each timepoint, compute the responsibility-weighted covariance matrix:

$$\hat{\Sigma}_{m,t} = \frac{\sum_c r_{cm} (\theta_{c,t} - \bar\theta_{m,t})(\theta_{c,t} - \bar\theta_{m,t})^T}{\sum_c r_{cm}}$$

Eigendecompose $\hat{\Sigma}_{m,t}$ separately for each subpopulation. Track $\lambda_{1,m,t}$, the largest eigenvalue of subpopulation $m$, over time along its lineage.

The bifurcation signal is expected in the subpopulation lineage that traces back to the normal proliferating state at day 0. Subpopulations that emerge post-bifurcation (CP and TR attractors) should show stable or decreasing $\lambda_{1,m,t}$ as they consolidate into new attractors.

**Step 4c: Critical slowing down fit.** Fit the critical slowing down curve to the $\lambda_{1,m,t}$ time series of the normal-lineage subpopulation:

$$\lambda_{1,m,t} \approx \frac{A}{t^* - t} + \lambda_0$$

Estimate $t^*$, the predicted bifurcation time. The bifurcation time should precede the timepoint at which the normal-lineage subpopulation splits into identifiable CP and TR successor subpopulations in the GMM lineage graph.

**Step 4d: Subpopulation switching.** Cells that switch subpopulation membership between timepoints (high $r_{cm}$ at $t$, high $r_{cm'}$ at $t+1$ for $m \neq m'$) are identified from the optimal transport map. High switching rates indicate cells in transition regions between attractors. These cells should have intermediate bifurcation proximity scores and ambiguous fate assignment, which is a testable prediction.

### Step 5: Soft eigenvector validation

Identify the eigenvector $v_{1,m,t}$ associated with $\lambda_{1,m,t}$ for the normal-lineage subpopulation. Map it back to gene space: for each gene $g$, compute its loading as the dot product of the gene's topic weight vector with $v_{1,m,t}$.

Run gene set enrichment analysis on genes ranked by loading magnitude. Primary test: enrichment for MSigDB Hallmark Hypoxia gene set. Secondary: HIF1α ChIP-seq target genes from ENCODE.

**Confirmation criterion**: top-loaded genes are significantly enriched for HIF1α targets (FDR \< 0.05, NES \> 1.5).

**Refutation criterion**: top-loaded genes are enriched for a different program (e.g. p53 response, inflammatory signaling) with no HIF1α enrichment.

### Step 6: Fate prediction

Compute per-cell bifurcation proximity scores at each timepoint using the normal-lineage soft eigenvector:

$$\text{ProxBif}_{c,t} = (\theta_{c,t} - \bar\theta_{m,t})^T v_{1,m,t}$$

where $m$ is the normal-lineage subpopulation and cell weights are responsibility-weighted. Cells assigned to post-bifurcation subpopulations (CP or TR) receive scores based on their topic space distance from the bifurcation point rather than from the current subpopulation mean.

Test whether $\text{ProxBif}_{c,t}$ at day 1 and day 2 predicts final assignment to CP versus TR fate at day 7. Additionally test whether cells with high subpopulation switching rates at day 1-2 show intermediate ProxBif scores, consistent with being in the transition region at bifurcation. Compare predictive performance (AUROC) against standard senescence scoring methods (ssGSEA with SenMayo, p21 expression, p16 expression).

**Confirmation criterion**: ProxBif at day 1-2 predicts day 7 fate with AUROC \> 0.70, outperforming classical markers. High-switching cells cluster near the bifurcation point in topic space.

**Refutation criterion**: ProxBif performs no better than classical markers. Switching cells are distributed randomly rather than near the bifurcation point.

### Step 7: Optimal transport flow

Compute Sinkhorn optimal transport maps between consecutive timepoints in topic space using the Python Optimal Transport (POT) library. Visualize flow vectors on the topic space manifold. Confirm that flow trajectories diverge toward CP versus TR attractors in a manner consistent with bifurcation proximity scores.

### Step 8: Negative controls

Apply the same pipeline to: - Quiescent cells: should show stable $\lambda_{1,t}$ with no critical slowing down signature - Nutlin condition: should show bifurcation dynamics but along a different eigenvector direction (p53 axis rather than pseudohypoxia axis)

These controls are essential. If the critical slowing down signature appears in quiescence, the result is artifactual.

------------------------------------------------------------------------

## 8. What Confirmation and Refutation Look Like

| Component | Confirmation | Refutation |
|------------------------|------------------------|------------------------|
| Eigenvalue growth | Anomalous $\lambda_{1,m,t}$ growth in normal-lineage subpopulation preceding day 7 markers | Flat or noisy spectrum; growth in all subpopulations equally |
| Subpopulation proliferation | $M_t$ increases after bifurcation, not before | $M_t$ constant or increases before eigenvalue growth |
| Soft eigenvector identity | HIF1α target enrichment on $v_{1,m,t}$ | Different pathway enriched |
| Fate prediction | ProxBif AUROC \> 0.70 at day 1-2; switching cells cluster near bifurcation point | ProxBif no better than p21; switching cells randomly distributed |
| Negative controls | No slowing down in quiescence; different eigenvector in nutlin | Slowing down in quiescence (artifact) |
| Dynamic CTM | Topics drift coherently; subpopulation lineage graph is interpretable | Topics unstable; lineage graph disconnected |

Refutation of the pseudohypoxia hypothesis does not refute the framework. The framework survives if critical slowing down exists along any biologically coherent direction. The pseudohypoxia claim is an additional bet on top of the framework.

------------------------------------------------------------------------

## 9. Tools and Implementation

| Task                           | Tool                             |
|--------------------------------|----------------------------------|
| Preprocessing                  | Scanpy (Python)                  |
| Batch correction               | Harmony                          |
| Topic model (CTM)              | Pyro or scholar                  |
| Subpopulation decomposition    | Gaussian mixture model (sklearn) |
| Subpopulation lineage tracking | Python Optimal Transport (POT)   |
| Optimal transport flow         | Python Optimal Transport (POT)   |
| Gene set enrichment            | GSEApy                           |
| Visualization                  | matplotlib, seaborn, scanpy      |

All code will be deposited on GitHub with reproducible Snakemake pipeline. Data is publicly available on GEO (GSE223128).

------------------------------------------------------------------------

## 10. Significance

If the framework is confirmed, it offers three things the field currently lacks.

**A pre-senescence biomarker**: bifurcation proximity is detectable before irreversible commitment, enabling intervention at the relevant window.

**A resolution to the reversibility debate**: reversibility is a question of attractor geometry, not marker expression. Cells before the bifurcation retain the normal attractor and can recover. Cells past it cannot, not because of epigenetic locks per se, but because the dynamical object they would relax toward no longer exists.

**A principled mathematical foundation**: the framework connects single-cell transcriptomics to dynamical systems theory in a way that is falsifiable, generalizable beyond senescence, and grounded in first principles rather than heuristic pipeline choices.

------------------------------------------------------------------------

## 11. Limitations and Honest Uncertainties

The topic model approximation is linear where the true manifold is curved. Transition cells may be poorly represented. The number of topics $K$ requires principled selection and the result is sensitive to this choice. The pseudohypoxia hypothesis is biologically motivated but unconfirmed — it could be wrong while the framework is right. GSE223128 uses fibroblasts under oxidative stress; generalizability to other cell types and senescence inducers requires validation on additional datasets. Critical slowing down theory was developed for simple dynamical systems — its application to high-dimensional transcriptional systems involves assumptions that are not fully validated.

Subpopulation tracking introduces additional uncertainty. Gaussian mixture model components are not guaranteed to correspond to biologically coherent attractors. Component number selection by BIC may underfit early timepoints where subpopulation structure is subtle, or overfit late timepoints where structure is clear. The optimal transport lineage graph assumes smooth population flow between timepoints; if cells undergo rapid state transitions between sampling points, the correspondence may be unreliable. The minimum number of cells per subpopulation per timepoint needed for reliable covariance estimation is approximately $10K$ where $K$ is the number of topics — with sparse subpopulations at early timepoints, covariance estimates will be noisy and eigenvalue dynamics unreliable. This should be assessed by bootstrap resampling of the covariance estimates at each timepoint.

These are not reasons to abandon the approach. They are specific, addressable limitations that define the scope of the initial paper and the roadmap for follow-up work.

------------------------------------------------------------------------

*Proposal drafted: March 2026*\
*Status: Pre-analysis. All claims are hypotheses pending computational validation on GSE223128.*

------------------------------------------------------------------------