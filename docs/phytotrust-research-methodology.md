# PhytoTRUST: Study-Aware Evidence Transferability Learning

## Research-gap analysis and prespecified methodological design

**Scope:** scientific design only; no algorithm has been implemented or evaluated.  
**Literature cut-off:** 21 October 2026.  
**Intended use:** predict the distribution of a future, independent study's comparable bioactivity measurement—not drug efficacy, safety, or clinical success.

## Executive verdict

The proposed contribution is scientifically plausible but **not yet established as novel or feasible at full scope**.

The component ideas are established: hierarchical random-effects models, heavy-tailed likelihoods, multi-assay learning, domain generalization, study-grouped validation, conformal prediction, and aleatoric/epistemic uncertainty decomposition all have substantial prior art. A defensible PhytoTRUST contribution would therefore be their carefully constrained application to a different estimand:

> Given a compound–target pair's evidence available before time \(\tau\), predict the full distribution of a comparable measurement from a previously unseen independent study, while integrating rather than pretending to know the unseen study effect.

The potential gap is not “another bioactivity predictor.” It is the absence, in the reviewed literature, of a benchmarked model that jointly:

1. treats publication/study as the transfer domain;
2. separates compound–target potency from between-study and within-assay variation;
3. prohibits future-study outcomes and post-result protocol annotations from its inputs;
4. derives disagreement risk from a continuous posterior predictive distribution;
5. validates on completely unseen, later studies; and
6. audits natural products as a prespecified subgroup without using natural-product status as a shortcut feature.

This is a **candidate research gap**, not a novelty claim. Before publication, it requires a formal systematic review and expert review by assay scientists and statisticians.

The current PhytoTruth development set can support a small, blind-transfer pilot with study random intercepts. It cannot identify the proposed protocol-conditioned mean effects or study-by-chemistry random slopes: 112 of 141 development studies are singletons, only 10 studies contain at least three eligible events, the median compound has one prior measurement from one prior document, and key protocol fields are mostly missing. The full method should not be implemented until the data gate in this document is met.

---

## 1. Review method

Searches covered PubMed, publisher pages, DOI records, and primary conference proceedings. Query families combined *bioactivity/IC50/assay variability/reproducibility*, *multi-assay/multitask/QSAR*, *domain generalization/group shift*, *hierarchical/random effects*, *uncertainty/conformal/calibration*, and *natural products*. Bibliographic details and DOIs were checked against PubMed or publisher records where available. The review is critical and targeted, not PRISMA-complete; it should not be described as a systematic review.

### 1.1 Critical state of the art

#### A. Heterogeneous bioactivity is a measurement problem, not only a modeling problem

Kramer et al. estimated experimental uncertainty in heterogeneous public \(K_i\) data, while Kalliokoski et al. showed that mixed public IC50 values carry substantial variability and documented annotation errors and scarce protocol detail [1,2]. The latter reported a mixed-IC50 standard deviation of about 0.68 log units and found that more than 90% of initially paired records failed independence/quality filters. An AstraZeneca database analysis later found lower variability under controlled internal practice, emphasizing that provenance and laboratory process matter [3]. Apparent activity can also arise from compound-dependent assay interference rather than intended target biology [6].

These studies justify modeling study heterogeneity explicitly. They do **not** show that future study disagreement is predictable from chemical structure.

Pharmacogenomic comparisons provide the same warning at another biological scale: large studies have produced discordant drug-response measurements [7], while independent re-profiling showed improved reproducibility when laboratory and analytic practices were harmonized [8]. The implication is not that one paper “wins,” but that the data-generating process and preprocessing can dominate transfer.

#### B. Standard QSAR and multi-assay models estimate activity, not study transferability

Profile-QSAR and its successors exploit activity profiles across related assays [13–16]. All-Assay-Max2 pQSAR scaled this idea to thousands of assays. Bayesian matrix completion likewise borrows information across chemicals and assay endpoints and quantifies uncertainty [18]. Multitask neural networks and large ChEMBL comparisons demonstrate benefits from shared representation under some splits [11,12,17].

This literature overlaps strongly with any proposal to pool assays or targets. PhytoTRUST cannot claim novelty for multitask learning, matrix completion, assay embeddings, or cross-assay information sharing. The distinction must be that assays/studies are modeled as explicit sources of measurement variation and that predictions target a **new-study measurement distribution**, not a missing entry in a fixed assay matrix.

#### C. Validation strategy can dominate apparent performance

Time-split validation better approximates prospective QSAR performance than random splitting in Sheridan's study [9]. Wallach and Heifets showed that common ligand benchmarks can reward memorization [10]. MoleculeNet standardized scaffold-aware molecular benchmarks [11], and Mayr et al. used nested cluster validation at ChEMBL scale [12]. Domain-generalization and WILDS work formalize evaluation under unseen domains and natural distribution shifts [25,26]. Group DRO targets worst-group performance, but requires meaningful training groups and does not identify unseen group effects [27].

Therefore, random molecule splits are unacceptable here. Publication, molecule-event history, target, and time must be controlled simultaneously.

#### D. Uncertainty methods are useful but do not automatically survive study shift

Conformal prediction can provide distribution-free marginal coverage under exchangeability and has been applied to QSAR [19,20]. Test-time dropout plus conformal calibration, Bayesian graph convolution, neural uncertainty comparisons, and large-scale UQ benchmarks provide practical uncertainty estimators [21–24]. These works repeatedly show that uncertainty quality is empirical and can degrade out of domain. Conformal validity also depends on the exchangeability unit: calibrating on molecules does not guarantee coverage for an unseen study cluster.

PhytoTRUST must report proper scoring rules and study-clustered coverage, interval width, and subgroup behavior. An ensemble spread is not a confidence interval unless calibrated and validated as such.

#### E. Natural-product databases identify provenance, not transferability

NPASS provides quantitative activity and species-source information [29,30], and COCONUT 2.0 provides a curated open natural-product collection [31]. Natural-product membership can define an audited subgroup through standardized structure matching. It does not establish assay reliability, mechanism, efficacy, or a distinct statistical data-generating process. Natural-product status should therefore not be a predictor in the primary model; chemical-domain distance and evidence coverage should be modeled directly.

### 1.2 Verified bibliography and relevance

| # | Peer-reviewed work | DOI / record | Contribution and limit for PhytoTRUST |
|---:|---|---|---|
| 1 | Kramer C, Kalliokoski T, Gedeck P, Vulpetti A. *The Experimental Uncertainty of Heterogeneous Public Ki Data*. J Med Chem. 2012;55:5165–5173. | [10.1021/jm300131x](https://doi.org/10.1021/jm300131x) | Quantifies heterogeneous public-data error; does not predict unseen-study shifts. |
| 2 | Kalliokoski T, Kramer C, Vulpetti A, Gedeck P. *Comparability of Mixed IC50 Data—A Statistical Analysis*. PLoS ONE. 2013;8:e61007. | [10.1371/journal.pone.0061007](https://doi.org/10.1371/journal.pone.0061007) | Direct evidence of mixed-IC50 variability, filtering difficulty, and protocol scarcity. |
| 3 | Kramer C, Dahl G, Tyrchan C, Ulander J. *A Comprehensive Company Database Analysis of Biological Assay Variability*. Drug Discov Today. 2016;21:1213–1221. | [10.1016/j.drudis.2016.03.015](https://doi.org/10.1016/j.drudis.2016.03.015) | Controlled internal repeatability; warns against assuming public and controlled data have equal noise. |
| 4 | Fourches D, Muratov EN, Tropsha A. *Trust, But Verify: On the Importance of Chemical Structure Curation in Cheminformatics and QSAR Modeling Research*. J Chem Inf Model. 2010;50:1189–1204. | [10.1021/ci100176x](https://doi.org/10.1021/ci100176x) | Establishes structure/data curation as prerequisite. |
| 5 | Fourches D, Muratov EN, Tropsha A. *Trust, but Verify II: A Practical Guide to Chemogenomics Data Curation*. J Chem Inf Model. 2016;56:1243–1252. | [10.1021/acs.jcim.6b00129](https://doi.org/10.1021/acs.jcim.6b00129) | Chemogenomic curation guidance; not a study-transfer model. |
| 6 | Thorne N, Auld DS, Inglese J. *Apparent Activity in High-Throughput Screening: Origins of Compound-Dependent Assay Interference*. Curr Opin Chem Biol. 2010;14:315–324. | [10.1016/j.cbpa.2010.03.020](https://doi.org/10.1016/j.cbpa.2010.03.020) | Mechanistic assay-interference taxonomy; supports chemistry × detection interactions. |
| 7 | Haibe-Kains B et al. *Inconsistency in Large Pharmacogenomic Studies*. Nature. 2013;504:389–393. | [10.1038/nature12831](https://doi.org/10.1038/nature12831) | Shows cross-study drug-response discordance. |
| 8 | Haverty PM et al. *Reproducible Pharmacogenomic Profiling of Cancer Cell Line Panels*. Nature. 2016;533:333–337. | [10.1038/nature17987](https://doi.org/10.1038/nature17987) | Shows harmonized practices can improve agreement; relevant counterpoint to inevitable-noise claims. |
| 9 | Sheridan RP. *Time-Split Cross-Validation as a Method for Estimating the Goodness of Prospective Prediction*. J Chem Inf Model. 2013;53:783–790. | [10.1021/ci400084k](https://doi.org/10.1021/ci400084k) | Basis for temporal evaluation. |
| 10 | Wallach I, Heifets A. *Most Ligand-Based Classification Benchmarks Reward Memorization Rather than Generalization*. J Chem Inf Model. 2018;58:916–932. | [10.1021/acs.jcim.7b00403](https://doi.org/10.1021/acs.jcim.7b00403) | Shows benchmark leakage/similarity bias. |
| 11 | Wu Z et al. *MoleculeNet: A Benchmark for Molecular Machine Learning*. Chem Sci. 2018;9:513–530. | [10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A) | Standardized molecular tasks/splits; study transfer is not its estimand. |
| 12 | Mayr A et al. *Large-Scale Comparison of Machine Learning Methods for Drug Target Prediction on ChEMBL*. Chem Sci. 2018;9:5441–5451. | [10.1039/C8SC00148K](https://doi.org/10.1039/C8SC00148K) | Strong cluster-nested ML benchmark; does not decompose study noise. |
| 13 | Martin E et al. *Profile-QSAR: A Novel Meta-QSAR Method that Combines Activities across the Kinase Family*. J Chem Inf Model. 2011;51:1942–1956. | [10.1021/ci1005004](https://doi.org/10.1021/ci1005004) | Cross-assay activity profiles; major overlap with multi-assay representation. |
| 14 | Martin EJ et al. *Profile-QSAR 2.0: Kinase Virtual Screening Accuracy Comparable to Four-Concentration IC50s for Realistically Novel Compounds*. J Chem Inf Model. 2017;57:2077–2088. | [10.1021/acs.jcim.7b00166](https://doi.org/10.1021/acs.jcim.7b00166) | Prospective-leaning chemical novelty; not new-study variance prediction. |
| 15 | Martin EJ et al. *All-Assay-Max2 pQSAR: Activity Predictions as Accurate as Four-Concentration IC50s for 8558 Novartis Assays*. J Chem Inf Model. 2019;59:4450–4459. | [10.1021/acs.jcim.9b00375](https://doi.org/10.1021/acs.jcim.9b00375) | Massively multitask assay modeling; an essential baseline if data scale permits. |
| 16 | Martin EJ, Zhu X-W. *Collaborative Profile-QSAR: A Natural Platform for Building Collaborative Models among Competing Companies*. J Chem Inf Model. 2021. | [10.1021/acs.jcim.0c01342](https://doi.org/10.1021/acs.jcim.0c01342) | Collaborative cross-assay modeling; does not solve unseen-study identifiability. |
| 17 | Xu Y, Ma J, Liaw A, Sheridan RP, Svetnik V. *Demystifying Multitask Deep Neural Networks for Quantitative Structure–Activity Relationships*. J Chem Inf Model. 2017;57:2490–2504. | [10.1021/acs.jcim.7b00087](https://doi.org/10.1021/acs.jcim.7b00087) | Clarifies multitask gains and limits. |
| 18 | Jin B et al. *Bayesian Matrix Completion for Hypothesis Testing*. J R Stat Soc C. 2023;72:254–270. | [10.1093/jrsssc/qlac005](https://doi.org/10.1093/jrsssc/qlac005) | Hierarchical sparse chemical–assay prediction and UQ; fixed endpoint matrix rather than new-study transfer. |
| 19 | Norinder U, Carlsson L, Boyer S, Eklund M. *Introducing Conformal Prediction in Predictive Modeling for Regulatory Purposes*. Regul Toxicol Pharmacol. 2015;71:279–284. | [10.1016/j.yrtph.2014.12.021](https://doi.org/10.1016/j.yrtph.2014.12.021) | Validity framework under exchangeability. |
| 20 | Bosc N et al. *Large Scale Comparison of QSAR and Conformal Prediction Methods and Their Applications in Drug Discovery*. J Cheminform. 2019;11:4. | [10.1186/s13321-018-0325-4](https://doi.org/10.1186/s13321-018-0325-4) | Temporal ChEMBL comparison; supports prospective calibration tests. |
| 21 | Cortés-Ciriano I, Bender A. *Reliable Prediction Errors for Deep Neural Networks Using Test-Time Dropout*. J Chem Inf Model. 2019;59:3330–3339. | [10.1021/acs.jcim.9b00297](https://doi.org/10.1021/acs.jcim.9b00297) | Conformalized molecular intervals; calibration unit is not automatically the study. |
| 22 | Scalia G et al. *Evaluating Scalable Uncertainty Estimation Methods for Deep Learning-Based Molecular Property Prediction*. J Chem Inf Model. 2020;60:2697–2717. | [10.1021/acs.jcim.9b00975](https://doi.org/10.1021/acs.jcim.9b00975) | Compares ensembles, bootstrapping, and MC dropout; emphasizes OOD difficulty. |
| 23 | Hirschfeld L et al. *Uncertainty Quantification Using Neural Networks for Molecular Property Prediction*. J Chem Inf Model. 2020;60:3770–3780. | [10.1021/acs.jcim.0c00502](https://doi.org/10.1021/acs.jcim.0c00502) | Molecular UQ benchmark; not hierarchical study uncertainty. |
| 24 | Ryu S, Kwon Y, Kim WY. *A Bayesian Graph Convolutional Network for Reliable Prediction of Molecular Properties with Uncertainty Quantification*. Chem Sci. 2019;10:8438–8446. | [10.1039/C9SC01992H](https://doi.org/10.1039/C9SC01992H) | Separates data/model uncertainty in a molecular model; no publication effects. |
| 25 | Zhou K et al. *Domain Generalization: A Survey*. IEEE TPAMI. 2023;45:4396–4415. | [10.1109/TPAMI.2022.3195549](https://doi.org/10.1109/TPAMI.2022.3195549) | Taxonomy of unseen-domain methods; chemistry-specific measurement structure is absent. |
| 26 | Koh PW et al. *WILDS: A Benchmark of in-the-Wild Distribution Shifts*. ICML/PMLR 139. 2021. | [Proceedings record](https://proceedings.mlr.press/v139/koh21a.html) (no DOI located) | Establishes domain-shift evaluation discipline. |
| 27 | Sagawa S et al. *Distributionally Robust Neural Networks for Group Shifts*. ICLR. 2020. | [OpenReview record](https://openreview.net/forum?id=ryxGuJrFvS) (no DOI located) | Worst-group training baseline; requires sufficiently sampled, meaningful study groups. |
| 28 | Yang K et al. *Analyzing Learned Molecular Representations for Property Prediction*. J Chem Inf Model. 2019;59:3370–3388. | [10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237) | Strong message-passing/descriptor benchmark across public and industrial data. |
| 29 | Zeng X et al. *NPASS: Natural Product Activity and Species Source Database*. Nucleic Acids Res. 2018;46:D1217–D1222. | [10.1093/nar/gkx1026](https://doi.org/10.1093/nar/gkx1026) | Structure/activity/source provenance for NP verification. |
| 30 | Zeng X et al. *NPASS Database Update 2023*. Nucleic Acids Res. 2023;51:D621–D628. | [10.1093/nar/gkac1069](https://doi.org/10.1093/nar/gkac1069) | Updated NP coverage; membership is not reliability evidence. |
| 31 | Chandrasekhar V et al. *COCONUT 2.0: A Comprehensive Overhaul and Curation of the Collection of Open Natural Products Database*. Nucleic Acids Res. 2025;53:D634–D643. | [10.1093/nar/gkae1063](https://doi.org/10.1093/nar/gkae1063) | Open NP structure provenance and curation. |
| 32 | Gaulton A et al. *The ChEMBL Database in 2017*. Nucleic Acids Res. 2017;45:D945–D954. | [10.1093/nar/gkw1074](https://doi.org/10.1093/nar/gkw1074) | Primary source for ChEMBL's curated bioactivity model. |

---

## 2. Specific research gap and bounded novelty claim

### 2.1 What current methods do not adequately address

1. **Mean prediction is conflated with transferability.** Standard QSAR asks for expected activity. A good mean predictor can still be badly calibrated for a new laboratory or protocol.
2. **Assay identity is often a task label, not a stochastic source of error.** Multi-assay methods learn correlations among observed assays. They generally do not integrate an unobserved future-study effect when scoring a new study.
3. **Heterogeneity analyses are mainly descriptive.** IC50/Ki reproducibility studies quantify aggregate variation, but do not produce compound-specific new-study predictive distributions conditioned only on historical evidence.
4. **Generic domain generalization ignores the measurement hierarchy.** A publication contains correlated compounds and assay choices; molecules, targets, assays, documents, and time are crossed or nested. Treating rows as independent domains gives anticonservative inference.
5. **Molecular UQ is not study UQ.** Chemical-distance or ensemble variance can detect chemical extrapolation while missing an unseen laboratory-wide shift.
6. **Protocol-informed retrospective models can leak.** ChEMBL assay descriptions are deposited with results. They are not automatically evidence that a field was known before measurement.
7. **Natural products are underpowered and confounded.** NP membership correlates with scaffold, physicochemical range, target choice, evidence age, and reporting practice. A binary NP feature can encode dataset construction rather than biology.

### 2.2 Proposed differentiator

PhytoTRUST would define and benchmark **new-study posterior prediction** with explicit information timing. The core differentiator is the combination of:

- a historical evidence-set encoder;
- a hierarchical heavy-tailed measurement model;
- study-level random intercepts and, only if identifiable, low-rank study-by-chemistry effects;
- integration over the random effect of a completely unseen study;
- separately reported epistemic, between-study, and residual uncertainty;
- study-level calibration and evaluation; and
- a derived probability/magnitude of disagreement instead of a directly tuned binary heuristic.

Each component has precedent. The potentially publishable contribution is the estimand, leakage contract, integrated design, natural-product audit, and empirical benchmark. “First,” “novel,” or “superior” must not be used unless a formal review and independent experiment support it.

### 2.3 Novelty comparison

| Method family | Learns potency | Shares across assays | Explicit unseen-study random effect | Timing-safe historical evidence | Study-clustered UQ | NP subgroup audit |
|---|---:|---:|---:|---:|---:|---:|
| Single-task QSAR | Yes | No | No | Usually not central | Rarely | Rarely |
| Multitask / pQSAR | Yes | Yes | Usually no; observed assay is a task | Variable | Usually molecule-level | No |
| Matrix completion | Yes | Yes | No; predicts cells in a defined matrix | Variable | Bayesian in some models | No |
| Mixed-effects assay analysis | Sometimes | Limited | Yes | Usually retrospective | Yes | No |
| Generic domain generalization / group DRO | Yes | Across groups | No probabilistic study decomposition | Depends on design | Not inherent | No |
| Molecular conformal/UQ | Yes | Optional | No | Depends on split | Usually row/molecule-level | Rarely |
| **PhytoTRUST proposal** | Yes, as latent historical location | Yes | **Integrated out for new studies** | **Explicit filtration** | **Study is calibration/resampling unit** | **Prespecified, not a feature** |

---

## 3. Proposed methodology

### 3.1 Prediction tasks and information sets

Let \(c\) denote compound, \(t\) target/context, \(s\) independent study/publication, \(a\) assay, and \(r\) replicate/record. Let \(Y_{ctsar}\) be an exact, normalized pChEMBL value for one compatible endpoint/context.

Define the information filtration immediately before a future study result at time \(\tau\):

\[
\mathcal F^-_{ct,\tau}=\{X_c,Z_t,H_{ct,<\tau},M_{<\tau}\},
\]

where:

- \(X_c\): standardized structure and molecular representation;
- \(Z_t\): target and prespecified compatible assay-context descriptors;
- \(H_{ct,<\tau}\): all allowed earlier measurements with activity, assay/document identifiers, dates, and explicitly reported metadata;
- \(M_{<\tau}\): database/version and curation information available before \(\tau\).

Two estimands must remain separate:

1. **Blind-study transfer (primary with current data):** predict \(p(Y_{new}\mid\mathcal F^-_{ct,\tau})\), integrating over the unknown future protocol and study effect.
2. **Registered-protocol transfer (future extension):** predict \(p(Y_{new}\mid\mathcal F^-_{ct,\tau},Q^-_{new})\), where \(Q^-_{new}\) is a timestamped assay plan known before results. Current ChEMBL descriptions do not prove such timing and cannot be used for this estimand.

The primary use case is a previously measured compound–target pair entering a new independent study. A cold-compound task with no \(H_{ct,<\tau}\) is a distinct secondary task.

### 3.2 Hierarchical generative model

Historical evidence is represented as a set, not only a median:

\[
h_{ct}=\rho\left(\left\{e_\omega(Y_j,Q_j,\Delta t_j,P_j) : j\in H_{ct,<\tau}\right\}\right),
\]

where \(P_j\) contains provenance and missingness indicators, \(e_\omega\) is a record encoder, and \(\rho\) is a permutation-invariant robust pooling operator. A median/MAD implementation is the mandatory simple version; a learned DeepSets/attention encoder is optional only with adequate sample size.

The latent historical location is

\[
\mu_{ct}=\operatorname{median}(H_{ct,<\tau}) + f_\theta(X_c,Z_t,h_{ct}).
\]

For an observed training study:

\[
Y_{ctsar}\sim t_\nu\!\left(
\mu_{ct}+g_\psi(Q^-_s,Z_t)+u_{st}+v_s^\top W\phi(X_c,Z_t),\;
\sigma_{ctsar}
\right),
\]

\[
u_{st}\sim\mathcal N(0,\tau^2_t),\qquad
v_s\sim\mathcal N(0,\Sigma_v),
\]

\[
\log\sigma_{ctsar}=q_\eta(X_c,Z_t,h_{ct},Q^-_s,R_{ct}),
\]

where:

- \(g_\psi\) models only pre-result protocol fields;
- \(u_{st}\) is a target-specific study shift;
- \(v_s^\top W\phi(\cdot)\) is an optional low-rank study-by-chemistry interaction;
- \(R_{ct}\) contains prior evidence coverage and missingness;
- Student-\(t\) residuals reduce sensitivity to remaining transcription errors/outliers.

For a genuinely unseen study, \(u_{new,t}\) and \(v_{new}\) are **not estimated from its outcomes**. They are integrated over their learned population distributions. If a “held-out study” contains several test compounds, using some of those outcomes to estimate its effect is transductive adaptation and must be reported as a separate sequential-update task, never mixed into blind transfer.

### 3.3 Stage-gated complexity

The model must be developed in stages:

- **M0:** robust target/context hierarchical model with study random intercept only.
- **M1:** M0 plus prior-history set features and heteroscedastic residual scale.
- **M2:** M1 plus prospectively timestamped protocol predictors.
- **M3:** M2 plus low-rank study-by-chemistry random slopes.

M2 is prohibited with the current retrospective held-out ChEMBL descriptions. M3 is prohibited until enough multi-compound studies exist to identify within-study response patterns.

### 3.4 Objective

For training studies \(\mathcal S_{train}\), use a study-balanced proper likelihood:

\[
\mathcal L_{NLL}=
\frac{1}{|\mathcal S_{train}|}
\sum_{s\in\mathcal S_{train}}
\frac{1}{n_s}\sum_{i\in s}
-\log p_\Theta(Y_i\mid \mathcal F_i^-)
+\lambda\,\Omega(\Theta).
\]

Equal study weighting prevents a single large publication from defining the solution. Hierarchical priors/penalties \(\Omega\) regularize study variance and optional interaction rank. A group-DRO/CVaR study loss may be evaluated as a baseline or sensitivity analysis, not assumed superior.

Hyperparameters are selected by inner study-grouped, forward-time validation using mean study-level negative log predictive density (NLPD), with CRPS as a co-primary diagnostic. Binary AUROC must not drive model selection.

### 3.5 Derived disagreement outcomes

Let \(A_{ct}=\operatorname{median}(H_{ct,<\tau})\) and \(\delta=1.0\) pChEMBL. The model returns:

\[
E[D_{new}]=E\left[|Y_{new}-A_{ct}|\mid\mathcal F^-\right]
\]

and

\[
P(D_{new}\ge\delta)=
P\left(|Y_{new}-A_{ct}|\ge1.0\mid\mathcal F^-\right).
\]

Thus binary disagreement is a transparent functional of the continuous predictive distribution, not an independently tuned classifier. The model also returns predictive quantiles and a decomposition of uncertainty.

### 3.6 Uncertainty decomposition

Using the law of total variance, report:

\[
\operatorname{Var}(Y_{new}\mid\mathcal F^-)=
\underbrace{\operatorname{Var}_{\Theta}\!\left[E(Y_{new}\mid\Theta,\mathcal F^-)\right]}_{\text{epistemic}}
+\underbrace{E_{\Theta}[\tau_t^2+\phi^\top W^\top\Sigma_vW\phi]}_{\text{between-study}}
+\underbrace{E_{\Theta}\!\left[\frac{\nu}{\nu-2}\sigma^2\right]}_{\text{residual/within-assay}},
\]

for \(\nu>2\). If protocol fields are missing, their integration contributes a separately labeled **protocol-missingness** component when estimable. Components are variance-model terms, not causal attributions.

Posterior/deep-ensemble uncertainty may represent epistemic uncertainty. Study-level split conformal calibration can be applied on untouched calibration studies to check or adjust predictive intervals. Any conformal claim must specify whether coverage is marginal over studies, targets, or observations and must report finite-cluster limitations.

### 3.7 Algorithm pseudocode

```text
INPUT:
  normalized exact measurements with compound, target/context, study, assay, date
  standardized structures and NPASS/COCONUT provenance
  protocol fields with explicit availability timestamps

1. Freeze targets, compatible contexts, endpoint, threshold, and study definition.
2. Build rolling-origin events:
     H = evidence strictly earlier than held-out study s
     y = measurements from s
   Reject assay/document overlap, incompatible contexts, duplicated literature values,
   and any feature not provably available before y.
3. Partition outer folds by publication time and connected study–compound components.
4. In each outer fold:
     a. Fit M0–M3 only where the data gate permits.
     b. Select hyperparameters in inner study-disjoint temporal folds by study-macro NLPD.
     c. Integrate unseen-study random effects; do not infer them from held-out outcomes.
     d. Calibrate intervals using calibration studies only.
     e. Predict the full distribution for each untouched study event.
5. Derive expected absolute shift and P(shift >= 1.0) from that distribution.
6. Score study-macro NLPD, CRPS, MAE, calibration, interval score, and binary metrics.
7. Repeat analyses for verified natural products without fitting or tuning on subgroup labels.
8. Bootstrap/permutation-resample studies, never individual rows, for inference.
9. Freeze one final candidate before acquiring a genuinely new independent test set.
```

---

## 4. Falsifiable scientific hypotheses

### Primary hypothesis

**H1:** On unseen, later publications, a stage-appropriate PhytoTRUST model will improve study-macro probabilistic accuracy—lower NLPD and CRPS—relative to (a) a historical-median model with global residual variance and (b) a robust target/context Student-\(t\) random-intercept model, without worsening calibrated interval score.

Support requires a paired study-cluster 95% confidence interval for the NLPD and CRPS difference entirely below zero in a new, untouched outer evaluation. AUROC improvement alone is insufficient.

### Mechanistic hypothesis

**H2:** Explicitly integrating unseen-study variance prevents study-level batch shifts from being misattributed to compound structure or historical potency. Therefore the gain, if real, will be larger in studies with multiple compounds and documented protocol differences than in singleton studies.

This is falsified if the random-effect model does not improve proper scores/calibration, or if gains vanish under leave-target-out validation.

### Natural-product hypothesis

**H3:** For structure-verified natural products, history-aware uncertainty and chemical-domain distance will improve calibration and interval score relative to a global model, even if discrimination of binary disagreement does not improve.

This is falsified if NP coverage error, normalized interval score, or NLPD is no better than the baseline in a prespecified untouched NP subgroup. No minimum AUROC is promised.

### Protocol hypothesis

**H4:** Prospectively timestamped protocol variables explain transferable between-study variation beyond generic missingness indicators.

This cannot be tested with current post-deposition ChEMBL descriptions. It becomes testable only with protocol information demonstrably recorded before outcomes. Missingness indicators alone do not count as support.

---

## 5. Experimental validation plan

### 5.1 Data sources

1. **ChEMBL:** exact IC50 pChEMBL measurements, structures, target/context, assay, document, organism, BAO format, descriptions, and dates.
2. **NPASS and COCONUT:** exact standardized-structure matches defining natural-product provenance; no reliability label is inferred.
3. **External replication candidates:** BindingDB or PubChem BioAssay records linked to original publications, after independent curation and deduplication. They must not be treated as independent if sourced from the same publication or copied from ChEMBL.
4. **Preferred prospective source:** preregistered assay plans and blinded replicate measurements from collaborating laboratories. This is required to test registered-protocol transfer cleanly.

### 5.2 Eligibility

- Exact numeric IC50 only for the primary benchmark; Ki and other endpoints are separate strata.
- Same target identity, organism, target type, endpoint, and prespecified compatible assay context.
- Standardized structures with stereochemistry retained; salts handled reproducibly.
- At least one genuinely earlier independent document for repeat-compound transfer.
- Assay and document IDs disjoint between history and outcome.
- Suspected duplicated literature values, unit errors, ambiguous targets, and incompatible biological systems excluded with reason codes.
- Protocol fields remain unknown when unreported; no NLP inference of experimental facts.

### 5.3 Nested leakage-resistant validation

**Outer test:** rolling publication-time blocks. Every held-out publication is unseen. Connected components of publication and constructed event prevent the same outcome event appearing in development. A compound's earlier measurements are allowed only as the intended historical input; later outcomes are forbidden.

**Inner selection:** forward-chaining, study-disjoint folds. All preprocessing, representation fitting, target pooling, calibration, and hyperparameter selection occur inside the training fold.

**Stress tests:**

- leave-one-target-out;
- leave-one-study-out sensitivity;
- chemical-scaffold distance strata;
- studies with one versus multiple compounds;
- protocol-completeness strata;
- natural products, locked before scoring;
- cold-compound transfer, reported separately from repeat-compound transfer.

The already examined 319-event evaluation partition, 90-case independent cohort, and 94-case expanded cohort cannot be used for method selection or advertised as new independent validation.

### 5.4 Strong baselines

1. Historical median with global target-specific residual distribution.
2. Historical median with robust MAD/document-range uncertainty.
3. Robust Student-\(t\) target/context hierarchical model with study random intercept (M0).
4. Existing PhytoTruth logistic disagreement model and magnitude RF, frozen.
5. Single-task ECFP Random Forest/XGBoost regressor.
6. Multitask molecular model with assay/target labels.
7. pQSAR-like profile baseline where scale permits.
8. Bayesian matrix-factorization/completion baseline.
9. Deep/bootstrapped ensemble with no explicit study effects.
10. Group-DRO study model and ordinary empirical-risk minimization with the same architecture.
11. Study-cluster conformalized baseline.

Every complex method must beat M0, not only the potency median.

### 5.5 Metrics

**Primary probabilistic metrics:** study-macro NLPD and CRPS.  
**Point accuracy:** study-macro MAE and RMSE.  
**Calibration/sharpness:** empirical coverage at 50%, 80%, 90%, and 95%; coverage error; mean width; weighted/normalized interval score; calibration plots by study and target.  
**Uncertainty utility:** correlation of predicted scale with absolute error, risk–coverage curves, and error retention after abstention.  
**Derived binary outcome:** AUROC, PR-AUC, Brier score, calibration intercept/slope, and Precision@K for \(|\Delta pChEMBL|\ge1\).  
**Robustness:** worst-target and lower-decile study performance, plus natural-product subgroup metrics.

Micro-averaged metrics are secondary because large publications otherwise dominate.

### 5.6 Ablations

- remove study random effects;
- Gaussian versus Student-\(t\) likelihood;
- homoscedastic versus evidence-dependent scale;
- summary statistics versus evidence-set encoder;
- remove molecular representation;
- remove prior assay/provenance metadata;
- protocol values versus protocol missingness only, when prospectively available;
- random intercept versus low-rank study-by-chemistry interaction;
- ordinary ERM versus equal-study weighting versus group DRO;
- ensemble UQ versus hierarchical UQ versus combined;
- without cluster-conformal calibration;
- remove target pooling;
- NP status as evaluation-only versus an intentionally tested sensitivity model containing it (never primary).

### 5.7 Statistical analysis

- Prespecify primary candidate, metrics, and direction before outer-test access.
- Use paired hierarchical bootstrap resampling of studies within target; never bootstrap measurement rows independently.
- Report 95% confidence intervals for metric differences, not only separate model intervals.
- Use study-level paired permutation or Wilcoxon tests for study-macro losses when assumptions are reasonable.
- Correct the family of secondary comparisons with Holm's procedure.
- Report effect sizes and confidence intervals; a nominal \(p<0.05\) alone does not establish usefulness.
- Use simulation-based power analysis based only on development variance components before setting acquisition targets.
- Natural-product analyses are confirmatory only if their sample-size target is frozen and met; otherwise label them exploratory.

### 5.8 Success and stop rules

Proceed to a new independent test only if, in development-only nested validation:

1. PhytoTRUST improves both NLPD and CRPS over M0 with paired study-cluster intervals excluding zero;
2. interval coverage is near nominal without excessive width according to interval score;
3. direction is consistent under time, study, and target stress tests;
4. no severe NP calibration degradation is observed; and
5. variance components are stable across folds and not collapsed to zero or dominated by one study.

If these fail, the defensible output is a calibrated hierarchical baseline or a statement that the variation is not predictable from available pre-study information.

---

## 6. Feasibility audit of existing PhytoTruth data

The audit below reads existing artifacts only; no dataset or model was changed.

### 6.1 Available development evidence

| Quantity | Existing pre-2020 development partition |
|---|---:|
| Events | 212 |
| Unique molecules | 211 |
| Targets | 5 |
| Held-out studies/publications | 141 |
| Binary disagreements | 50 |
| Verified natural products | 47 |
| NP disagreements | 14 |
| Singleton held-out studies | 112/141 |
| Studies with at least 3 eligible events | 10 |
| Largest study | 25 events |
| Median prior measurements per event | 1 |
| Median prior documents per event | 1 |

The 319-event 2020+ partition is already examined and remains locked. It contains a 131-compound study and cannot be treated as fresh evidence. The first independent cohort has 90 cases but only seven disagreements and five verified NPs. The expanded cohort has 94 cases, 24 disagreements, 58 studies, 20 verified NPs, and seven NP disagreements; it has also already been evaluated and cannot be used for model development or fresh validation.

### 6.2 Protocol completeness

In the 212-event development partition:

- descriptions are present, but prior structured assay parameters are absent;
- mean direct-field completeness is 0.354;
- methodology is detected/reported for 17.0% of events;
- detection method for 25.5%; endpoint text for 19.8%;
- temperature for 3.3%; pH for 2.4%; duration for 9.0%.

In the 94-event expansion, temperature, pH, strain, parameter type, and several biological-system comparisons are unknown for essentially all events; 17 of 24 positive disagreements remain unresolved because protocol metadata are incomplete. “Unknown” cannot be recoded as “same.”

### 6.3 Identifiability problems

1. **Study effects versus target/time:** only five development targets and uneven time coverage make these effects strongly confounded.
2. **Random slopes:** singleton studies provide no within-study compound contrast. With only 10 studies containing three or more events, study-by-chemistry interactions are not identifiable.
3. **Heteroscedasticity:** one prior measurement for the median event provides almost no compound-specific historical variance information.
4. **Protocol effects:** descriptions deposited with results are not demonstrably pre-result covariates; most decisive fields are missing.
5. **Publication is not laboratory:** author overlap, shared protocols, and copied database values can violate independence.
6. **Selection bias:** repeated compounds and published IC50 values are not a random sample of experiments; extremes may be preferentially reported.
7. **NP confounding:** 47 development NPs and 14 positives are too few to distinguish NP-specific behavior from target, scaffold, age, and coverage.
8. **Truth is latent:** disagreement does not identify which study is correct. Without reference standards or controlled replication, “true potency,” study bias, and residual error are only partially identifiable.

### 6.4 Data gate for implementation

The current data are sufficient only for M0/M1 feasibility analysis with strong shrinkage. Before implementing M2/M3, acquire or curate approximately:

- at least 20 targets with compatible repeated endpoints;
- at least 200 independent studies in development, with enough multi-compound clusters for variance estimation;
- at least 30–50 studies containing five or more shared compound–target measurements across studies;
- substantially more than 100 verified NP events and enough NP disagreements for simulation-based subgroup power;
- timestamped pre-result protocol fields, especially substrate/concentration, construct/variant, cell system, detection technology, incubation time, temperature, pH, and laboratory;
- reference compounds or within-study controls linking studies; and
- one newly acquired, label-quarantined external or prospective evaluation set.

These are design targets, not guarantees. Final sample sizes must come from simulation-based power and parameter-recovery studies using development-only variance estimates.

### 6.5 Feasibility verdict

**Full PhytoTRUST: not currently identifiable and therefore not ready to implement.**

**Defensible next research step:** curate a study-by-compound matrix with multi-compound independent studies and timestamped protocols, then perform parameter-recovery simulations and fit M0 only. Escalate to M1–M3 solely when prespecified identifiability checks pass.

If additional data cannot be obtained, the scientifically honest endpoint is a study-aware uncertainty baseline that widens predictions for unseen studies. It may improve calibration while offering little or no discrimination of which compound will disagree. That would still be a useful negative result, but not evidence that assay disagreement is predictably transferable.

---

## 7. Risks and safeguards

| Risk | Consequence | Safeguard |
|---|---|---|
| Post-result protocol leakage | Inflated apparent transfer | Timestamp contract; blind-transfer primary estimand |
| One large publication dominates | Learns batch identity | Equal-study loss; study-cluster inference |
| Overparameterized hierarchy | Unidentified variance/slopes | Stage gates, shrinkage, parameter-recovery simulation |
| Protocol missingness treated as absence | False scientific facts | Explicit unknown category and provenance span |
| Natural-product shortcut | Dataset-label discrimination | NP status excluded from primary predictors |
| Chemical memorization | Optimistic accuracy | temporal, scaffold-distance, and target stress tests |
| Transductive held-out adaptation | Leaks other test outcomes | separate blind and sequential-update tasks |
| Ensemble spread called confidence | Misleading uncertainty | proper calibration and coverage tests |
| Binary threshold instability | Information loss | continuous magnitude/distribution primary |
| Database copies called independent | Pseudoreplication | publication/assay/activity provenance and duplicate audit |
| Negative result hidden | Unsupported novelty claim | frozen hypotheses and mandatory baseline reporting |

## Final methodological position

PhytoTRUST should be framed as **study-aware predictive measurement modeling**, not a claim to discover biological truth. Its most defensible innovation would be an explicitly timed, hierarchical predictive distribution for an unseen study and a benchmark that exposes whether molecular/history signals add anything beyond honest between-study uncertainty. The current evidence supports designing that experiment, but not claiming that the proposed method will work, is novel in all components, or is ready to build.
