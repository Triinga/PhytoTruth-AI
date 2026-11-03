# PhytoTRUST-M1: evidence-conditioned study transferability

## Status and scope

PhytoTRUST-M1 is a development-only methodological experiment for predicting the
distribution of a comparable IC50 pChEMBL measurement in a previously unseen
publication. It does not estimate drug efficacy, safety, clinical success, or an
unobservable biological “truth.” The implementation is valid and reproducible;
methodological novelty and superiority remain empirical questions.

The experiment uses only the connected CHEMBL222 development partition and assay
metadata for its 67 development assay identifiers. The evaluation history, future
index, and quarantined labels are prohibited inputs. Natural-product membership is
used only for subgroup analysis.

## Relationship to prior work

M1 does not claim novelty for hierarchical regression, Student-t likelihoods,
heteroscedastic regression, multi-assay modeling, domain generalization, or
uncertainty quantification. Important overlap includes:

- Kramer et al., heterogeneous public bioactivity uncertainty,
  DOI `10.1021/jm300131x`.
- Kalliokoski et al., comparability and uncertainty of public IC50 values,
  DOI `10.1371/journal.pone.0061007`.
- Sheridan, temporal validation for prospective molecular prediction,
  DOI `10.1021/ci400084k`.
- Mayr et al., cluster-nested ChEMBL-scale validation,
  DOI `10.1039/C8SC00148K`.
- Kendall and Gal, decomposition of predictive uncertainty,
  DOI `10.48550/arXiv.1703.04977`.
- Kuleshov et al., calibrated regression uncertainty,
  DOI `10.48550/arXiv.1807.00263`.

The potentially differentiating, testable contribution is narrower: a strict
information-timed new-study estimand in which an historical evidence set controls
both the expected transfer correction and residual uncertainty, while the random
effect of an entirely unseen publication is integrated rather than estimated from
its outcomes. The benchmark tests whether that integration adds value beyond M0,
homoscedastic Student-t regression, median potency, and a nonlinear tabular model.
The present experiment cannot establish that this combination is novel.

## Information filtration

For compound \(c\), target/context \(t\), held-out publication \(s\), and time
\(\tau_s\), define

\[
\mathcal F^-_{c,t,s} = \{Y_{c,t,j}, Q_j, P_j, X_c: j \ne s\},
\]

for retrospective study-grouped validation and

\[
\mathcal F^-_{c,t,s} = \{Y_{c,t,j}, Q_j, P_j, X_c:
\tau_j < \tau_s\}
\]

for forward-temporal validation. Here \(Q_j\) contains literal protocol facts
reported for allowed historical assays, \(P_j\) contains their provenance and
missingness, and \(X_c\) contains deterministic molecular descriptors.

No activity value, assay ID, publication ID, or protocol field belonging to the
held-out publication is included in \(\mathcal F^-\). Same-compound historical
measurements are intentionally allowed because transfer for a previously measured
compound is the declared estimand.

## Historical evidence representation

The robust historical center is

\[
m_c = \operatorname{median}_{j\in\mathcal F^-} Y_{c,t,j}.
\]

The fixed, prespecified location feature vector contains historical mean-minus-
median, historical MAD, log evidence count, scaled molecular weight, and scaled
logP. The scale feature vector contains historical MAD, inverse square-root evidence
count, historical protocol missingness, protocol-category diversity, and historical
year span.

Protocol features come from deterministic literal rules applied to ChEMBL assay
descriptions and direct fields. Unreported fields remain missing. No absence of a
text match is interpreted as proof that a protocol characteristic was absent.

## Generative model

For transfer event \(i=(c,t,s)\):

\[
\mu_i = m_c + \alpha + \beta^\top \widetilde x_i,
\]

\[
\log \sigma_i = \log \sigma_0 + \gamma^\top \widetilde r_i,
\]

\[
u_s \sim \mathcal N(0,\tau^2),
\qquad
Y_i \mid u_s \sim t_\nu(\mu_i + u_s, \sigma_i),
\]

with \(\nu=4\). Tildes denote standardization using training-fold statistics only.
For blind transfer, \(u_{new}\) is integrated over
\(\mathcal N(0,\tau^2)\); it is never estimated from held-out outcomes.

The predictive variance, when \(\nu>2\), is

\[
\operatorname{Var}(Y_{new}\mid\mathcal F^-)
= \tau^2 + \sigma_i^2\frac{\nu}{\nu-2}.
\]

The first term is labeled between-study uncertainty and the second
evidence-conditioned residual uncertainty. These are model components, not causal
attributions. Prediction quantiles and likelihoods integrate the Normal random
effect numerically by 15-node Gauss-Hermite quadrature.

## Objective and regularization

Let \(S_{train}\) be the training publications and \(n_s\) their event counts. M1
minimizes the equal-publication-weighted marginal negative log likelihood

\[
\mathcal L(\Theta)=
\frac{1}{|S_{train}|}\sum_{s\in S_{train}}
\frac{-1}{n_s}\log
\int \prod_{i\in s}
p(Y_i\mid u_s,\mathcal F^-_i,\Theta)p(u_s)\,du_s
+ \frac{\|\beta\|_2^2}{2a^2N}
+ \frac{\|\gamma\|_2^2}{2b^2N},
\]

where \(a=0.35\), \(b=0.5\), and \(N\) is the number of training events. These
values and all feature lists were frozen before final validation. Equal publication
weighting prevents a large compound series from defining the fit.

## Baselines and ablations

- Median potency: historical median with a training-only Student-t scale.
- Student-t regression: historical median plus global correction, no study effect.
- XGBoost: shallow, strongly regularized residual regression using the same
  prediction-time features. Its uncertainty scale is estimated from internal
  publication-grouped cross-fitted residuals.
- M0: homoscedastic Student-t model with an unseen-publication random intercept.
- M1 without protocol features: removes historical protocol missingness/diversity.
- Homoscedastic M1: removes all evidence-conditioned scale features.
- M1 without study variance: fixes \(\tau=0\).

These ablations identify whether any gain is attributable to historical protocol
signals, heteroscedasticity, or explicit study variance rather than the shared
historical median.

## Validation

Study-grouped validation leaves out a whole publication, activity IDs, and assay
IDs. It is retrospective because other development publications may occur later.
Forward-temporal validation trains only on strictly earlier publication years.
Publication is the resampling unit for paired bootstrap confidence intervals.

Reported metrics are MAE, RMSE, negative log likelihood, 80% and 90% predictive-
interval coverage and width, and mean absolute coverage error. Natural products and
one-versus-multiple prior-study support are prespecified subgroups.

M1 is considered development-supported only if it improves NLL over the strongest
baseline in both validation schemes with both publication-bootstrap 95% intervals
excluding zero, without worse interval calibration or natural-product degradation.
Failure prevents access to the quarantined final test.

## Algorithm pseudocode

```text
INPUT: development measurements D, development assay metadata A
ASSERT: all rows are CHEMBL222 development rows before 2020
ASSERT: A contains exactly the assay IDs referenced by D
AGGREGATE D to compound-publication medians, preserving provenance
EXTRACT literal protocol tokens and explicit missingness from A

FOR each outer validation fold:
    separate whole held-out publications
    if temporal: discard every contemporaneous or later training publication
    for each eligible compound-publication training edge:
        construct history excluding that publication
        compute robust evidence, molecular, protocol, and missingness features
    for each held-out edge:
        compute the same features from allowed training history only
    verify held-out publication/activity/assay IDs do not enter training/history
    fit median, Student-t, XGBoost, M0, M1, and ablations
    integrate unseen-study effect for M0/M1 predictive distributions
    save paired predictions and uncertainty components

COMPUTE metrics by scheme and subgroup
BOOTSTRAP paired differences by held-out publication
FIT final development-only M1 for parameter reporting
WRITE hashes, predictions, metrics, figures, and scientific report
DO NOT open quarantined final outcomes unless the frozen support rule passes
```

## Assumptions and identifiability risks

1. Publications approximate independent studies, which may fail when laboratories,
   protocols, or copied values span documents.
2. The compatible ChEMBL context reduces but cannot eliminate protocol variation.
3. Literal descriptions are incomplete and retrospective; they are not registered
   experimental plans.
4. Only 29 compounds connect development studies and one compound appears in three
   studies. Mean and scale effects are consequently weakly identifiable.
5. Historical protocol missingness may reflect reporting practice rather than
   experimental reliability.
6. Exponential heteroscedastic models can extrapolate sharply when temporal feature
   values lie outside their narrow training range. This risk is evaluated rather
   than hidden by post-validation clipping.
7. Natural-product results are exploratory because only three connected verified
   natural products are available.

## Interpretation rule

Positive results would justify one new untouched evaluation, not a novelty claim.
Negative results indicate that available historical evidence or protocol metadata
does not support reliable evidence-conditioned transferability learning. Such a
result must not be repaired by tuning on the quarantined cohort.
