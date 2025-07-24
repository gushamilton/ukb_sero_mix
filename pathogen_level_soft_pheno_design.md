# Pathogen-level Soft Seropositivity Probabilities

This note documents the design decisions for **script 04** (to be added next) which will generate per-pathogen posterior infection probabilities by combining multiple antigen measurements collected in UK Biobank Multiplex Serology.

---

## 1  Motivation

The DKFZ/UKB validation paper provides hard cut-offs and logical rules (e.g. “positive for ≥2 of 4 antigens means EBV-seropositive”).  These rules ignore:

* continuous antibody titre information
* uncertainty near the cut-offs
* correlation between antigens

Script 02 already fits *per-antigen* two-component skew-t mixtures and outputs a soft probability
\( p_{ik}=P(\text{seropositive}\mid\text{antigen }k) \).

The new goal is to move **within a pathogen**, combining its antigens into a *single* posterior
\[ p_{j}^{\text{pathogen}} = P(\text{truly infected}\mid \mathbf y_j) \]
for each individual \(j\), without borrowing information across pathogens.

---

## 2  Statistical model

### 2.1  Latent-class mixture

For a pathogen with \(m\) antigens and participant \(j\), let

* \(\mathbf y_j = (y_{j1},\ldots,y_{jm})\) be the log-MFI values (columns `*_IgG_raw`).
* \(z_j \in \{0,1\}\) denote the unobserved true infection state.

We assume a **two-class multivariate mixture** with *full* covariance matrices:

\[
\mathbf y_j \mid z_j = 0 \;\sim\; \mathcal N(\boldsymbol\mu_0, \mathbf\Sigma_0),\qquad
\mathbf y_j \mid z_j = 1 \;\sim\; \mathcal N(\boldsymbol\mu_1, \mathbf\Sigma_1).
\]

Key points
1. \(\mathbf\Sigma_z\) is unrestricted → antigen correlations are learned, no independence assumption needed.
2. A normal mixture is sufficient; heavy tails were already handled at per-antigen stage.
3. The mixture weight \(\pi = P(z=1)\) is estimated from the data (no cross-pathogen prior).

### 2.2  Initialisation using Script 02 output

Unsupervised EM can stumble, so we seed with the per-antigen probabilities:

1. Compute a provisional infection label
   \[ \hat z_j^{(0)} = \mathbb 1\Bigl[\frac1m\sum_{k=1}^m p_{jk} > 0.5 \Bigr]. \]
2. Estimate initial means/covariances within those two groups.
3. Provide this classification to `Mclust` via `initialization=list(classification=…)`.

This preserves the information already extracted while allowing the multivariate model to refine it.

---

## 3  Pathogens to model

| Pathogen | Antigen set |
|----------|-------------|
| EBV      | `ebv_vca`, `ebv_ebna1`, `ebv_zebra`, `ebv_ead` |
| CMV      | `cmv_pp150`, `cmv_pp52`, `cmv_pp28` |
| HHV-6    | `hhv6_ie1a`, `hhv6_ie1b`, `hhv6_p101k` |
| HBV      | `hbv_hbc`, `hbv_hbe` |
| HCV      | `hcv_core`, `hcv_ns3` |
| HPV-16   | `hpv16_l1`, `hpv16_e6`, `hpv16_e7` |
| C. trachomatis | `ct_pgp3`, `ct_mompa`, `ct_mompd`, `ct_tarpf1`, `ct_tarpf2`, `ct_porb` |
| H. pylori | `hp_caga`, `hp_vaca`, `hp_omp`, `hp_groel`, `hp_catalase`, `hp_urea` |

The script will silently skip antigens missing from the phenotype table (e.g. if mixture fitting failed upstream).

---

## 4  Outputs

For each pathogen and each individual:

* `pathogen_sero_soft` – posterior mean probability \(p_{j}^{\text{pathogen}}\)
* `pathogen_sero_hard` – MAP call (`1` if probability ≥ 0.5)
* `pathogen_w_soft`    – weight = `p_sero_soft`
* `pathogen_w_hard`    – weight = `2·|p_sero_soft − 0.5|`

Files written (mirroring Script 02 conventions):

```
quickdraws_input/phenotypes_pathogen.tsv            # FID/IID + soft & hard calls
quickdraws_input/<pathogen>_p_soft.weights          # per-sample soft weights
quickdraws_input/<pathogen>_sero_hard.weights       # per-sample hard weights
pathogen_mixture_model_fits.rds                    # list of fitted Mclust objects
```

No manifest editing is performed automatically; you can add the new phenotypes and weight files to downstream analyses as desired.

---

## 5  Software & resource notes

* R ≥ 4.0, packages: **`tidyverse`, `mclust`, `furrr`, `glue`, `data.table`** (all installed via `pacman` like Script 02).
* Parallelisation: `future::plan(multisession)` with `availableCores() − 1` workers.
* Memory: each mixture fit involves at most a few thousand rows × up to 7 dimensions → negligible on RAP.

---

## 6  Next steps

1. Commit this design document.
2. Implement **`04_pathogen_level_phenotypes.R`** following the above specification.
3. Re-run on RAP; inspect convergence diagnostics and variance explained. 