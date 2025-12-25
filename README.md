# PolynomialChaos

MATLAB experiments for building polynomial chaos surrogates of lognormal random processes and variables. The scripts walk through Karhunen–Loève expansions (KLE), Monte Carlo sampling, and both intrusive and non-intrusive PCE workflows for 1-D, 2-D, and 3-D domains with exponential covariance kernels.

<img width="770" height="840" alt="Screenshot 2025-12-25 at 5 45 18 PM" src="https://github.com/user-attachments/assets/396fb329-e54d-4760-8ea6-c3b9b67c36f8" />


## Repository Structure

| Path | Highlights |
| --- | --- |
| `KLE/` | 1-D exponential covariance eigen-analysis. `ExpCov_Rog_diffB.m` finds eigenvalues/eigenfunctions and writes `g_x.dat`, while helper scripts such as `graphical_roots.m`, `kle_omg_lam.m`, and the `*_check` files verify spectra and orthogonality. |
| `Lognormal/` | Lognormal random variable and field studies. Includes non-intrusive reconstructions (`LNProcess_NISP.m`, `PCE_1Dmesh.m`, `PCE_2var_1dmesh.m`), intrusive/Galerkin derivations (`PCE_LNRV_Intrusive_MCS.m`, `PCE_LNRV_Galerkin_MCS.m`), and Monte Carlo validation scripts (`Log_Exponential_MCS.m`, `PCE_LNRV1D_MCS.m`, `PCE_LNRV2D_MCS.m`). |
| `2DLNProcess/` | Two-dimensional KLE and lognormal field generation. `Exponential_TwoDimension.m` builds tensor-product eigen-pairs and writes `eigfun_2D.dat`; `LNProcess_2D_NISP.m` consumes that file to run a non-intrusive stochastic projection. |
| `3Deigen_wave/` | Prototype for 3-D fields. `lambda.m` combines 1-D spectra (see `lambda_114.mat`, `lambda_200.mat`, etc.) into multi-index products and plots decay/energy of the resulting eigenvalues. |
| `misc/` | Small supporting experiments comparing Monte Carlo and kernel density estimates (`MCSvsKDE.m`), computing scalar PCE coefficients (`PCEwithMCS.m`, `check_pcvariation.m`), and storing test data (`MCS_pcvariation.mat`). |

## Requirements

- MATLAB R2018b or later is recommended; Symbolic Math Toolbox is required for scripts that declare `syms` (e.g., `PCE_LNRV_Galerkin_MCS.m`).
- The Statistics and Machine Learning Toolbox is used for `ksdensity`, `lognpdf`, and related routines.
- Scripts read/write plain-text `.dat` files (`g_x.dat`, `eigfun_2D.dat`) and `.mat` files. Keep them in the same folder as the corresponding `.m` files to avoid path issues.

## Typical Workflows

- **1-D lognormal field via KLE → NISP (`KLE/` + `Lognormal/`).** Run `KLE/ExpCov_Rog_diffB.m` after setting the correlation length (`b`), KLE truncation, and grid. The script saves `g_x.dat`, which `Lognormal/LNProcess_NISP.m` then loads to project a Monte Carlo ensemble onto Hermite polynomials and compare marginal PDFs.
- **2-D tensor-product field (`2DLNProcess/`).** Execute `Exponential_TwoDimension.m` to compute sorted tensor-product eigenvalues (`lambda_multi`) and to write `eigfun_2D.dat`. Afterwards, `LNProcess_2D_NISP.m` rebuilds the truncated lognormal expansion, computes statistics from 50k samples, and plots the resulting PCE coefficients.
- **3-D eigenvalue bookkeeping (`3Deigen_wave/`).** Use `lambda.m` to combine `lambda_114.mat`/`lambda_200.mat` spectra, inspect cumulative energy content, and visualize the eigenfunction normalization factors stored in `mult_*.mat`. Adjust the `omega_*` vectors inside the script if you generate alternative correlation lengths.
- **Standalone PCE examples (`Lognormal/`, `misc/`).** Scripts such as `PCE_2var_1dmesh.m` and `PCE_LNRV_Intrusive_MCS.m` show how to compute coefficients analytically or via Monte Carlo for low-dimensional systems. `misc/MCSvsKDE.m` demonstrates validating Monte Carlo histograms with analytical PDFs and kernel density estimates.

## Extending or Adapting

1. Tune correlation length (`b`), variance (`sigma_g`), and truncation order (`KLE_dim`, `ord_in`) in the KLE scripts to match your input random field.
2. Regenerate `g_x.dat` or `eigfun_2D.dat` whenever the stochastic basis changes so that downstream lognormal/PCE scripts pick up the new eigenfunctions.
3. To switch between intrusive vs non-intrusive formulations, reuse the Hermite basis definitions shown in the `Lognormal/` scripts—only the coefficient estimation block needs to change.
4. Use the plotting sections (e.g., cumulative energy plots) as quick diagnostics for how many modes are needed to capture the variance you care about before running expensive Monte Carlo studies.

## Reference

> **[Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems](https://doi.org/10.22215/etd/2023-15817)**

**Authors:** Sudhi Sharma Padillath Vasudevan
**Institution:** Carleton University (2023)  
**DOI:** [10.22215/etd/2023-15817](https://doi.org/10.22215/etd/2023-15817)

<details>
<summary><b>Click to expand BibTeX citation</b></summary>

```bibtex
@phdthesis{vasudevan2023scalable,
  title={Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems},
  author={Vasudevan, Padillath and Sharma, Sudhi},
  year={2023},
  school={Carleton University},
  doi={10.22215/etd/2023-15817}
}
\```
</details>

## Questions ?

Contact : Sudhi Sharma P V  
Email: sudhisharmapadillath@gmail.com
