# Bayesian Econometric Model: Predicting S&P 500 Returns with Macroeconomic Indicators

## Overview

This project investigates whether **Bayesian regression** improves the
prediction of monthly S&P 500 equity returns compared to classical
Ordinary Least Squares (OLS) regression, using macroeconomic indicators
as predictors. It also quantifies the uncertainty in macro-return
relationships through posterior distributions and credible intervals.

---

## Research Questions

- Can Bayesian regression improve prediction of equity market returns
  vs. traditional OLS?
- Do macroeconomic indicators have probabilistic predictive power for
  S&P 500 returns?
- How uncertain are the relationships between macro variables and
  equity returns?

---

## Project Structure
```
bayesian-macro-sp500/
│
├── bayesian_macro_sp500.R        # Full analysis script
├── bayesian_macro_sp500.Rmd      # R Markdown report (knittable)
├── README.md                     # This file
│
├── outputs/
│   ├── correlation_heatmap.png
│   ├── posterior_distributions.png
│   ├── posterior_predictive_check.png
│   ├── forecast_comparison.png
│   ├── coefficient_comparison.png
│   ├── model_performance.csv
│   └── forecast_results.csv
```

---

## Data Sources

All data is pulled automatically via the `fredr` package from the
[Federal Reserve Bank of St. Louis (FRED)](https://fred.stlouisfed.org/).

| Variable | FRED Series ID | Description |
|---|---|---|
| S&P 500 Index | `SP500` | Monthly closing price |
| Inflation | `CPIAUCSL` | CPI, All Urban Consumers (YoY % change) |
| Unemployment | `UNRATE` | Civilian Unemployment Rate (%) |
| Interest Rate | `FEDFUNDS` | Federal Funds Effective Rate (%) |
| Industrial Production | `INDPRO` | Industrial Production Index (YoY % change) |

**Sample period:** January 2000 – December 2024 (monthly frequency)

---

## Methodology

### Feature Engineering
- S&P 500 monthly returns computed as: `(P_t / P_{t-1} - 1) * 100`
- Inflation computed as trailing 12-month CPI change
- Industrial production computed as trailing 12-month change
- Train / test split: 80% training, 20% holdout

### Model 1 — OLS Regression
Classical linear regression of S&P 500 returns on the four macro
predictors. Provides point estimates, standard errors, and 95%
prediction intervals.

### Model 2A — Bayesian Regression (rstanarm)
Bayesian linear regression estimated via MCMC (Stan backend).

**Priors:**
- Coefficients: `Normal(0, 2.5)` — weakly informative, autoscaled
- Intercept: `Normal(0, 10)`

**MCMC settings:**
- 4 chains, 4000 iterations, 1000 warmup draws
- Convergence checked via R-hat (< 1.01) and effective sample size

### Model 2B — Bayesian Regression (brms)
More flexible Bayesian regression with explicit prior specification
via the `brms` interface.

**Priors:**
- Coefficients: `Normal(0, 2.5)`
- Intercept: `Normal(0, 10)`
- Sigma: `Exponential(1)`

---

## Key Outputs

| Output | Description |
|---|---|
| Posterior distributions | Full probability distributions over each macro coefficient |
| 95% Credible intervals | Bayesian uncertainty bounds on coefficients |
| Posterior predictive check | Simulated vs actual return distributions |
| MCMC trace plots | Chain mixing and convergence diagnostics |
| Forecast comparison plot | OLS vs Bayesian predictions with uncertainty bands |
| Coefficient comparison plot | Side-by-side OLS CIs vs Bayesian credible intervals |
| RMSE / MAE table | Out-of-sample point forecast performance |
| LOO-IC | Leave-one-out cross-validation score for Bayesian model |

---

## Requirements

### R Version
R >= 4.2.0 recommended

### Packages
```r
install.packages(c(
  "tidyverse",
  "fredr",
  "rstanarm",
  "brms",
  "bayesplot",
  "loo",
  "lubridate",
  "forecast",
  "Metrics",
  "patchwork",
  "scales",
  "broom",
  "knitr",
  "kableExtra"
))
```

### Stan / C++ Toolchain
`rstanarm` and `brms` both require a working C++ compiler via RStan.

**Windows:** Install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)

**Mac:** Install Xcode Command Line Tools:
```bash
xcode-select --install
```

**Linux:** Install `r-base-dev` via your package manager.

Full setup guide: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

---

## Setup & Usage

### 1. Get a FRED API Key
Register for free at:
https://fred.stlouisfed.org/docs/api/api_key.html

### 2. Add your API key
In both `bayesian_macro_sp500.R` and `bayesian_macro_sp500.Rmd`,
replace:
```r
fredr_set_key("YOUR_FRED_API_KEY_HERE")
```
with your actual key:
```r
fredr_set_key("abcd1234yourkey")
```

### 3. Run the script
Open `bayesian_macro_sp500.R` in RStudio and run it top to bottom,
or source it:
```r
source("bayesian_macro_sp500.R")
```

### 4. Knit the report
Open `bayesian_macro_sp500.Rmd` in RStudio and click **Knit**, or:
```r
rmarkdown::render("bayesian_macro_sp500.Rmd")
```

---

## Runtime

| Step | Approximate Time |
|---|---|
| Data pull (FRED) | < 1 minute |
| OLS model | < 1 second |
| rstanarm (4 chains × 4000 iter) | 2 – 5 minutes |
| brms (4 chains × 4000 iter) | 5 – 15 minutes |

Times vary by machine. Using `cores = parallel::detectCores() - 1`
in `brms` will parallelize chains and reduce wall time significantly.

---

## Known Issues & Troubleshooting

**`object 'value' not found` in pivot_longer**
Namespace conflict between `rstanarm`/`rstan` and `tidyr`. Use
explicit namespacing:
```r
tidyr::pivot_longer(cols = dplyr::everything(), ...)
```
Or use the base R `lapply` alternative provided in the script.

**`object 'p_coefs' not found`**
The coefficient comparison block failed silently upstream. Run each
block individually to find where it breaks — most likely the
`pivot_longer` conflict above.

**Stan compilation errors on first run**
Ensure your C++ toolchain is correctly installed. Run
`pkgbuild::has_build_tools(debug = TRUE)` to verify.

**FRED series returns no data**
The `SP500` series on FRED requires a free API key and may have
a short lag. Verify your key is active at
https://fred.stlouisfed.org/docs/api/api_key.html

---

## Interpretation Notes

- **Wider credible intervals vs OLS** are expected and desirable —
  they reflect genuine uncertainty the frequentist intervals
  understate.
- **R-hat < 1.01** on all parameters indicates good MCMC convergence.
- **LOO-IC** is preferred over AIC/BIC for Bayesian model comparison.
- Point forecast RMSE/MAE may be similar between OLS and Bayesian —
  the Bayesian advantage is in uncertainty quantification, not
  necessarily raw accuracy.

---

## License

MIT License. Free to use, modify, and distribute with attribution.

---

## Author

Sara El Fellah  
