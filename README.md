# netsim

An R package for simulating city-scale contact networks calibrated to egocentric survey data, using Stochastic Approximation of Networks (SAN) from the `ergm.ego` framework.

## Overview

`netsim` constructs synthetic population contact networks that match empirical targets derived from egocentric network surveys. Given a set of sampled egos (with demographic attributes) and population-level mixing statistics (mixing matrices and edgewise shared partners by contact mode), it:

1. Builds a basis network matching the target attribute mixing matrix
2. Rewires edges via SAN to match triangulation targets (edgewise shared partners)
3. Outputs a full `network` object for downstream epidemic simulation

Contact modes (household, work, school, other) are handled independently and then combined into a single multiplex network.

## Scripts

| File | Description |
|------|-------------|
| `city_sim_singlescript.R` | Self-contained implementation of `citysim()` — the main simulation function |
| `R/citysim.R` | Package version of `citysim()` |
| `R/basis_network_create.R` | `CreateNetworkBasis()` — builds initial network matching mixing matrix targets |
| `R/get_basis_inputs.R` | Preprocesses `eldemog2.Rdata` travel survey data into `basis_inputs.Rdata` |
| `R/find_transitivity_value.R` | Helper to calibrate transitivity targets |
| `R/modify_mm.R` | Utilities for adjusting mixing matrices |
| `R/SIR.R` | Simple SIR epidemic simulation on a contact network |
| `R/example.R` | Usage examples: `citysim()` from an egor object and from precomputed sufficient statistics |

## Data

| File | Description | Size |
|------|-------------|------|
| `eldemog2.Rdata` | Travel survey contact data (raw input) | 254 MB — not included |
| `prem_mixing_matrices.xlsx` | PREM age-structured contact matrices by setting | external |
| `basis_inputs.Rdata` | Preprocessed sufficient statistics from `get_basis_inputs.R` | regenerated |

`eldemog2.Rdata` is excluded from this repository. Run `R/get_basis_inputs.R` (pointing `path` to your copy of the PREM mixing matrices) to regenerate `basis_inputs.Rdata`.

## Usage

```r
# Install the package
devtools::install(".")
# or load in-place
devtools::load_all(".")

# Use case 1: from a fully observed egocentric network
library(ergm.ego)
data(faux.mesa.high)
egonet <- as.egor(faux.mesa.high)
net <- citysim(egonet = egonet, N_pop = 10000, xnames = c("Grade", "Sex"))

# Use case 2: from precomputed sufficient statistics
load("basis_inputs.Rdata")   # -> sampled_nodes, samstats
nets <- list()
for (mode in names(samstats)) {
  nets[[mode]] <- citysim(egos = sampled_nodes, samstats = samstats[[mode]], N_pop = 100000)
}
```

## Requirements

- R (≥ 4.0)
- Packages: `ergm.ego`, `ergm`, `intergraph`, `igraph`, `data.table`, `readxl`
