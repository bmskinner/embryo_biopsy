# Embryo mosaicism analysis scripts

These scripts generate the model embryos and biopsies described in Skinner _et al_. Explaining the counter-intuitive effectiveness of trophectoderm biopsy for PGT-A using computational modelling.

They are designed to be run on a computer with multiple cores. Set the number of
cores using `parameters.R`. Note that Windows can only use a single core.

## Prerequisites

The scripts use the following R packages from CRAN:

```
fs
parallel
tidyverse
patchwork
data.table
svglite
```

We also use `tessera`, installable from Github:

```
library(devtools)
devtools::install_github("bmskinner/tessera")
```


## Running the scripts

The main scripts:

- `makeAllCombinations.R` : creates the model embryos and takes biopsies. Biopsies are saved in `./data/raw/`
- `analyseCombos.R` : aggregates the raw values for figure generation. Outputs are saved in `./data/aggregates/`
- `makeRankOrderData.R`: creates pools of embryos for testing rank order. Outputs are saved in `./data/`
- `makeTwoBiopsyData.R`: creates embryos and takes two biopsies. Outputs are saved in `./data/two_biopsy/`
- `Figure_xxxx.R` : creates the figures from the paper. Outputs are saved in `./figure/`

Supplementary scripts:

- `functions.R`: common functions used throughout
- `parameters.R`: global variables used throughout