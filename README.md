# Embryo mosaicism analysis scripts

These scripts generate the model embryos and biopsies described in xxxx.

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

There are three main scripts:

- `makeCombos.R` : creates the model embryos and takes biopsies. Biopsies are saved in `./data/raw/`
- `analyseCombos.R` : aggregates the raw values for figure generation. Outputs are saved in `./data/aggregates/`
- `Figure_xxxx.R` : create the summary figures. Outputs are saved in `./figure/`

Clone this repo and invoke the scripts using Rscript _e.g_:

```
Rscript makeCombos.R
```

Adapt to pass to your computing cluster as appropriate.