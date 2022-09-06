# Embryo mosaicism analysis scripts

These scripts generate the model embryos and biopsies described in xxxx.

They are designed to be run on a computer with multiple cores. Set the number of
cores using `parameters.R`. Note that Windows can only use a single core.

There are three main scripts:

- makeCombos.R : creates the model embryos and takes biopsies. Biopsies are saved in ./data/raw/
- analyseCombos.R : aggregates the raw values for figure generation. Outputs are saved in ./data/aggregates/
- makeFigures.R : creates summary figures. Outputs are saved in ./figure/