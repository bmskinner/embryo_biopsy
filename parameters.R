# Shared variables
ANEUPLOIDY.RANGE <- seq(0, 1, 0.01)
DISPERSAL.RANGE <- seq(0, 1, 0.01)
N.REPLICATES <- 100
BIOPSY.SIZES <- c(3:10, 15, 20, 25, 30)
EMBRYO.SIZES <- c(100, 150, 200, 250)
N.CORES <- ifelse(Sys.info()["sysname"] == "Windows", 1, 40)

RAW.DATA.PATH <- "data/raw"
AGGREGATE.DATA.PATH <- "data/aggregates"
TWO.BIOPSY.DATA.PATH <- "data/two_biopsy"

# Where figures should be saved
FIGURE.OUTPUT.DIR <- "figure"
