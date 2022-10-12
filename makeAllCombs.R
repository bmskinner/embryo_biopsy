# Make and save all combinations of embryos of different sizes

library(tessera)
library(parallel)

# Global variables
ANEUPLOIDY.RANGE <- seq(0, 1, 0.01)
DISPERSAL.RANGE <- seq(0, 1, 0.01)
N.REPLICATES <- 100
BIOPSY.SIZES <- c(3:10, 15, 20, 25, 30)
EMBRYO.SIZES <- c(100, 150, 200, 250)
N.CORES <- 40

combs <- expand.grid(
  b = BIOPSY.SIZES,
  s = 1:N.REPLICATES
)

# Loop with regular file write to avoid using too much memory
for (embryo.size in EMBRYO.SIZES) {
  for (a in ANEUPLOIDY.RANGE) {
    for (d in DISPERSAL.RANGE) {
      raw.part.file <- paste0("data/raw/raw_values_e", embryo.size, "_a", a, "_d", d, ".csv")

      if (!file.exists(raw.part.file)) {
        # cat("Aneuploidy:", a, "\tDispersal:", d, "\n")
        get.biopsies <- function(b, s) {
          e <- tessera::Embryo(
            n.cells = embryo.size,
            prop.aneuploid = a,
            dispersal = d,
            concordance = 1, rng.seed = s
          )
          o <- tessera::takeAllBiopsies(e, biopsy.size = b, chromosome = 1)

          c(a, d, b, embryo.size, s, o)
        }

        r <- do.call(rbind, mcmapply(get.biopsies, b = combs$b, s = combs$s, SIMPLIFY = F, mc.cores = N.CORES))
        colnames(r) <- c("Aneuploidy", "Dispersal", "Biopsy_size", "Embryo_size", "Seed", paste0("V", 1:N.CORES))

        write.table(r, raw.part.file, append = F, col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  }
}
