# Explore the blastocyt data
library(parallel)
source("fibonacci.R")

conditions = expand.grid("disps"     = seq(0, 1, 0.2), 
                         "aneuploids"= seq(0, 1, 0.2), 
                         "samples"   = seq(1, 20, 1), 
                         "cells"     = seq(50, 200, 50), 
                         "reps"      = seq(1, 50, 1))

# Create a blastocyst with the given characteristics and sample it
run.simulation = function(disp, aneuploid, cells, replicate, cells.per.sample){
  # cat("Simulating", disp, "dispersion,", cells, "cells,", aneuploid, "aneuploids,", cells.per.sample, "cells per sample, rep",replicate,"\n")
  d = create.blastocyst(n.cells=cells, prop.aneuploid = aneuploid, dispersion = disp)
  result = make.samples(d, n.cells.per.sample=cells.per.sample)
  return(result)
}

set.seed(42)
conditions$output = mcmapply(run.simulation,
                           disp            = conditions$disps,
                           aneuploid       = conditions$aneuploids,
                           cells           = conditions$cells,
                           cells.per.sample= conditions$samples,
                           replicate       = conditions$reps,
                           SIMPLIFY = F,
                           mc.cores = 60)

# Convert list to charactar for export
conditions$output = vapply(conditions$output, paste, collapse = ", ", character(1L))

write.table(conditions, file="conditions.tsv", quote = T, col.names = T, row.names = F, sep = "\t")

