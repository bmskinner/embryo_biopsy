# Explore the blastocyt data
library(parallel)
source("fibonacci.R")

# Create the possible conditions to test
conditions = expand.grid("disps"     = seq(0, 1, 0.05), 
                         "aneuploids"= seq(0, 1, 0.05), 
                         "samples"   = seq(1, 20, 1), 
                         "cells"     = seq(50, 300, 50), 
                         "reps"      = seq(1, 100, 1))

# Create a blastocyst with the given characteristics and sample it
run.simulation = function(disp, aneuploid, cells, replicate, cells.per.sample, index, total){
  if(index%%20000==0){
    cat(format(Sys.time(), "%Y-%m-%d %X"), ": Simulation",index, "of", total, "(", format(index/total, nsmall=1, digits = 3), "%)\n")
  }

  d = create.blastocyst(n.cells=cells, prop.aneuploid = aneuploid, dispersion = disp)
  result = make.samples(d, n.cells.per.sample=cells.per.sample)
  return(result)
}

n.cores = ifelse(Sys.info()["sysname"]=="Windows", 1, 60) # parallel only works on Unix-like
set.seed(42)

cat(format(Sys.time(), "%Y-%m-%d %X"), ": Running ", nrow(conditions), "simulations on ", n.cores, "cores\n")
conditions$output = mcmapply(run.simulation,
                           disp            = conditions$disps,
                           aneuploid       = conditions$aneuploids,
                           cells           = conditions$cells,
                           cells.per.sample= conditions$samples,
                           replicate       = conditions$reps,
                           index           = seq(1:nrow(conditions)),
                           total           = nrow(conditions),
                           SIMPLIFY = F,
                           mc.cores = n.cores)

# Convert list to character for export. vapply can be faster than sapply if output format can
# be pre-specified
conditions$output = vapply(conditions$output, paste, collapse = ", ", character(1L))

# Save the results to file
write.table(conditions, file="conditions.tsv", quote = T, col.names = T, row.names = F, sep = "\t")

