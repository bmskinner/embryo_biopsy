# Make and save all combinations of embryos

library(tessera)
library(parallel)
library(data.table)

# Global variables
ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.01)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)

out.file = "data/all.combos.csv"

if(!file.exists(out.file)){
  
  # Write the column header
  cat("Aneuploidy\tDispersal\tBiopsy_size\tSeed\t", paste0("V", 1:200, collapse = "\t"), "\n", 
      file = out.file, append = F)
  
  combs = expand.grid(b = BIOPSY.SIZES,
                      s = 1:N.REPLICATES)
  
  # Loop with regular file write to avoid using too much memory
  for(a in ANEUPLOIDY.RANGE){
    
    for(d in DISPERSAL.RANGE){
      cat("Aneuploidy:", a, "\tDispersal:", d, "\n")
      
      get.biopsies = function(b, s){
        e = tessera::Embryo(n.cells = 200,
                            prop.aneuploid = a,
                            dispersal = d,
                            concordance = 1, rng.seed = s)
        o = tessera::takeAllBiopsies(e, biopsy.size = b, chromosome = 1)
        
        c(a, d, b, s, o)
      }
      
      r = do.call(rbind, mcmapply(get.biopsies, combs$b, combs$s, SIMPLIFY = F, mc.cores = 20))
      write.table(r, out.file, append = T, col.names = F, row.names = F, sep = "\t")
    }
  }
}



