# Create two cell biopsies

library(tessera)
library(parallel)
library(tidyverse)
library(data.table)
library(fs)

source("parameters.R")
source("functions.R")

# Create the possible embryo, biopsy and seed sizes
ad.combo = expand.grid(embryo.size  = EMBRYO.SIZES,
                       biopsy.size  = BIOPSY.SIZES,
                       seed         = 1:N.REPLICATES)

# For a given embryo, take all possible two biopsies
take.two.biopsies = function(embryo, biopsy.size, chromosome){
  
  # cat("Taking two biopsies\n")
  getNeighboursToExclude = function(embryo, cell.index){
    # get the indexes of direct neighbours - these will be in the biopsy
    neighbours = tessera::getNeighbouringCellIndexes(embryo, cell.index)
    
    # Get indexes of neighbours of neighbours so we don't have any surrounding biopsies 
    # including the missing cells
    neighbours2 = vapply(neighbours, tessera::getNeighbouringCellIndexes, FUN.VALUE = integer(6), embryo=embryo)
    return(unique(c(cell.index, neighbours, neighbours2)))
  }
  
  all.cells = 1:length(embryo)
  
  # Get all biopsies in a vector indexed by cell
  all.biopsies = vapply(all.cells, function(i){
    tessera::takeBiopsy(embryo, biopsy.size = biopsy.size, chromosome = chromosome, index.cell = i)
  }, FUN.VALUE = integer(1))
  
  # Now find the biopsy combinations to look up
  unlist(sapply(1:length(embryo), function(i){
    # We now need to exclude neighbouring cells for second biopsy
    # Cells that are not neighbours or neighbours of neighbours
    keep.list = all.cells[!all.cells %in% getNeighboursToExclude(embryo, i)]
    
    vapply(keep.list, function(j){ (all.biopsies[i] + all.biopsies[j])/2 }, FUN.VALUE = numeric(1))
  }), use.names = FALSE)
}

# For a combination of aneuploidy and dispersal, calculate all two biopsy values
make.two.biopsy.data = function(aneuploidy, dispersal){
  # Output for the data
  part.file = paste0(TWO.BIOPSY.DATA.PATH, "/two_biopsy_values_a", aneuploidy, "_d", dispersal, ".csv")
  
  if(!file.exists(part.file)){
    # cat("Making for aneuploidy dispersal combo\n")
    # Generate the biopsies from an embryo of the given size
    get.biopsies = function(s, embryo.size, biopsy.size){
      # cat("Getting biopsies\n")
      e = tessera::Embryo(n.cells = embryo.size,
                          prop.aneuploid = aneuploidy,
                          dispersal = dispersal,
                          concordance = 1, rng.seed = s)
      b = take.two.biopsies(e, biopsy.size = biopsy.size, chromosome = 1)
      
      list(Aneuploidy = aneuploidy, Dispersal = dispersal, 
           Embryo_size = embryo.size, Biopsy_size = biopsy.size, n_aneuploid = b)
    }

    # Run all combinations and combine data
    # Note - do not parallel here because we use mcmapply in the main loop
    r = do.call(rbind, mapply(get.biopsies, 
                                s           = ad.combo$seed, 
                                embryo.size = ad.combo$embryo.size,
                                biopsy.size = ad.combo$biopsy.size,
                                SIMPLIFY = F))
    # cat("Made aneuploidy dispersal combo\n")
    
    r = r %>% 
      as.data.frame %>%
      tidyr::unnest_longer(., col=n_aneuploid, simplify = T)
    # cat("Unnested\n")
      
    r$Aneuploidy  = unlist(r$Aneuploidy, use.names = FALSE)
    r$Dispersal   = unlist(r$Dispersal, use.names = FALSE)
    r$Embryo_size = unlist(r$Embryo_size, use.names = FALSE)
    r$Biopsy_size = unlist(r$Biopsy_size, use.names = FALSE)
    
    r = r %>% 
      dplyr::group_by(Aneuploidy, Dispersal, Embryo_size, Biopsy_size, n_aneuploid) %>%
      dplyr::summarise(Total_biopsies = dplyr::n())
    
    # Export data to file and cleanup
    write.table(r, part.file, append = F, col.names = T, row.names = F, sep = "\t", quote = F)
    rm(r)
    gc()
  }
}

# Make output folder if needed
if(!fs::dir_exists(TWO.BIOPSY.DATA.PATH)){
  fs::dir_create(TWO.BIOPSY.DATA.PATH, recursive = T)
}

# Make combinations of parameters to run
combinations = expand.grid(aneuploidies = ANEUPLOIDY.RANGE,
                           dispersals   = DISPERSAL.RANGE)

# Generate the data files
mcmapply(make.two.biopsy.data, 
         aneuploidy  = combinations$aneuploidies,
         dispersal   = combinations$dispersals,
         mc.cores    = N.CORES)


