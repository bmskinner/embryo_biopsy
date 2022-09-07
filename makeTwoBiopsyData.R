# Create two cell biopsies

library(tessera)
library(parallel)
library(tidyverse)
library(data.table)
library(fs)

source("parameters.R")
source("functions.R")

# For a given embryo, take all possible two biopsies
take.two.biopsies = function(embryo, biopsy.size, chromosome){
  
  # cat("Taking two biopsies\n")
  getNeighboursToExclude = function(embryo, cell.index){
    # get the direct neighbours - these will be in the biopsy
    neighbours = tessera::getNeighbouringCellIndexes(embryo, cell.index)
    
    # Get neighbours of neighbours so we don't have any surrounding biopsies 
    # including the missing cells
    neighbours2 = sapply(neighbours, tessera::getNeighbouringCellIndexes, embryo=embryo)
    return(unique(c(cell.index, neighbours, neighbours2)))
  }
  
  all.cells = 1:length(embryo)
  
  # Get all biopsies in a vector indexed by cell
  all.biopsies = sapply(all.cells, function(i){
    tessera::takeBiopsy(embryo, biopsy.size = biopsy.size, chromosome = chromosome, index.cell = i)
  })
  
  # Now find the biopsy combinations to look up
  unlist(sapply(1:length(embryo), function(i){
    # We now need to exclude neighbouring cells for second biopsy
    # Cells that are not neighbours or neighbours of neighbours
    keep.list = all.cells[!all.cells %in% getNeighboursToExclude(embryo, i)]
    
    sapply(keep.list, function(j){ (all.biopsies[i] + all.biopsies[j])/2 })
  }))
}

# For a combination of aneuploidy and dispersal, calculate all two biopsy values
make.two.biopsy.data = function(aneuploidy, dispersal){
  
  part.file = paste0(TWO.BIOPSY.DATA.PATH, "/two_biopsy_values_a", aneuploidy, "_d", dispersal, ".csv")
  
  if(!file.exists(part.file)){
    
    get.biopsies = function(s, embryo.size, biopsy.size){
      e = tessera::Embryo(n.cells = embryo.size,
                          prop.aneuploid = aneuploidy,
                          dispersal = dispersal,
                          concordance = 1, rng.seed = s)
      b = take.two.biopsies(e, biopsy.size = biopsy.size, chromosome = 1)
      
      list(Aneuploidy = aneuploidy, Dispersal = dispersal, 
           Embryo_size = embryo.size, Biopsy_size = biopsy.size, n_aneuploid = b)
    }
    
    ad.combo = expand.grid(embryo.size  = EMBRYO.SIZES,
                           biopsy.size  = BIOPSY.SIZES,
                           seed         = 1:N.REPLICATES)
    
    r = do.call(rbind, mcmapply(get.biopsies, 
                                s           = ad.combo$seed, 
                                embryo.size = ad.combo$embryo.size,
                                biopsy.size = ad.combo$biopsy.size,
                                mc.cores = N.CORES, SIMPLIFY = F)) %>% 
      as.data.frame %>%
      tidyr::unnest_longer(., col=n_aneuploid, simplify = T)
      
    r$Aneuploidy  = unlist(r$Aneuploidy)
    r$Dispersal   = unlist(r$Dispersal)
    r$Embryo_size = unlist(r$Embryo_size)
    r$Biopsy_size = unlist(r$Biopsy_size)
    
    r = r %>% 
      dplyr::group_by(Aneuploidy, Dispersal, Embryo_size, Biopsy_size, n_aneuploid) %>%
      dplyr::summarise(Total_biopsies = dplyr::n())
    
    write.table(r, part.file, append = F, col.names = T, row.names = F, sep = "\t", quote = F)
    rm(r)
    gc()
  }
}

if(!fs::dir_exists(TWO.BIOPSY.DATA.PATH)){
  fs::dir_create(TWO.BIOPSY.DATA.PATH, recursive = T)
}

combinations = expand.grid(aneuploidies = ANEUPLOIDY.RANGE,
                           dispersals   = DISPERSAL.RANGE)

# Make the combinations and generate the data files
# Note - do not parallel here because each make.two.biopsy.data call invokes
# mcmapply
mapply(make.two.biopsy.data, 
       aneuploidy  = combinations$aneuploidies,
       dispersal   = combinations$dispersals)

