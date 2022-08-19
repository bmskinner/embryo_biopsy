# Create two cell biopsies

library(tessera)
library(parallel)
library(tidyverse)

ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.01)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)
EMBRYO.SIZES = seq(100, 250, 50)
N.CORES = ifelse(Sys.info()["sysname"]=="Windows", 1, 5) 

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

  part.file = paste0("data/two_biopsy/two_biopsy_values_a", aneuploidy, "_d", dispersal, ".csv")
  
  if(!file.exists(part.file)){

    get.biopsies = function(s, embryo.size, biopsy.size){
      e = tessera::Embryo(n.cells = embryo.size,
                          prop.aneuploid = aneuploidy,
                          dispersal = dispersal,
                          concordance = 1, rng.seed = s)
      b = take.two.biopsies(e, biopsy.size = biopsy.size, chromosome = 1)
      
      list(Aneuploidy = aneuploidy, Dispersal = dispersal, 
           Embryo_size = embryo.size, Biopsy_size = biopsy.size, Seed = s, n_aneuploid = b)
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
      tidyr::unnest_longer(., col=n_aneuploid, simplify = T) %>%
      tidyr::unnest_longer(., col=Aneuploidy:Seed) %>%
      dplyr::group_by(Aneuploidy, Dispersal, Embryo_size, Biopsy_size, n_aneuploid) %>%
      dplyr::summarise(Total_biopsies = dplyr::n())

    write.table(r, part.file, append = F, col.names = T, row.names = F, sep = "\t", quote = F)
    rm(r)
    gc()
  }
}


combinations = expand.grid(aneuploidies = ANEUPLOIDY.RANGE,
                           dispersals   = DISPERSAL.RANGE)

# Make the combinations
mapply(make.two.biopsy.data, 
       aneuploidy = combinations$aneuploidies,
       dispersal  = combinations$dispersals)