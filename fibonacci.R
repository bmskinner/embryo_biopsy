# Fibonacci lattice model
# library(tidyverse)
# library(plotly)
# library(assertthat)

# Create a sphere of evenly spaced points 
# n.points - the number of points to create
create.blank.sphere = function(n.points){
  # Make a sphere of evenly spaced points using the Fibonacci spiral
  indices = seq(0, n.points-1, 1)+0.5
  phi = acos(pmin(pmax( 1-2*indices/n.points,-1.0),1.0)) # constrain to avoid rounding errors
  theta = pi * (1 + sqrt(5)) * indices
  
  x = cos(theta)*sin(phi)
  y = sin(theta)*sin(phi)
  z = cos(phi)
  d = as.data.frame(cbind(x, y, z))
  
  # Create distance matrix for each point
  # Set the 6 closest points to be neighbours
  for(i in 1:nrow(d)){
    # current row
    dist = sqrt( (d$x-x[i])**2 + (d$y-y[i])**2 + (d$z-z[i])**2) # distance between points
    d[[paste0("d", i)]] = dist
    # calculate distance from this row to current row
    d[[paste0("isNeighbour", i)]] = dist > 0 & dist <= max(head(sort(dist), n=7))
  }
  d$isSeed = FALSE
  
  # assertthat::assert_that(nrow(d)==n.points, msg=paste("Expected", n.points, "cells, found", nrow(d)))
  return(d)
}

# Test if any of the neighbouring cells have an aneuploidy
# d - the blastocyst
# index - the cell to test
# Returns true if any of the closest cells are aneuploid
has.adjacent.aneuploid = function(d, index){
  adj.list = d[[paste0("isNeighbour", index)]]
  return(any(adj.list & d$isSeed))
}

# Make a blastocyst matrix
# n.cells - the number of cells in the blastocyst
# prop.aneuploid - the proportion of aneuploid cells (0-1)
# dispersion - the dispersion of the aneuploid cells (0-1)
# Returns a blastocyst matrix
create.blastocyst = function(n.cells, prop.aneuploid, dispersion){

  d = create.blank.sphere(n.cells)
  
  if(prop.aneuploid==0) return(d)
  
  n.aneuploid = ceiling(max(1, n.cells * prop.aneuploid))
  
  # Create seeds for aneuploid regions
  n.seeds = ceiling(max(1, n.aneuploid * dispersion))
  n.to.make = n.seeds
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isSeed[seed]) next
    # if(has.adjacent.aneuploid(d, seed)) next # spread seeds out (fails when all cells have at least one neighbour)
    d$isSeed[seed] = T
    n.to.make = n.to.make-1L
  }
  # assertthat::assert_that(sum(d$isSeed)==n.seeds, 
                          # msg = paste("Expected", n.seeds, "seeds, found", sum(d$isSeed)))
  # Grow the seeds into neighbouring cells
  n.expands = n.aneuploid - n.seeds
  n.to.make = n.expands
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isSeed[seed]) next # skip cells already aneuploid
    if(!has.adjacent.aneuploid(d, seed)) next # only grow into areas with an existing aneuploidy
    d$isSeed[seed] = T
    n.to.make = n.to.make-1
  }
  # assertthat::assert_that(sum(d$isSeed)==n.aneuploid, 
  #                         msg = paste("Expected", n.aneuploid, "aneuploids, found", sum(d$isSeed)))
  return(d)
}

# Take a sample from a blastocyst. The cell at the given index is taken, 
# plus the closest n neighbouring cells where n = n.sampled.cells-1.
# d - the blastocyst matrix
# n.sampled.cells - the number of cells to biopsy
# seed.sample - the index of the cell to begin biopsying
sample.blastocyst = function(d, n.sampled.cells, seed.sample){
  
  sample.list = d[[paste0("d", seed.sample)]]
  
  d$isSampled = d[[paste0("d", seed.sample)]] <= max(head(sort(sample.list), n=n.sampled.cells))
  return(d)
}

# Find the number of aneuploid cells in all biopsies of 
# a given size for the given blastocyst
# d - the blastocyst matrix
# n.cells.per.sample - the number of cells to take in each biopsy
make.samples = function(d, n.cells.per.sample){
  result = c()
  for(i in 1:nrow(d)){ # sample each cell in turn, so we get every cell
    f = sample.blastocyst(d, n.cells.per.sample, i)
    result = c(result, sum(f[f$isSampled,]$isSeed))
  }
  return(result)
}


