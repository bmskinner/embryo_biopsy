# Fibonacci lattice model
# library(assertthat)

# The number of cells to be considered neighbours to a given cell
N_NEIGHBOURS=6

# Create a sphere of evenly spaced points
# n.points - the number of points to create
create.blank.sphere = function(n.points){
  # Make a sphere of evenly spaced points using the Fibonacci spiral
  indices = seq(0, n.points-1, 1)+0.5
  # Map to a sphere via cylindrical equal area projection
  phi = acos(pmin(pmax( 1-2*indices/n.points,-1.0),1.0)) # constrain to avoid rounding errors
  theta = pi * (1 + sqrt(5)) * indices

  # Transform spherical coordinates to Cartesian
  x = cos(theta)*sin(phi)
  y = sin(theta)*sin(phi)
  z = cos(phi)
  d = as.data.frame(cbind(x, y, z))

  # Create distance matrix for each point
  # Set the N_NEIGHBOURS closest points to be neighbours
  for(i in 1:nrow(d)){
    dist = sqrt( (d$x-x[i])**2 + (d$y-y[i])**2 + (d$z-z[i])**2) # distance between points
    d[[paste0("d", i)]] = dist # create a column to store the distances
    # A point is a neighbour if it is not this point, and it is in the list of closest points
    d[[paste0("isNeighbour", i)]] = dist > 0 & dist <= max(head(sort(dist), n=N_NEIGHBOURS+1))
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

# Count the number of potential seed locations
# (cells that are not adjacent to a seed)
# d - the blastocyst matrix
# Returns the number of potential seeds
count.empty.blocks = function(d){
  n = 0
  for(i in 1:nrow(d)){
    if(!d$isSeed[i] & !has.adjacent.aneuploid(d, i)) n = n+1
  }
  return(n)
}

# Make a blastocyst matrix
# n.cells - the number of cells in the blastocyst
# prop.aneuploid - the proportion of aneuploid cells (0-1)
# dispersion - the dispersion of the aneuploid cells (0-1)
# Returns a blastocyst matrix
create.blastocyst = function(n.cells, prop.aneuploid, dispersion){

  d = create.blank.sphere(n.cells)

  # Shortcut the easy cases
  if(prop.aneuploid==0) return(d)

  if(prop.aneuploid==1){
    d$isSeed=T
    return(d)
  }

  # We must have an integer value of at least one aneuploid cell
  n.aneuploid = ceiling(max(1, n.cells * prop.aneuploid))

  # The approach for dispersal is to set seed cells which will
  # grow into separate aneuploid patches. The more dispersion, the more
  # initial seeds.

  # Choose number of seeds for aneuploid regions
  n.seeds = ceiling(max(1, n.aneuploid * dispersion))
  n.to.make = n.seeds

  # We can disperse up to a certain number of initial blocks with
  # no aneuploid neighbours. After this, every cell will have at least
  # one aneuploid neighbour. We stop a bit before this to make the maths simpler.
  initial.blocks = max(1,floor(n.cells/N_NEIGHBOURS))

  # Disperse seeds as much as possible
  while(initial.blocks>0 & n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isSeed[seed]) next
    if(has.adjacent.aneuploid(d, seed)) next # spread seeds out
    d$isSeed[seed] = T
    n.to.make = n.to.make-1L
    initial.blocks = initial.blocks-1L
  }

  # When all dispersed seeds have been added, add the remaining seeds randomly
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isSeed[seed]) next
    d$isSeed[seed] = T
    n.to.make = n.to.make-1L
  }
  # assertthat::assert_that(sum(d$isSeed)==n.seeds,
                          # msg = paste("Expected", n.seeds, "seeds, found", sum(d$isSeed)))

  # Grow the seeds into neighbouring cells for remaining aneuploid cells
  n.to.make = n.aneuploid - n.seeds
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isSeed[seed]) next # skip cells already aneuploid
    if(!has.adjacent.aneuploid(d, seed)) next # only grow next to existing aneuploid
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


