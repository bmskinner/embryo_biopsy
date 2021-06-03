# Fibonacci lattice model
library(tidyverse)
library(plotly)
library(assertthat)

# Create a sphere of evenly spaced points 
create.blank.sphere = function(n.points){
  # Make a sphere of evenly spaced points using the Fibonacci spiral
  indices = seq(0, n.points, 1)+0.5
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
  return(d)
}

# Test if any of the neighbouring cells have an aneuploidy
has.adjacent.aneuploid = function(d, index){
  adj.list = d[[paste0("isNeighbour", index)]]
  return(any(adj.list & d$isSeed))
}

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
  assertthat::assert_that(sum(d$isSeed)==n.seeds, 
                          msg = paste("Expected", n.seeds, "seeds, found", sum(d$isSeed)))
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
  assertthat::assert_that(sum(d$isSeed)==n.aneuploid, 
                          msg = paste("Expected", n.aneuploid, "aneuploids, found", sum(d$isSeed)))
  return(d)
}

# Take a sample from a given blastocyst of the given index
sample.blastocyst = function(d, n.sampled.cells, seed.sample){
  
  sample.list = d[[paste0("d", seed.sample)]]
  
  d$isSampled = d[[paste0("d", seed.sample)]] <= max(head(sort(sample.list), n=n.sampled.cells))
  return(d)
}

# Find the number of aneuploid cells found in a sample of 
# a given size for the given blastocyst
make.samples = function(d, n.cells.per.sample){
  result = c()
  for(i in 1:nrow(d)){ # sample each cell in turn, so we get every cell
    f = sample.blastocyst(d, n.cells.per.sample, i)
    result = c(result, sum(f[f$isSampled,]$isSeed))
  }
  return(result)
}


disps = seq(0, 1, 0.1)
cells = seq(50, 200, 10)
aneus = seq(0, 1, 0.1)
conditions = expand.grid("disps"=disps, "aneuploids"=aneus)

# Create a blastocyst with the given characteristics and sample it
run.simulation = function(disp, aneuploid, cells){
  cat("Simulating", disp, "dispersion,", aneuploid, "aneuploids\n")
  d = create.blastocyst(n.cells=cells, prop.aneuploid = aneuploid, dispersion = disp)
  sample.results = make.samples(d, n.cells.per.sample=5)
  return(mean(sample.results))
}


conditions$output = mapply(run.simulation, conditions$disps, conditions$aneuploids, cells=200)
ggplot(conditions, aes(x=disps, y=aneuploids, fill=output))+
  geom_tile()+
  scale_fill_viridis_c()

# Single sample
# d = create.blastocyst(n.cells=200, prop.aneuploid = 0.1, dispersion = 0.3)
# sample.results = make.samples(d, n.samples=2000, n.cells.per.sample=5)


plot_ly(x=d$x, y=d$y, z=d$z, type="scatter3d",
        mode="markers", 
        color=d$isSeed, 
        colors = c("#00FF00", "#FF0000"))

