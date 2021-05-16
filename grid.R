# Ideas for using a static grid of cells rather than random
library(tidyverse)
library(matrixStats)
library(assertthat)

# Test if the 8-connected block around the point
# has the same value.
# Return true if all points have the same value
block.has.value = function(m, x, y){
  
  value = m[y, x]
  x.left  = ifelse(x==1, ncol(m), x-1)
  x.right = ifelse(x==ncol(m), 1, x+1)
  y.below = ifelse(y==1, nrow(m), y-1)
  y.above = ifelse(y==nrow(m), 1, y+1)
  
  # print(paste("above:", y.above, "below:", y.below, "left:", x.left, "right:", x.right))
  
  return(  m[y.above, x.left ]==value & 
           m[y.above, x      ]==value & 
           m[y.above, x.right]==value &
           m[y,       x.left ]==value & 
           m[y,       x.right]==value &
           m[y.below, x.left ]==value & 
           m[y.below, x      ]==value & 
           m[y.below, x.right]==value)
}

# Get all points in the matrix with the given value
# which are not completely surrounded by points with
# the same value
# Returns a data frame with x and y columns
get.points = function(m, value){
  result = data.frame('x'=NA, 'y'=NA)
  for(x in 1:ncol(m)){
    for(y in 1:nrow(m)){
      # We don't need to include in seeds values that are 
      # completely surrounded by points with the same value
      if(m[y, x]==value & !block.has.value(m, x, y)){
        result = rbind(result, c('x'=x, 'y'=y))
      }
    }
  }
  return(na.omit(result))
}

# Make a seed at a random location in a matrix
# Returns the modified matrix
make.seed = function(m, value){
  x = sample(1:ncol(m), size=1)
  y = sample(1:nrow(m), size=1)
  if(m[y, x]!=value){
    m[y, x] <- value
    return(m)
  } else make.seed(m, value)
}

# Choose a random point that is a seed
# Returns a point
pick.seed = function(m, value){
  all.seeds = get.points(m, value)
  # print("All possible seeds:")
  # print(all.seeds)
  random.row = sample(1:nrow(all.seeds), size=1)
  # print(paste("Row:", random.row))
  x = all.seeds$x[random.row]
  y = all.seeds$y[random.row]
  
  if(nrow(all.seeds)==0) print(m)
  # cat("x", x, "y", y, "\n")
  if(m[y, x]==value) return(list("x" = x, "y"=y))
  return(pick.seed(m, value))
}

# Pick an 8 connected point in a matrix
# that has the given value
# Returns a point
pick.adjacent = function(m, point, value){
  
  # Wrap edges
  x = point[['x']]
  y = point[['y']]
  x.left  = ifelse(x==1, ncol(m), x-1)
  x.right = ifelse(x==ncol(m), 1, x+1)
  y.below = ifelse(y==1, nrow(m), y-1)
  y.above = ifelse(y==nrow(m), 1, y+1)
  
  # Find all valid points
  possible.points = c()
  
  if(m[y.above, x.left ]==value) possible.points = c(possible.points, x.left , y.above)
  if(m[y.above, x      ]==value) possible.points = c(possible.points, x      , y.above)
  if(m[y.above, x.right]==value) possible.points = c(possible.points, x.right, y.above)
  
  if(m[y      , x.left ]==value) possible.points = c(possible.points, x.left , y      )
  if(m[y      , x.right]==value) possible.points = c(possible.points, x.right, y      )
  
  if(m[y.below, x.left ]==value) possible.points = c(possible.points, x.left , y.below)
  if(m[y.below, x      ]==value) possible.points = c(possible.points, x      , y.below)
  if(m[y.below, x.right]==value) possible.points = c(possible.points, x.right, y.below)
  
  # Choose a random point from those that are valid
  new.point = sample(1:length(possible.points)/2, size=1)*2
  
  return(list('x'=possible.points[new.point-1], 'y'=possible.points[new.point]))
}

# Set a value in a matrix
# Returns the modified matrix
set.value = function(m, point, value){
  m[point[['y']], point[['x']]] <- value
  return(m)
}

# Add aneuploid cells to a matrix
# Returns the modified matrix
make.aneuploids = function(cell.matrix, dispersion, n.cells){
  
  w = ncol(cell.matrix)
  h = nrow(cell.matrix)
  total.cells = w * h
  
  # If there are no aneuploid cells, return the original matrix
  if(n.cells == 0 ){
    return(cell.matrix)
  }
  
  # If all the cells are aneuploid, just return a matrix directly
  if(n.cells == total.cells){
    return(matrix(F, nrow=h, ncol=w))
  }
  
  cell.type = F # default to setting aneuploid
  if(n.cells > total.cells/2){ # if more than half the cells should be aneuploid,
    cell.type=T                # make the matrix aneuploid, and set euploids instead
    cell.matrix = matrix(F, nrow=h, ncol=w)
    n.cells = total.cells - n.cells # only perform the smaller operation
  }
  
  # Must have a minimum of 1 seed
  n.seeds = max(1, floor(dispersion * n.cells))
  
  # Set up seed cells randomly distributed
  seeds.to.add = n.seeds
  while(seeds.to.add>0){
    x = sample(1:w, size=1)
    y = sample(1:h, size=1)
    if(cell.matrix[y, x]!=cell.type){
      cell.matrix[y, x] = cell.type
      seeds.to.add = seeds.to.add - 1
    }
  }

  # cat(paste(n.seeds, " seeds created\n"))
  
  n.remaining = n.cells - n.seeds

  # Expand the seeds randomly for all remaining cells
  while(n.remaining>0){
    
    # Choose a random seed
    point = pick.seed(cell.matrix, cell.type)

    # Choose a random valid 8-connected point adjacent to the seed
    new.point = pick.adjacent(cell.matrix, point, !cell.type)
    
    # Update the matrix
    cell.matrix[new.point[['y']], new.point[['x']]] = cell.type
    n.remaining = n.remaining-1
  }
  
  return(cell.matrix)
}

# Calculate the number of euploids in the 5 cells at
# the centre of the matrix
make.count = function(prop.aneuploids, dispersion, make.chart=F, dim.x=10, dim.y=10, ...){
  
  # How many cells in the embryo at this stage?

  n.aneuploids = floor(prop.aneuploids * dim.x * dim.y)

  # True is euploid
  cells = matrix(T, nrow=dim.y, ncol=dim.x)
  
  result = make.aneuploids(cells, dispersion, n.aneuploids)
  
  data = as.data.frame(result)
  data$y = row(data)
  filt = tidyr::pivot_longer(data,cols=(-y), names_to="x", names_pattern="(\\d+)")
  filt$x = as.integer(filt$x)
  filt$y = filt$y[,1]
  
  # Sample 5 cells from the centre
  x.centre = ceiling(dim.x/2)
  y.centre = ceiling(dim.y/2)
  count = as.integer(result[y.centre, x.centre])+
    as.integer(result[y.centre+1, x.centre])+
    as.integer(result[y.centre-1, x.centre])+
    as.integer(result[y.centre, x.centre+1])+
    as.integer(result[y.centre, x.centre-1])
  
  p=NULL
  if(make.chart){
    p = ggplot(filt, aes(x=x, y=y, col=value))+
      geom_point(size=6)+
      geom_segment(aes(x =x.centre-1.5 , y = y.centre+0.5, xend = x.centre+1.5, yend = y.centre+0.5), col="black", size=2)+
      geom_segment(aes(x =x.centre-1.5 , y = y.centre-0.5, xend = x.centre+1.5, yend = y.centre-0.5), col="black", size=2)+
      geom_segment(aes(x =x.centre-0.5 , y = y.centre-1.5, xend = x.centre-0.5, yend = y.centre+1.5), col="black", size=2)+
      geom_segment(aes(x =x.centre+0.5 , y = y.centre-1.5, xend = x.centre+0.5, yend = y.centre+1.5), col="black", size=2)+
      geom_segment(aes(x =x.centre-1.5 , y = y.centre-0.5, xend = x.centre-1.5, yend = y.centre+0.5), col="black", size=2)+
      geom_segment(aes(x =x.centre+1.5 , y = y.centre-0.5, xend = x.centre+1.5, yend = y.centre+0.5), col="black", size=2)+
      geom_segment(aes(x =x.centre-0.5 , y = y.centre+1.5, xend = x.centre+0.5, yend = y.centre+1.5), col="black", size=2)+
      geom_segment(aes(x =x.centre-0.5 , y = y.centre-1.5, xend = x.centre+0.5, yend = y.centre-1.5), col="black", size=2)+
      scale_colour_manual(values = c("FALSE"="red", "TRUE"="darkgreen"))+
      theme_minimal()+
      labs(title=paste("Example plot:", count, "euploid cells in biopsy"))+
      theme(axis.text = element_blank(), 
            axis.title = element_blank(),
            axis.line = element_blank(),
            legend.position = "none",
            panel.grid = element_blank())

  }
  
  return(list("plot"=p, "n.euploid"=count))
}


####################################
# Running the analysis
####################################

# Make a grid of parameters
props = seq(0.0, 1.0, 0.2)
disps = seq(0.1, 1, 0.2) # dispersion - if 0, all one clump, if 1 all cells individual
combs = expand.grid(props, disps)

# Prepare result data frame
result = data.frame("prop"=NA, "disp"=NA, "e"=NA)

iterations = 100 # number of simulations per parameter combination

for(i in 1:nrow(combs)){
  cat("Combination", i, "Disp:", combs$Var2[i], "Prop:", combs$Var1[i],"\n")
  ptm <- proc.time()
  counts.and.plots = lapply(1:iterations, make.count,
                            prop.aneuploids=combs$Var1[i],
                            dispersion=combs$Var2[i], make.chart=T, dim.x=10, dim.y=10)
  time.taken = proc.time() - ptm
  cat("Combination", i, "took ", time.taken[['elapsed']], "seconds\n")
  # Extract counts and plots separately
  counts = sapply(1:iterations, function(i) counts.and.plots[[i]][['n.euploid']])
  plots  = lapply(1:iterations, function(i) counts.and.plots[[i]][['plot']])

  temp = data.frame("prop"=rep(combs$Var1[i], iterations),
                    "disp"=rep(combs$Var2[i], iterations),
                    "e" = counts)
  result = rbind(result, temp)
}

# Plot the results
result = na.omit(result)
ggplot(result, aes(x=e))+
  geom_bar()+
  facet_grid(prop~disp)+
  labs(title=paste("Cols = dispersal (0=clumped, 1=dispersed), rows = proportion aneuploid cells"),
       x="Number of euploid cells in 5 cell biopsy",
       y="Simulations in which this occurred")



####################################
# TESTS
####################################

create.test.matrix = function(){
  dim.x = 10
  dim.y = 10
  cell.matrix = matrix(T, nrow=dim.y, ncol=dim.x)
  cell.matrix[4, 4]=F
  cell.matrix[4, 5]=F
  cell.matrix[4, 6]=F
  cell.matrix[5, 4]=F
  cell.matrix[5, 5]=F
  cell.matrix[5, 6]=F
  cell.matrix[6, 4]=F
  cell.matrix[6, 5]=F
  cell.matrix[6, 6]=F
  
  cell.matrix[1, 1]=F
  cell.matrix[1, 2]=F
  cell.matrix[1, 10]=F
  cell.matrix[2, 1]=F
  cell.matrix[2, 2]=F
  cell.matrix[2, 10]=F
  cell.matrix[10, 1]=F
  cell.matrix[10, 2]=F
  cell.matrix[10, 10]=F
  
  return(cell.matrix)
}

# Test that make seed can fill the entire grid
# and produces the expected number of seeds
test.make.seed = function(){
  dim.x = 10
  dim.y = 10
  test.iterations = dim.x*dim.y
  for(n.aneuploid in 1:test.iterations){
    cell.matrix = matrix(T, nrow=dim.y, ncol=dim.x)
    for(i in 1:n.aneuploid){
      cell.matrix = make.seed(cell.matrix, F)
    }
    sum.aneuploids = sum(matrixStats::rowCounts(cell.matrix, value=F))
    assertthat::assert_that(are_equal(n.aneuploid, sum.aneuploids),
                            msg=paste("Exp:", n.aneuploid, "Act:", sum.aneuploids))
  }
}
# test.make.seed()

# Test that pick.seed can find all possible seeds
# Note this will not work until it has been corrected for
# the block filter
test.pick.seed = function(){
  dim.x = 10
  dim.y = 10
  test.iterations = dim.x*dim.y
  for(n.aneuploid in 1:test.iterations){
    cell.matrix = matrix(T, nrow=dim.y, ncol=dim.x)
    
    # Create the aneuploid cells
    for(i in 1:n.aneuploid){
      cell.matrix = make.seed(cell.matrix, F)
    }
    
    # Create a test matrix
    result.matrix = matrix(T, nrow=dim.y, ncol=dim.x)
    
    for(i in 1:1000){
      point = pick.seed(cell.matrix, F)
      result.matrix = set.value(result.matrix, point, F)
      if(i%%100==0){
        if(are_equal(cell.matrix, result.matrix)){ break}
      }
    }
    
    if(!are_equal(cell.matrix, result.matrix)){
      print(cell.matrix)
      print(result.matrix)
    }
    
    assertthat::assert_that(assertthat::are_equal(cell.matrix, result.matrix),
                            msg="Result matrix is not equal to cell matrix")
  }
}
# test.pick.seed()

test.block.has.value = function(){
  cell.matrix = create.test.matrix()
  print(cell.matrix)
  for(w in 1:ncol(cell.matrix)){
    for(h in 1:nrow(cell.matrix)){
      if((w==5 & h==5) | (w==1 & h==1)){
        assertthat::assert_that(block.has.value(cell.matrix, w, h),
                                msg = paste("Entire block does not have value in x", w, "y", h))
      }
      # There are also blocks with T values
      # else{
      #   assertthat::assert_that(!block.has.value(cell.matrix, w, h), 
      #                           msg = paste("Entire block has value in x", w, "y", h))
      # }
      
    }
    
  }
}
test.block.has.value()


test.pick.adjacent = function(){
  cell.matrix = create.test.matrix()
  
  print(cell.matrix)
  
  point = pick.adjacent(cell.matrix, list('x'=4, 'y'=4), T)
  print(paste("Point chosen: x:", point[['x']], "y:", point[['y']]))
}
# test.pick.adjacent()

test.get.points = function(){
  cell.matrix = create.test.matrix()
  print(cell.matrix)
  points = get.points(cell.matrix, F)
  print(points)
}
# test.get.points()
