# SPDE Animal Movement Sim Script
.libPaths(c('/zfs/users/joe.watson/joe.watson/rlibs2'))

# Simulate movement according to SPDE. For conservative system - link to IPP

# Load packages
library(tidyverse)
#library(magick)
library(sp)
library(spatstat)
library(raster)
#library(raster)
#library(inlabru)
library(maptools)
library(rgeos)
#library(gganimate)
library(abind)
#library(Rfast)
library(mgcv)
library(parallel)
library(foreach)
library(doParallel)

#library(dismo)

test_code=F

sigma <- sqrt(100)
mean_x <- 50
mean_y <- 50
sigma_observer <- sqrt(200)
sigma_Brownian <- sqrt(2)
sigma_Brownian_observer_lowbias <- sqrt(8)
sigma_Brownian_observer_highbias <- sqrt(2)

# Make SpatialPixelsDataFrame for both potential and long-run density
animal_potential <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))#raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
domain_boundary <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))
animal_potential <- makegrid(animal_potential,n=10000)
animal_UD <- animal_potential
animal_potential$z <- log(sqrt(1/(2*pi*sigma^2))) + ((-1/(2*sigma^2))*((animal_potential$x1-50)^2 + (animal_potential$x2-50)^2) )#exp( log((1/(2*pi*sigma^2))) + ((-1/(2*sigma^2))*((tc$x1-50)^2 + (tc$x2-50)^2) ) )
animal_UD$z <- exp(2*(log(sqrt(1/(2*pi*sigma^2))) + (-1/(2*sigma^2))*((animal_UD$x1-50)^2 + (animal_UD$x2-50)^2) )/ sigma_Brownian^2)#exp( 2 * exp( log((1/(2*pi*sigma^2))) + ((-1/(2*sigma^2))*((tc$x1-50)^2 + (tc$x2-50)^2) ) ) / sigma_Brownian^2 )
animal_UD$z <- animal_UD$z / sum(animal_UD$z)
animal_potential <- SpatialPointsDataFrame(coords=cbind(animal_UD$x1, animal_UD$x2), data=data.frame(z=animal_potential$z))
animal_potential <- as(animal_potential, 'SpatialPixelsDataFrame')
plot(animal_potential)
animal_UD <- SpatialPointsDataFrame(coords=cbind(animal_UD$x1, animal_UD$x2), data=data.frame(z=animal_UD$z))
animal_UD <- as(animal_UD, 'SpatialPixelsDataFrame')
plot(animal_UD, main='Long-run density of the animal')
#ggplot() + gg(animal_UD) + ggtitle('Long-run density of the animal') + xlab('x') + ylab('y') + scale_fill_viridis_c(name = expression(pi[m](x,y)))

# Make SpatialPixelsDataFrame for both potential and long-run density for observers
observer_potential_lowbias <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))#raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
domain_boundary <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))
observer_potential <- makegrid(domain_boundary,n=10000)
observer_UD_lowbias <- observer_potential
observer_UD_highbias <- observer_potential

#exp( 2 * ( log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((boundary_grid$x2-100)^2) ) ) / sigma_B^2 )
observer_potential$z <- log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((observer_potential$x2-100)^2) )#exp( log((1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((tc2$x2-100)^2) ) )
observer_UD_lowbias$z <- exp( 2 * ( log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((observer_UD_lowbias$x2-100)^2) ) ) / sigma_Brownian_observer_lowbias^2 )  #exp( 2 * exp( log((1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((tc2$x2-100)^2) ) ) / sigma_Brownian^2 )
observer_UD_lowbias$z <- observer_UD_lowbias$z / sum(observer_UD_lowbias$z)
observer_UD_lowbias <- SpatialPointsDataFrame(coords=cbind(observer_UD_lowbias$x1, observer_UD_lowbias$x2), data=data.frame(z=observer_UD_lowbias$z))
observer_UD_lowbias <- as(observer_UD_lowbias, 'SpatialPixelsDataFrame')
plot(observer_UD_lowbias)
observer_UD_highbias$z <- exp( 2 * ( log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((observer_UD_highbias$x2-100)^2) ) ) / sigma_Brownian_observer_highbias^2 )  #exp( 2 * exp( log((1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((tc2$x2-100)^2) ) ) / sigma_Brownian^2 )
observer_UD_highbias$z <- observer_UD_highbias$z / sum(observer_UD_highbias$z)
observer_UD_highbias <- SpatialPointsDataFrame(coords=cbind(observer_UD_highbias$x1, observer_UD_highbias$x2), data=data.frame(z=observer_UD_highbias$z))
observer_UD_highbias <- as(observer_UD_highbias, 'SpatialPixelsDataFrame')
plot(observer_UD_lowbias)
plot(observer_UD_highbias)
#ggplot() + gg(observer_UD_highbias) + ggtitle('Long-run density of the observers') + xlab('x') + ylab('y') + scale_fill_viridis_c(name = expression(pi[xi](x,y)))

# write function for solving roots of quadratic
quad_fun <- function(x, theta){
  return((theta[1] + theta[2]*x[1] + theta[3]*x[2] + theta[4]*(x[1]^2) + theta[5]*(x[2]^2))^2 )
}
quad_fun_solve <- function(theta){
  return(optim(par = c(50,50), quad_fun, theta=theta))
}

# write function for computing proportion of area of intersection of a set of circles located
# at observer_coords, each with radius radii
prop_area_int_circ <- function(radii, observer_coords)
{
  d_circles <- apply(observer_coords, c(2), dist)
  # this returns the matrix of 0.5n(n-1) pairwise distances
  n_circles <- dim(d_circles)[1]
  res <- 2*radii^2*acos(d_circles/(2*radii))-(0.5)*d_circles*sqrt(4*radii^2-d_circles^2)
  # returns matrix of areal intersections
  res[is.na(res)] <- 0
  # need to divide mean areal intersection by area of circle
  return(mean(res)/(pi*radii^2))
}

potential <- function(x, y, sigma, mean){ return( log(sqrt(1/(2*pi*sigma^2))) + ((-1/(2*sigma^2))*((x-mean[1])^2 + (y-mean[2])^2) ) )  }
deriv_potential <- deriv( (P ~ log(sqrt(1/(2*pi*sigma^2))) + ((-1/(2*sigma^2))*((x-mean_x)^2 + (y-mean_y)^2) ) ) ,
                          c('x','y'), func=T)
long_run <-  function(x, y, sigma, mean){exp(2*potential(x, y, sigma, mean)/sigma_Brownian^2)}
potential_observers <- function(y){log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((y-100)^2) ) }
deriv_potential_observers <- deriv( (P ~ log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((y-100)^2) ) ) ,
                                    c('x','y'), func=T)
#long_run_observers <-  function(x, y, sigma, mean){exp(2*potential_observers(y)/sigma_Brownian_observer^2)}
deriv_potential_fun <- function(x, y){return(attr(deriv_potential(x=x, y=y),"gradient"))}
deriv_potential_observers_fun <- function(x, y){return(attr(deriv_potential_observers(x=x, y=y),"gradient"))}

species <- setClass("species", slots=c(x="numeric", y="numeric"))
observers <- setClass("observers", slots=c(x="numeric", y="numeric", n='numeric', static="logical"))

center_domain <- gCentroid(as(observer_UD_highbias, 'SpatialPolygonsDataFrame'))

# define distance sampling function as a function of r (max distance)
distance_function <- function(d,r)
{
  #return(1-1.5*(d/r)+0.5*(d/r)^3)
  return(pmax(1-(d/r),0))
}

# function for simulating an animal and observers
go <- function (sp, sp2, n, sigma, sigma_observer, n_keep, Brillinger=T, r=1e9, dist_fun=function(x){1}, observers_UD) {
  # animal trajectory objects and initialise
  track <- matrix(0, nrow = n, ncol = 2)
  # randomly initialise the starting location according to stationary density
  track[1,] <- animal_UD@coords[ sample.int(n=length(animal_UD@data$z), size=1, prob=animal_UD@data$z) ,]
  #track[1,2] <- sp@y
  x_velocity_old <- 0
  y_velocity_old <- 0
  velocities <- matrix(0, nrow = n, ncol = 2)
  velocities[1,1] <- 0
  velocities[1,2] <- 0

  #browser()

  # observers trajectories objects and initialise
  track_observers <- array(0, dim=c(sp2@n, n, 2))
  # randomly initialise the starting location according to stationary density
  track_observers[,1,] <- observers_UD@coords[ sample.int(n=length(observers_UD@data$z), size=sp2@n, prob=observers_UD@data$z) ,]
  #track_observers[,1,2] <- sp2@y
  x_velocity_old_observers <- rep(0, sp2@n)
  y_velocity_old_observers <- rep(0, sp2@n)
  velocities_observers <- array(0, dim=c(sp2@n, n, 2))
  velocities_observers[,1,1] <- 0
  velocities_observers[,1,2] <- 0

  # which observers are static?
  static <- sp2@static
  for (step in 2:n) {
    lon_candidate <- -999999
    lat_candidate <- -999999

    x <- track[step-1,1]
    y <- track[step-1,2]
    gradient_potential <- deriv_potential_fun(x, y)

    #browser()
    lon_candidate_observers <- rep(-999999, sp2@n)
    lat_candidate_observers <- rep(-999999, sp2@n)

    x_observers <- track_observers[,step-1,1]
    y_observers <- track_observers[,step-1,2]
    gradient_potential_observers <- deriv_potential_observers_fun(x_observers, y_observers)

    # Compute distance between observers and individual
    dist_mat <- as.matrix(dist(rbind(cbind(x,y),cbind(x_observers, y_observers))))

    # first check if any observer has spotted individual
    # if euclidean distance between individual and observer is less
    # than r then sample bernoulli probability based on distance sampling function
    detected <- 0
    if(min(dist_mat[1,-1]) < r)
    {
      # Which observers are within r units away?
      possible_sighting_ind <- which(dist_mat[1,-1] < r)

      for(j in possible_sighting_ind)
      {
        # distance from observer?
        distance <- dist_mat[1,-1][j]
        # terminate if detected
        detected <- detected + rbernoulli(n=1, p = dist_fun(d=distance,r=r))
      }

      if(detected > 0){

        thinned_points <- track[seq(from=1, to=step-1, by = n_keep),]
        thinned_points_observers <- track_observers[,seq(from=1, to=step-1, by = n_keep),,drop=FALSE]
        # compute average overlap between observers' fields-of-view
        prop_intersection <- suppressWarnings(prop_area_int_circ(radii = r, track_observers[,c(1:(step-1)),,drop=FALSE]))
        return(list(track=track[c(1:(step-1)),], thinned_points=thinned_points, velocities=velocities[c(1:(step-1)),],
                    track_observers=track_observers[,c(1:(step-1)),,drop=FALSE],
                    thinned_points_observers=thinned_points_observers,
                    velocities_observers=velocities_observers[,c(1:(step-1)),,drop=FALSE],
                    detected_point=track[step-1,],
                    prop_intersection=prop_intersection))}
    }

    if(Brillinger)
    {
      while ( !gWithin(SpatialPoints(cbind(lon_candidate, lat_candidate)), domain_boundary) ) {
        x_velocity <- gradient_potential[1]
        y_velocity <- gradient_potential[2]

        lon_candidate <- track[step-1,1] + (sigma * rnorm(1)) + x_velocity
        lat_candidate <- track[step-1,2] + (sigma * rnorm(1)) + y_velocity
      }
    }

    track[step,1] <- lon_candidate
    track[step,2] <- lat_candidate

    velocities[step,1] <- x_velocity
    velocities[step,2] <- y_velocity

    # Repeat for observers
    #browser()
    if(Brillinger)
    {
      while ( !gWithin(SpatialPoints(cbind(lon_candidate_observers, lat_candidate_observers)), domain_boundary) ) {
        x_velocity_observers <- gradient_potential_observers[,1]
        y_velocity_observers <- gradient_potential_observers[,2]

        lon_candidate_observers <- track_observers[,step-1,1] + (sigma_observer * rnorm(sp2@n)) + x_velocity_observers
        lat_candidate_observers <- track_observers[,step-1,2] + (sigma_observer * rnorm(sp2@n)) + y_velocity_observers
      }
    }
    track_observers[,step,1] <- lon_candidate_observers
    track_observers[,step,2] <- lat_candidate_observers

    velocities_observers[,step,1] <- x_velocity_observers
    velocities_observers[,step,2] <- y_velocity_observers

    # for the static observers fix their location equal to the previous one
    track_observers[static,step,1] <- track_observers[static,(step-1),1]
    track_observers[static,step,2] <- track_observers[static,(step-1),2]

  }
  thinned_points <- track[seq(from=1, to=n, by = n_keep),]
  thinned_points_observers <- track_observers[,seq(from=1, to=n, by = n_keep),]
  # compute average overlap between observers' fields-of-view
  prop_intersection <- suppressWarnings(prop_area_int_circ(radii = r, track_observers))
  return(list(track=track, thinned_points=thinned_points, velocities=velocities,
              track_observers=track_observers, thinned_points_observers=thinned_points_observers,
              velocities_observers=velocities_observers,
              prop_intersection=prop_intersection))
}

animal <- species(x= 50, y =50)  # starts at the mode
one_boat <- observers(x=75,y=50, n=1, static=F)
two_boat_one_static <- observers(x=c(75,75,75),y=c(50,50,50), n=3, static=c(F,T,F))
twenty_boats <- observers(x=rep(75,20),y=rep(50,20), n=20, static=rep(F,20))

# Test the go function
if(test_code == T)
{
  simul <- go(animal,twenty_boats, 500, sigma=sigma_Brownian, sigma_observer = sigma_Brownian_observer_lowbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_lowbias)

  plot(observer_UD_lowbias)
  lines(simul$track, lwd=0.1, col="white")
  lines(simul$track_observers[1,,], lwd=10, col="green")
  lines(simul$track_observers[2,,], lwd=10, col="green")
  lines(simul$track_observers[3,,], lwd=10, col="green")
  lines(simul$track_observers[4,,], lwd=10, col="green")
  lines(simul$track_observers[5,,], lwd=10, col="green")
  lines(simul$track_observers[6,,], lwd=10, col="green")
  lines(simul$track_observers[7,,], lwd=10, col="green")
  lines(simul$track_observers[8,,], lwd=10, col="green")
  lines(simul$track_observers[9,,], lwd=10, col="green")
  lines(simul$track_observers[10,,], lwd=10, col="green")
  lines(simul$track_observers[11,,], lwd=10, col="green")
  lines(simul$track_observers[12,,], lwd=10, col="green")
  lines(simul$track_observers[13,,], lwd=10, col="green")
  lines(simul$track_observers[14,,], lwd=10, col="green")
  lines(simul$track_observers[15,,], lwd=10, col="green")
  lines(simul$track_observers[16,,], lwd=10, col="green")
  lines(simul$track_observers[17,,], lwd=10, col="green")
  lines(simul$track_observers[18,,], lwd=10, col="green")
  lines(simul$track_observers[19,,], lwd=10, col="green")
  lines(simul$track_observers[20,,], lwd=10, col="green")

  if(!is.null(simul$detected_point))
  {
    points(x=simul$detected_point[1],y=simul$detected_point[2], lwd=10, col='black')
  }

  # how fast is the average boat travelling?
  velocity_boats_lowbias <- rep(NA, 50)
  velocity_boats_highbias <- rep(NA, 50)
  velocity_animal <- rep(NA, 50)
  for(i in 1:50)
  {
    print(i)
    simul <- go(animal,twenty_boats, 500, sigma=sigma_Brownian, sigma_observer = sigma_Brownian_observer_lowbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_lowbias)
    simul2 <- go(animal,twenty_boats, 500, sigma=sigma_Brownian, sigma_observer = sigma_Brownian_observer_highbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_highbias)
    if(dim(simul$track_observers)[2] > 1) #movement happened
    {
      velocity_boats_lowbias[i] <- mean(sqrt(apply(apply(simul$track_observers, c(1,3), FUN = diff)^2, c(1,2), sum)))
      velocity_animal[i] <- mean(sqrt(apply(diff(simul$track)^2, 1, sum)))
    }
    if(dim(simul2$track_observers)[2] > 1) #movement happened
    {
      velocity_boats_highbias[i] <- mean(sqrt(apply(apply(simul2$track_observers, c(1,3), FUN = diff)^2, c(1,2), sum)))
      velocity_animal[i] <- mean(sqrt(apply(diff(simul2$track)^2, 1, sum)))
    }
  }
  hist(velocity_boats_highbias)
  hist(velocity_boats_lowbias)
  hist(velocity_animal)
  summary(velocity_boats_highbias)
  summary(velocity_boats_lowbias)
  summary(velocity_animal)
  mean(velocity_boats_highbias/velocity_animal, na.rm=T) # same speed at ~1.75 units per time step
  mean(velocity_boats_lowbias/velocity_animal, na.rm=T) # boat is 1.85 times faster at ~3.5 units per time step


  simul2 <- go(animal, two_boat_one_static, 500, sigma=sigma_Brownian, sigma_observer = sigma_Brownian_observer_highbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_highbias)
  plot(observer_UD_highbias)
  lines(simul2$track, lwd=0.1, col="red")
  lines(simul2$track_observers[1,,], lwd=10, col="green")
  lines(simul2$track_observers[2,,], lwd=10, col="yellow")
  lines(simul2$track_observers[3,,], lwd=10, col="pink")
  if(!is.null(simul2$detected_point))
  {
    points(x=simul2$detected_point[1],y=simul2$detected_point[2], lwd=10, col='black')
  }
}

# Now write a script to run this simulation M times,
# i.e. M independent encounters or periods
# store the detection locations and store the
# observer tracks as a SpatialPoints object
simulator <- function(M, sp, sp2, n, sigma, sigma_observer, n_keep, Brillinger=T, r=1e9, dist_fun=function(x){1}, observers_UD)
{
  observer_tracks <- vector(mode = "list", length = sp2@n)
  observer_tracks <- rep(list(observer_tracks), M)
  encounter_locations <- NULL
  animal_tracks <- NULL
  prop_intersection <- rep(NA, M)

  for(i in 1:M)
  {
    #print(paste('iteration', i, 'of', M))
    result <- go(sp, sp2, n, sigma, sigma_observer, n_keep,
                 Brillinger, r, dist_fun, observers_UD)
    for(j in 1:sp2@n)
    {
      #browser()
      observer_tracks[[i]][[j]] <- matrix(result$track_observers[j,,], ncol=2)#abind(observer_tracks, result$track_observers, along = 2)
      observer_tracks[[i]][[j]] <- SpatialPoints(coords = observer_tracks[[i]][[j]])
    }

    animal_tracks <- rbind(animal_tracks, result$thinned_points)
    prop_intersection[i] <- result$prop_intersection

    if(!is.null(result$detected_point))
    {
      encounter_locations <- rbind(encounter_locations, result$detected_point)
    }

  }
  #dim(observer_tracks) <- c(dim(observer_tracks)[1]*dim(observer_tracks)[2], dim(observer_tracks)[3])
  #observer_tracks <- SpatialPoints(coords = observer_tracks)
  animal_tracks <- SpatialPoints(coords = animal_tracks)
  return(list(observer_tracks=observer_tracks, encounter_locations=encounter_locations,
              animal_tracks=animal_tracks, prop_intersection=prop_intersection))
}


if(test_code==T)
{
  # Test the simulator function
  results <- simulator(M=10, animal, one_boat, n=500, sigma = sigma_Brownian, sigma_observer = sigma_Brownian_observer_lowbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_lowbias)

  plot(observer_UD_lowbias)
  lines(results$observer_tracks@coords, lwd=5, col="green")
  if(!is.null(results$encounter_locations))
  {
    points(x=results$encounter_locations[,1],y=results$encounter_locations[,2], lwd=10, col='black')
  }
  # Create an effort covariate - scaled by the distance function
  effort <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))#raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
  effort <- makegrid(effort,n=10000)
  #effort <- as(SpatialPoints(effort), 'ppp')

  effort <- as(SpatialPoints(effort), 'SpatialPixels')
  # compute distance matrix between pixels
  dist_mat_pixels <- readRDS('dist_mat_pixels.rds')#Dist(effort@coords)
  # map the observed tracks to the closest pixels
  # map the observed tracks to the closest pixels
  approx_tracks <- lapply(results$observer_tracks, FUN = function(x){lapply(x, FUN = function(y){over(y, effort)})}) #over(results$observer_tracks, effort)
  # compute the approx distance between each observer track point and the pixels
  # Need to loop to avoid memory crash
  effort_val <- matrix(0, nrow=dim(dist_mat_pixels)[1], ncol = length(approx_tracks))
  for(ind2 in 1:length(approx_tracks))
  {
    for(ind in 1:(dim(dist_mat_pixels)[1]))
    {
      effort_vals[ind, ind2] <- sum(1 - apply(matrix(sapply(approx_tracks[[ind2]],
                                                            FUN = function(x){distance_function(d=dist_mat_pixels[ind,x],r=10)}),ncol=length(approx_tracks[[ind2]])),
                                              1, function(y){prod(1-y)}))
    }
  }

  effort_vals <- apply(effort_vals, 1, sum)
  effort$effort <- log(effort_vals+1e-5)

  # now create ppp object for animal tracks
  animal_tracks <- as(results$animal_tracks, 'ppp')

  #ggplot() + gg(effort)
  #ggplot() + gg(results$observer_tracks)

  # Fit the model in spatstat and use kernel density estimate to estimate intensity
  effort_im <- as(as(effort, 'SpatialGridDataFrame'),'im')
  encounters_ppp <- as(SpatialPoints(results$encounter_locations),'ppp')
  window <- as(domain_boundary, 'owin')
  encounters_ppp$window <- window
  effort_predict <- effort_im
  effort_predict$v = rep(0, length(effort_predict$v))
  animal_tracks$window <- window
  # Smooth to get the long run density
  animal_longrun <- as(density(animal_tracks) ,'SpatialGridDataFrame')
  animal_longrun$v <- animal_longrun$v / sum(animal_longrun$v)

  model <- ppm(Q=encounters_ppp ~ offset(effort) + x+y+I(x^2)+I(y^2), data=list(effort=effort_im ))
  model2 <- ppm(Q=encounters_ppp ~  x+y+I(x^2)+I(y^2))
  plotting_coef <- coef(model)
  quad_fun_solve(theta = model$coef)$par
  quad_fun_solve(theta = model2$coef)$par

  #plotting_coef[2] <- 0
  # predict the normalized invariant density for the models
  pred_model <- predict(model, type='intensity',new.coef=plotting_coef, locations=effort_predict, covariates=data.frame(effort=effort_predict$v))
  pred_model$v <- pred_model$v / sum( pred_model$v)
  pred_model2 <- predict(model2, type='intensity', locations=effort_predict)
  pred_model2$v <- pred_model2$v / sum( pred_model2$v)

  multiplot(
    ggplot() + gg(as(pred_model,'SpatialGridDataFrame')) + theme(legend.position = "none"),
    ggplot() + gg(as(pred_model2,'SpatialGridDataFrame')) + theme(legend.position = "none"),
    ggplot() + gg(observer_UD_lowbias) + theme(legend.position = "none"),
    ggplot() + gg(animal_longrun) + theme(legend.position = "none"),
    layout = matrix(c(1:4), nrow=2,ncol=2, byrow=T))

  # Now try with 1 mobile and 20 static observers
  one_boat_twenty_static <- observers(x=rep(75,21),y=rep(50,21), n=21, static=c(F,rep(T, 20)))

  results2 <- simulator(M=10, animal, one_boat_twenty_static, n=500, sigma = sigma_Brownian, sigma_observer = sigma_Brownian_observer_lowbias, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD = observer_UD_lowbias)

  # Create an effort covariate - scaled by the distance function
  effort2 <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))#raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
  effort2 <- makegrid(effort2,n=10000)
  #effort <- as(SpatialPoints(effort), 'ppp')

  effort2 <- as(SpatialPoints(effort2), 'SpatialPixels')

  # map the observed tracks to the closest pixels
  #approx_tracks2 <- over(results2$observer_tracks, effort2)
  # compute the approx distance between each observer track point and the pixels
  # Need to loop to avoid memory crash
  effort_vals2 <- rep(0, dim(dist_mat_pixels)[1])

  # map the observed tracks to the closest pixels
  approx_tracks <- lapply(results2$observer_tracks, FUN = function(x){lapply(x, FUN = function(y){over(y, effort2)})}) #over(results$observer_tracks, effort)
  # compute the approx distance between each observer track point and the pixels
  # Need to loop to avoid memory crash
  effort_vals <- matrix(0, nrow=dim(dist_mat_pixels)[1], ncol = length(approx_tracks))
  for(ind2 in 1:length(approx_tracks))
  {
    for(ind in 1:(dim(dist_mat_pixels)[1]))
    {
      effort_vals[ind, ind2] <- sum(1 - apply(matrix(sapply(approx_tracks[[ind2]],
                                             FUN = function(x){distance_function(d=dist_mat_pixels[ind,x],r=10)}),ncol=length(approx_tracks[[ind2]])),
                                      1, function(y){prod(1-y)}))
    }
  }

  effort_vals <- apply(effort_vals, 1, sum)
  effort2$effort <- log(effort_vals+1e-5)

  # now create ppp object for animal tracks
  animal_tracks2 <- as(results2$animal_tracks, 'ppp')

  ggplot() + gg(effort2)
  #ggplot() + gg(results2$observer_tracks)

  # Fit the model in spatstat and use kernel density estimate to estimate intensity
  effort_im2 <- as(as(effort2, 'SpatialGridDataFrame'),'im')
  encounters_ppp2 <- as(SpatialPoints(results2$encounter_locations),'ppp')
  encounters_ppp2$window <- window
  effort_predict2 <- effort_im2
  effort_predict2$v = rep(0, length(effort_predict2$v))
  animal_tracks2$window <- window
  # Smooth to get the long run density
  animal_longrun2 <- as(density(animal_tracks2) ,'SpatialGridDataFrame')
  animal_longrun2$v <- animal_longrun2$v / sum(animal_longrun2$v)

  model3 <- ppm(Q=encounters_ppp2 ~ offset(effort) + x+y+I(x^2)+I(y^2), data=list(effort=effort_im2 ))
  model4 <- ppm(Q=encounters_ppp2 ~  x+y+I(x^2)+I(y^2))
  plotting_coef2 <- coef(model3)
  #plotting_coef[2] <- 0
  # predict the normalized invariant density for the models
  pred_model3 <- predict(model3, type='intensity',new.coef=plotting_coef2, locations=effort_predict2, covariates=data.frame(effort=effort_predict2$v))
  pred_model3$v <- pred_model3$v / sum( pred_model3$v)
  pred_model4 <- predict(model4, type='intensity', locations=effort_predict2)
  pred_model4$v <- pred_model4$v / sum( pred_model4$v)

  multiplot(
    ggplot() + gg(as(pred_model3,'SpatialGridDataFrame')) + theme(legend.position = "none"),
    ggplot() + gg(as(pred_model4,'SpatialGridDataFrame')) + theme(legend.position = "none"),
    ggplot() + gg(observer_UD_lowbias) + theme(legend.position = "none"),
    ggplot() + gg(animal_UD) + theme(legend.position = "none"),
    layout = matrix(c(1:4), nrow=2,ncol=2, byrow=T))

  #setwd("~/ownCloud/Whale Project/Encounters Simulations")
  #saveRDS(object = list(results=results, results2=results2), 'Encounters_results.rds')
}


######## Write the simulation study script ########
# This requires simulating the possible methods taken by the analysts
# set-up 4 different observers
twenty_static <- observers(x=rep(75,20),y=rep(50,20), n=20, static=c(rep(T, 20)))
one_boat_twenty_static <- observers(x=rep(75,21),y=rep(50,21), n=21, static=c(F,rep(T, 20)))
one_boat <- observers(x=rep(75,1),y=rep(50,1), n=1, static=F)
twenty_boats <- observers(x=rep(75,20),y=rep(50,20), n=20, static=rep(F, 20))
observers_sim <- list(twenty_boats)#one_boat, twenty_static, twenty_boats, one_boat_twenty_static)
# set-up different observer biases
sigma_Brownian_observer_sim <- c(sqrt(8))#, sqrt(2))
# set-up different numbers of trips
n_trips_sim <- c(150, 300)
# do we assume observer effort or perfect detectability? Truth remains same with a linear decaying detection probability surface
detect_surface_sim <- c(T)#c(F, T)
# What is the assumed detection range of the observers? The true detection function remains same at 10
detect_range_sim <- 10#c(2, 10, 50)
# How many replications per simulation setting?
n_rep_sim <- 50

boundary_polygons <- SpatialPolygons(Srl=list(Polygons(srl=list(Polygon(cbind(c(0,100,100,0),c(0,0,100,100)))), ID='1')))#raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
boundary_grid <- makegrid(boundary_polygons, n=10000)
boundary_grid_pixels <- as(SpatialPoints(boundary_grid), 'SpatialPixels')
# compute distance matrix between pixels
dist_mat_pixels <- readRDS('dist_mat_pixels.rds')#Dist(boundary_grid_pixels@coords)

# Write the simulation function that runs simulator and go and then compiles the result
# as if an analyst had the data. Then, compare the analysts results with the truth
simulation_study <- function(observers_val, sigma_Brownian_observer_val, n_trips_val,
                             distance_function_assumed, detect_range_val, n_rep_val,
                             sp, n, sigma, n_keep, Brillinger=T, observers_UD
)
{

  results <- simulator(M=n_trips_val, animal, observers_val, n=500, sigma = sigma_Brownian, sigma_observer = sigma_Brownian_observer_val, n_keep = 10, Brillinger=T, r=10, dist_fun = distance_function, observers_UD=observers_UD)

  # Create an effort covariate - scaled by the distance function
  effort <- boundary_grid_pixels

  # map the observed tracks to the closest pixels
  approx_tracks <- lapply(results$observer_tracks, FUN = function(x){lapply(x, FUN = function(y){over(y, effort)})}) #over(results$observer_tracks, effort)
  # compute the approx distance between each observer track point and the pixels
  # Need to loop to avoid memory crash
  effort_vals <- matrix(0, nrow=dim(dist_mat_pixels)[1], ncol = length(approx_tracks))
  for(ind2 in 1:length(approx_tracks))
  {
    for(ind in 1:(dim(dist_mat_pixels)[1]))
    {
      effort_vals[ind, ind2] <- sum(1 - apply(matrix(sapply(approx_tracks[[ind2]],
                                                            FUN = function(x){distance_function(d=dist_mat_pixels[ind,x],r=10)}),ncol=length(approx_tracks[[ind2]])),
                                              1, function(y){prod(1-y)}))
    }
  }

  effort_vals <- apply(effort_vals, 1, sum)
  effort$effort <- log(effort_vals+1e-5)

  # now create ppp object for animal tracks
  animal_tracks <- as(results$animal_tracks, 'ppp')

  # Fit the model in spatstat and use kernel density estimate to estimate intensity
  effort_im <- as(as(effort, 'SpatialGridDataFrame'),'im')
  encounters_ppp <- as(SpatialPoints(results$encounter_locations),'ppp')
  window <- as(domain_boundary, 'owin')
  encounters_ppp$window <- window
  effort_predict <- effort_im
  effort_predict$v = rep(0, length(effort_predict$v))
  animal_tracks$window <- window
  # Smooth to get the long run density
  animal_longrun <- as(density(animal_tracks) ,'SpatialGridDataFrame')
  animal_longrun$v <- animal_longrun$v / sum(animal_longrun$v)

  library(mgcv)
  model <- ppm(Q=encounters_ppp ~ offset(effort) + x+y+I(x^2)+I(y^2), data=list(effort=effort_im ))
  model2 <- ppm(Q=encounters_ppp ~  x+y+I(x^2)+I(y^2))
  plotting_coef <- coef(model)
  estimated_UDcenter_model <- quad_fun_solve(theta = model$coef)$par
  estimated_UDcenter_model2 <- quad_fun_solve(theta = model2$coef)$par

  #plotting_coef[2] <- 0
  # predict the normalized invariant density for the models
  pred_model <- predict(model, type='intensity',new.coef=plotting_coef, locations=effort_predict, covariates=data.frame(effort=effort_predict$v))
  pred_model$v <- pred_model$v / sum( pred_model$v)
  pred_model2 <- predict(model2, type='intensity', locations=effort_predict)
  pred_model2$v <- pred_model2$v / sum( pred_model2$v)

  # multiplot(
  #   ggplot() + gg(as(pred_model,'SpatialGridDataFrame')) + theme(legend.position = "none"),
  #   ggplot() + gg(as(pred_model2,'SpatialGridDataFrame')) + theme(legend.position = "none"),
  #   ggplot() + gg(tc_2) + theme(legend.position = "none"),
  #   ggplot() + gg(animal_longrun) + theme(legend.position = "none"),
  #   layout = matrix(c(1:4), nrow=2,ncol=2, byrow=T))

  # compute Prediction variance, Preidction Bias and UD center bias and variance
  PV <- mean((pred_model$v - animal_UD$z)^2)
  PB <- mean(pred_model$v - animal_UD$z)
  UD_center_bias <- estimated_UDcenter_model - c(50,50)
  #UD_center_variance <- (estimated_UDcenter_model - c(50,50))^2

  PV2 <- mean((pred_model2$v - animal_UD$z)^2)
  PB2 <- mean(pred_model2$v - animal_UD$z)
  UD_center_bias2 <- estimated_UDcenter_model2 - c(50,50)
  #UD_center_variance2 <- (estimated_UDcenter_model2 - c(50,50))^2
  prop_intersection <- mean(results$prop_intersection)

  return(c(PV, PV2, PB, PB2, UD_center_bias, UD_center_bias2, prop_intersection))

}

if(test_code)
{
  # test the simulation_study function
  sigma_Brownian_observer <- sigma_Brownian_observer_sim[1]
  # update the observer field
  tc2_2$z <- exp( 2 * ( log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((boundary_grid$x2-100)^2) ) ) / sigma_Brownian_observer^2 )  #exp( 2 * exp( log((1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((tc2$x2-100)^2) ) ) / sigma_Brownian^2 )
  tc2_2$z <- tc2_2$z / sum(tc2_2$z)

  long_run_observers <-  function(x, y, sigma, mean){exp(2*potential_observers(y)/sigma_Brownian_observer^2)}
  distance_function_assumed <- function(d,r)
  {
    return(ceiling(pmax(1-(d/r),0)))
  }

  results_test_sim <- simulation_study(observers_val = observers_sim[[1]],
                                       sigma_Brownian_observer_val = sigma_Brownian_observer_sim[1],
                                       n_trips_val=n_trips_sim[1],
                                       distance_function_assumed = distance_function_assumed,
                                       detect_range_val=detect_range_sim[1],
                                       n_rep_val=n_rep_sim[1],
                                       sp=animal, n=500, sigma=sigma_Brownian,
                                       n_keep=1000, Brillinger=T)
}


### Now conduct the simulation study
# create objects for storing results
Sim_Results <- array(0, dim = c(length(observers_sim), length(sigma_Brownian_observer_sim), length(n_trips_sim), length(detect_surface_sim), length(detect_range_sim), n_rep_sim, 9),
                     dimnames = list(Observers=c('twenty_boats'),
                                     Sigma_Observer=as.character(sigma_Brownian_observer_sim),
                                     N_Trips=as.character(n_trips_sim),
                                     Detect_surface=as.character(detect_surface_sim),
                                     Detect_Range=as.character(detect_range_sim),
                                     N_rep = as.character(1:n_rep_sim),
                                     Quantity=c('PV', 'PV2', 'PB', 'PB2', 'UD_center_bias_x',
                                                'UD_center_bias_y','UD_center_bias2_x',
                                                'UD_center_bias2_y','Prop_Intersection')))

# register parallel sockets
registerDoParallel(cores=25)
tc2_2 <- observer_UD_lowbias
# loop over all simulation settings
count_iter <- 1
max_iter <- prod(dim(Sim_Results)[-c(6,7)])
count1 <- 1
for(observers in observers_sim)
{
  count2 <- 1
  for(sigma_B in sigma_Brownian_observer_sim)
  {
    # update the observer field
    tc2_2$z <- exp( 2 * ( log(sqrt(1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((boundary_grid$x2-100)^2) ) ) / sigma_B^2 )  #exp( 2 * exp( log((1/(2*pi*sigma_observer^2))) + ((-1/(2*sigma_observer^2))*((tc2$x2-100)^2) ) ) / sigma_Brownian^2 )
    tc2_2$z <- tc2_2$z / sum(tc2_2$z)

    long_run_observers <-  function(x, y, sigma, mean){exp(2*potential_observers(y)/sigma_B^2)}
    count3 <- 1
    for(n_trips in n_trips_sim)
    {
      count4 <- 1
      for(detect_surface in detect_surface_sim)
      {
        # update the detection function assumed?
        if(detect_surface){
          distance_function_assumed <- function(d,r)
          {
            return(pmax(1-(d/r),0))
          }
        }
        if(!detect_surface){
          distance_function_assumed <- function(d,r)
          {
            return(ceiling(pmax(1-(d/r),0)))
          }
        }
        count5 <- 1
        for(detect_range in detect_range_sim)
        {
          count6 <- 1
          print(paste('iteration', count_iter, 'out of', max_iter))
          Sim_Results[count1,count2,count3,count4,count5,,] <-
            foreach(nrep=1:n_rep_sim, .combine='rbind') %dopar%
          {
            simulation_study(observers_val = observers,
                               sigma_Brownian_observer_val = sigma_B,
                               n_trips_val=n_trips,
                               distance_function_assumed = distance_function_assumed,
                               detect_range_val=detect_range,
                               n_rep_val=n_rep,
                               sp=animal, n=500, sigma=sigma_Brownian,
                               n_keep=1000, Brillinger=T, observers_UD = tc2_2)
            #count6 <- count6 + 1

          }
          count5 <- count5 + 1
          count_iter <- count_iter + 1
        }
        count4 <- count4 + 1
      }
      count3 <- count3 + 1
    }
    count2 <- count2 + 1
    # save intermediate results
    saveRDS(Sim_Results, 'Encounter_Sim_Results_New_overlapcorrect2.rds')
  }
  count1 <- count1 + 1
}

saveRDS(Sim_Results, 'Encounter_Sim_Results_New_overlapcorrect2.rds')
