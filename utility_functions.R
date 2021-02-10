# Utility functions for LGCP INLA - mainly mapping satellite rasters onto mesh

####### Interpolator function needed for INLAbru -- change elevation to the desired dataset
f.interp <- function(sp_points,sp_grid) {
  # turn SpatialPointsDataFrame into SpatialPoints object - i.e. extract coords:
  spp <- SpatialPoints(data.frame(x=coordinates(sp_points)[,1],
                                  y=coordinates(sp_points)[,2])) 
  # attach the appropriate coordinate reference system (CRS)
  proj4string(spp) <- CRS(proj4string(sp_points))
  # save this for later (gDistance function does not work in long/lat)
  projection = CRS(proj4string(sp_points))
  # change to the sp_grid's coordinate reference system (CRS)
  spp = spTransform(spp, CRSobj = CRS(proj4string(sp_grid)))
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, sp_grid) 
  #browser()
  # Detect NAs
  nbhd = 9 # neighbourhood size needed to remove all NAs
  NA_ind = as.matrix(is.na(v))
  print(paste('NA matrix is of dimension',dim(NA_ind)))
  if(sum(NA_ind)>0) # NAs are problematic
  {
    # Need to extract the coordinates of the sp_grid and put into SpatialPoints 
    grid_points = SpatialPoints(coords = coordinates(sp_grid),
                                proj4string = CRS(proj4string(sp_grid)))
    # Need to ensure we project away from long/lat
    grid_points = spTransform(grid_points, projection)
    spp = spTransform(spp, CRSobj = projection)
    #loop over the different variables in sp_grid and match missing values to closest non-missing polygon value.
    missing_vars = which(colSums(NA_ind)>0) # which variables have missing values
    count = 1 # progress report
    for(i in missing_vars)
    {
      print(paste(count,'out of',length(missing_vars),'missing variables processed')) 
      NAs_gone_ind = F
      while(NAs_gone_ind == F)
      {
        v2 = v[,i] # subset the ith variable
        NA_ind2 = which(is.na(v2))
        dists = gDistance( grid_points, spp[NA_ind2,], byid = T )
        closest_inds = apply(dists, 1, FUN = function(x){sort(x, decreasing = F, index.return=T)$ix[1:nbhd]})#pick 9 closest spots
        v[NA_ind2,i] = apply(as.matrix(closest_inds), 2, FUN = function(x){median(sp_grid@data[closest_inds,i], na.rm = T)}) # take the median of ith variable
        if(sum(is.na(v[,i])) == 0)
        {
          NAs_gone_ind = T # we've successfully removed the NAs!
        }
        if(sum(is.na(v[,i])) != 0)
        {
          nbhd = nbhd*2 # Double the neihhbourhood size to search
        }
      }
      count = count + 1
    }  
  }
  print(paste('How many NAs were found and removed',sum(NA_ind), 'out of',prod(dim(NA_ind))))
  print(paste('Neighbourhood size needed to remove NAs',nbhd))
  return(v)
} 

######## Interpolator function needed to compute the area-average of the gridded covariate over the dual mesh
f.interp.dmesh <- function(dual_mesh_poly,sp_grid,dmesh_dist = NULL) {
  # Step 1 - project the dual mesh polygons onto the CRS of the covariate grid
  dmesh_trans = spTransform(dual_mesh_poly, CRSobj = CRS(proj4string(sp_grid)))
  # Step 2 - aggregate the covariate values onto the dmesh by computing weighted average - removing NA covar cells.
  dmesh_covar = aggregate(sp_grid, by = dmesh_trans, FUN=mean, areaWeighted = T)
  dmesh_covar_robust = aggregate(sp_grid, by = dmesh_trans, FUN=median, na.rm=T)
  # Step 3 - Do we have NAs in the dual mesh?
  NA_ind = as.matrix(is.na(dmesh_covar@data))
  if(sum(NA_ind) > 0)
  {
    #print('NAs found')
    # Step 4 - If dmesh_dist matrix not supplied, compute it
    if (is.null(dmesh_dist))
    {
      dmesh_dist = gDistance(dual_mesh_poly,byid = T)
    }
    # Step 5 - loop over the different variables in sp_grid and match missing values to closest non-missing polygon value.
    missing_vars = which(colSums(NA_ind)>0)
    count=1 # progress report
    for (i in missing_vars)
    {
      print(paste(count,'out of',length(missing_vars),'missing variables processed')) 
      # find the closest non-missing polygon value and use it
      closest_ind = apply(dmesh_dist[NA_ind[,i],!NA_ind[,i]], 1, FUN = function(x){which(x == min(x))[1]})
      dmesh_covar@data[which(NA_ind[,i]),i] = dmesh_covar@data[which(!NA_ind[,i])[closest_ind],i]
      dmesh_covar_robust@data[which(NA_ind[,i]),i] = dmesh_covar_robust@data[which(!NA_ind[,i])[closest_ind],i]
      count = count + 1
    }
    
  }
  print(paste('How many NAs were found and removed',sum(NA_ind), 'out of',prod(dim(NA_ind))))
  return(list(weighted_mean = dmesh_covar, median = dmesh_covar_robust))
}

# create the dual mesh
# Load the Simpson dual mesh creator
inla.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}