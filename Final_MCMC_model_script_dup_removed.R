# MCMC model fitter
tic = proc.time()
no_samples = 1000
no_T = 6
seeds = sample.int(1e6, no_samples)
.libPaths('/zfs/users/joe.watson/joe.watson/rlibs')

# library.dynam('libgeos-3.7.0.so', rgeos, 
#               lib.loc='/zfs/users/joe.watson/joe.watson/opt/geos/lib',
#               verbose = getOption("verbose"),
#               file.ext = .Platform$dynlib.ext, ...)

# load the libraries
library(INLA)
#library(inlabru)
library(ggplot2)
library(sp)
library(rgeos)
#library(rgdal)
library(abind)
inla.setOption( num.threads = 10 ) 

########## ======== Load pre-compiled objects
MCMC_files = readRDS('MCMC_files_dup_removed.rds')
list2env(MCMC_files, .GlobalEnv)
rm(MCMC_files)

# read the pixels #
pixels_file = readRDS('pixels_meshpoly.rds')

# load the mesh
mesh = readRDS('mesh_transformed.rds')

# read the proj matrix
proj = readRDS('proj_meshpoly.rds')

# load the COAST_SIMP shapefile. Note we want a pretty polygon. We don't care if it 
# Does not match with the one used for computation.
# Not actually used!
COAST_simp = readRDS('COAST_transformed.rds')

# load the covariates data
covariates_pp = readRDS('covariates_pp_dup_removed.rds')

# load the model
simple.model = readRDS('simple.model.rds')

# load sampling intensity files
SI_WW = readRDS('SI_WW2_dup_removed.rds')
#SI_WW2 = SI_WW#readRDS('SI_WW2.rds')
#dmesh_barrier_transformed = readRDS('dmesh_barrier_transformed.rds')
dmesh = readRDS('dmesh_transformed.rds')

dmesh_dist = readRDS('dmesh_dist.rds')

E_J = array(0, dim = c(no_samples,mesh$n,no_T))
E_K = array(0, dim = c(no_samples,mesh$n,no_T))
E_L = array(0, dim = c(no_samples,mesh$n,no_T))

DIC_values = rep(0, times = no_samples)
####### ======== load functions

# rescale to equal the mean
rescale_fun = function(input, scale)
{
  if(length(scale)>1){# we should have a matrix input
    if(dim(input)[2] != length(scale))
    {
      stop('The length of the scale argument should match the number of cols of input')
    }
    for(i in 1:length(scale))
    {
      input[,i] = (input[,i] * scale[i]) / sum(input[,i], na.rm=T)
    }
  }
  if(length(scale)==1){# we should have a vector input
    if(!is.vector(input))
    {
      stop('The input should be a vector')
    }
    input = (input * scale) / sum(input, na.rm=T)}
  return(input)
}

######## Interpolator function needed to compute the area-average of the gridded covariate over the dual mesh
f.interp.dmesh <- function(dual_mesh_poly,sp_grid,dmesh_dist = NULL) {
  # Step 1 - project the dual mesh polygons onto the CRS of the covariate grid
  #dmesh_trans = spTransform(dual_mesh_poly, CRSobj = CRS(proj4string(sp_grid)))
  if( proj4string(sp_grid) != proj4string(dual_mesh_poly) )
  {
    stop('CRS not the same')
  }
  # Step 2 - aggregate the covariate values onto the dmesh by computing weighted average - removing NA covar cells.
  dmesh_covar = aggregate(sp_grid, by = dual_mesh_poly, FUN=mean, areaWeighted = T)
  dmesh_covar_robust = aggregate(sp_grid, by = dual_mesh_poly, FUN=median, na.rm=T)
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
  #print(paste('How many NAs were found and removed',sum(NA_ind), 'out of',prod(dim(NA_ind))))
  return(list(weighted_mean = dmesh_covar, median = dmesh_covar_robust))
}

# Create objects to store results

effort_layers = SI_WW$normalized_effort_layers

# standardise to sum to 1
for( i in 1:length(effort_layers))
{
  effort_layers[[i]]$field = (effort_layers[[i]]$field / 
                                sum(effort_layers[[i]]$field))
}



# run the simulations
for(i in 1:no_samples)
{

####### Step 1 load the search effort

  grid_length = length(SI_WW$normalized_effort_layers[[1]]@data$field)
  
  SI_temp = SI_WW$normalized_effort_layers[[1]]
  temp_df2 = data.frame(matrix(0, nrow = length(effort_layers[[1]]$field), ncol = 3*6))
  
  # compute the total boat hours (sum the months)
  boat_hours_monthly_port_J = SI_WW$total_boat_hours_per_period_per_port_J[[i]] 
  boat_hours_monthly_port_K = SI_WW$total_boat_hours_per_period_per_port_K[[i]] 
  boat_hours_monthly_port_L = SI_WW$total_boat_hours_per_period_per_port_L[[i]] 
  
  #SI_temp2@data = temp_df1
  
  SI_temp@data = temp_df2
  # for each port field normalize to sum to zero, then scale by total hours
  count = 1
  for(month_ind in 1:6)
  {
    
    for(port_ind in 1:16)
    {
      field_temp = effort_layers[[port_ind]]$field 
      
      SI_temp@data[,count] = SI_temp@data[,count] + field_temp *
        boat_hours_monthly_port_J[port_ind, month_ind]
      
      SI_temp@data[,(count+1)] = SI_temp@data[,(count+1)] + field_temp *
        boat_hours_monthly_port_K[port_ind, month_ind]
      
      SI_temp@data[,(count+2)] = SI_temp@data[,(count+2)] + field_temp *
        boat_hours_monthly_port_L[port_ind, month_ind]
      
    }
    count = count + 3 #6
  }

  # Map to dual mesh pixels - computing weighted average
  # Note that areaWeighted = TRUE always computes weighted average. Need to multiply by area of dmesh
  search_effort_dmesh = aggregate( SI_temp, dmesh, mean, areaWeighted=TRUE)@data

  pid <- sapply(slot(dmesh, "polygons"), function(x) slot(x, "ID")) 
  SI_dmesh = SpatialPolygonsDataFrame(dmesh,
                                      data.frame(search_effort_dmesh,#, SI_WW_complete_monthly_L_var),
                                                 row.names = pid))
  # fill in any gaps for completeness
  SI_dmesh = f.interp.dmesh(dmesh, SI_dmesh,dmesh_dist)$weighted_mean@data
  # give the columns names
  names_SI_dmesh = c(rep('SI_WW_complete_monthly_J',6), 
                     #rep('SI_WW_complete_monthly_J_var',6),
                     rep('SI_WW_complete_monthly_K',6), 
                     #rep('SI_WW_complete_monthly_K_var',6),
                     rep('SI_WW_complete_monthly_L',6))#,
  
  # scale the search effort density by the area of the pixels
  SI_dmesh = apply(SI_dmesh, 2 ,FUN=function(x){x*(w/mean(w))})
  
  # put the objects in the environment
  for(i2 in unique(names_SI_dmesh))
  {
    ind_name = which(names_SI_dmesh == i2)
    assign(i2, SI_dmesh[,ind_name])
  }
  rm(SI_dmesh)

  # remove effort on land
  SI_WW_complete_monthly_J[w==0,] = 0
  SI_WW_complete_monthly_K[w==0,] = 0
  SI_WW_complete_monthly_L[w==0,] = 0
  
  # rescale to equal total boat hours
  SI_WW_complete_monthly_J = rescale_fun(SI_WW_complete_monthly_J, colSums(boat_hours_monthly_port_J))
  SI_WW_complete_monthly_K = rescale_fun(SI_WW_complete_monthly_K, colSums(boat_hours_monthly_port_K))
  SI_WW_complete_monthly_L = rescale_fun(SI_WW_complete_monthly_L, colSums(boat_hours_monthly_port_L))
  
  # combine with DFO's search effort
  SI_complete_monthly_J = SI_WW_complete_monthly_J + SI_Brian_monthly_corrected_J
  SI_complete_monthly_K = SI_WW_complete_monthly_K + SI_Brian_monthly_corrected_K
  SI_complete_monthly_L = SI_WW_complete_monthly_L + SI_Brian_monthly_corrected_L
  
  ##### Step 2 Modify the INLA stack files
  covariates_pp_J$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_J))
  covariates_pp_K$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_K))
  covariates_pp_L$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_L))
  
  temp_names = intersect(names(covariates_pp_J), names(covariates_site_J))
  covariates_PP_join_J = rbind.data.frame(covariates_pp_J[,temp_names], 
                                          covariates_site_J[,temp_names])
  covariates_PP_join_J$POD_INLA = 1
  
  covariates_PP_join_K = rbind.data.frame(covariates_pp_K[,temp_names], 
                                          covariates_site_K[,temp_names])
  covariates_PP_join_K$POD_INLA = 2
  
  covariates_PP_join_L = rbind.data.frame(covariates_pp_L[,temp_names], 
                                          covariates_site_L[,temp_names])
  covariates_PP_join_L$POD_INLA = 3
  
  covariates_PP_join_J$Intercept = 1
  PP_stack_J = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                          A=list(1,A.pp_J), tag='pp',
                          effects=list(covariates_PP_join_J, 
                                       s_index)) # i are the indices for the SPDE 
  print('stack pp J complete')
  
  covariates_PP_join_K$Intercept = 1
  PP_stack_K = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                          A=list(1,A.pp_K), tag='pp',
                          effects=list(covariates_PP_join_K, 
                                       s_index)) # i are the indices for the SPDE 
  print('stack pp K complete')
  
  covariates_PP_join_L$Intercept = 1
  PP_stack_L = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                          A=list(1,A.pp_L, A.pp_L), tag='pp',
                          effects=list(covariates_PP_join_L, 
                                       s_index, s_index_L)) # i are the indices for the SPDE 
  print('stack pp L complete')
  
  # finally store the search effort for posterior predictive checks
  E_J[i,,] = matrix(covariates_pp_J$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F)
  E_K[i,,] = matrix(covariates_pp_K$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F)
  E_L[i,,] = matrix(covariates_pp_L$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F)
  
  Big_Stack_podSPDE = inla.stack(PP_stack_J, PP_stack_K, PP_stack_L)
  
  ###### Step 3 - Fit the Model
  simple_fit_covars4 = inla(y ~ -1 + f(simple.field, model = simple.model) + # shared by all pods
                                              #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                              f(simple.field.L, model = simple.model) + # L contrast
                                              factor(POD_INLA) +
                                              SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                              SSTminusspacetime + chlorominusspacetime + #depth +  
                                              f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                            data = inla.stack.data(Big_Stack_podSPDE),
                                            family = 'poisson',
                                            E = e,
                                            control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE)),
                                            control.compute = list(config=T, dic = T, cpo = F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                            control.mode = list(theta = c(3.66723160,  1.06620230,  5.32711111, -0.03530855,  1.25606437), restart=TRUE ), 
                                            control.inla = list(int.strategy = 'eb'),
                                            #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                            verbose = F)
  print(paste('Model number ',i,' out of  ',no_samples,' fitted',sep = ''))
  DIC_values[i] = simple_fit_covars4$dic$dic
  ####### Step 4 - Sample once from the posterior
  
  if(i==1)
  {
    temp_mode = simple_fit_covars4$misc$theta.mode
    saveRDS(temp_mode,'temp_mode.rds')
    samp = inla.posterior.sample(n = 1, simple_fit_covars4, seed = seeds[i]) 
    # reduce the size of samp by removing redundant information
    samp[[i]]$latent = samp[[i]]$latent[-c(grep('Predictor', rownames(samp[[i]]$latent), fixed = T)),]
  }
  if(i!=1)
  {
    samp[i] = inla.posterior.sample(n = 1, simple_fit_covars4, seed = seeds[i])
    # reduce the size of samp by removing redundant information
    samp[[i]]$latent = samp[[i]]$latent[-c(grep('Predictor', rownames(samp[[i]]$latent), fixed = T)),]
  }
  rm(simple_fit_covars4)
  print('Successfully sampled from posterior')
  saveRDS(samp,'MCMC_samples_finalmodel_dup_removed.rds')
  saveRDS(list(E_J = E_J, E_K = E_K, E_L = E_L), 'E_MCMC_dup_removed.rds')
  saveRDS(seeds,'seeds_dup_removed.rds')
  saveRDS(DIC_values, 'DIC_values_dup_removed.rds')
  print('Successfully saved intermediate results')
}
  
  print('finished model fitting')
  
  ###### Step 5 - Create the results file from the posterior results
  # define own plotting function
  local.plot.field = function(mesh, samp, E_J, E_K, E_L, poly, n_samples, nx = 150, ny = 150, 
                              no_T=1, funct = exp, mesh_data = NULL,
                              model = '0', pixels_file, ...){
    # define spatial pixels to plot over
    #pixels = pixels(mesh, mask = poly, nx = nx, ny = ny)
    
    # generate n_samples MCMC from model
    #samp = inla.posterior.sample(n = n_samples, fit, seed = 1234)
    
    # store hyperparameters and parameters (fixed effects)
    hyper_pars = matrix(0, nrow = n_samples, ncol = length(samp[[1]]$hyperpar))
    colnames(hyper_pars) = names(samp[[1]]$hyperpar)
    pars = matrix(0, nrow = n_samples, ncol = 9)
    names_pars = c("factor(POD_INLA)1","factor(POD_INLA)2","factor(POD_INLA)3",   
    "SSTmonthavg","SSTspatialavg","chloromonthavg","chlorospatialavg",
    "SSTminusspacetime","chlorominusspacetime")
    colnames(pars) = names_pars
    
    MONTH_INLA_J = matrix(0, nrow = n_samples, ncol = no_T)
    MONTH_INLA_K = matrix(0, nrow = n_samples, ncol = no_T)
    MONTH_INLA_L = matrix(0, nrow = n_samples, ncol = no_T)
    
    # create arrays to store results of particular covariates etc...
    mesh_field_J = array(0, dim = c(n_samples, no_T, mesh$n))
    mesh_field_K = array(0, dim = c(n_samples, no_T, mesh$n))
    mesh_field_L = array(0, dim = c(n_samples, no_T, mesh$n))
    
    # combined field (sum of intensities)
    mesh_field = array(0, dim = c(n_samples, no_T, mesh$n))
    
    if(model != '0') # covariates are included
    {
      mesh_covariates = array(0, dim = c(n_samples, no_T, mesh$n))
      
      pixels_covariates = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
      
      predicted_sightings = array(0, dim = c(n_samples, no_T, 3))
      
      pod_probabilities_J = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
      pod_probabilities_K = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
      pod_probabilities_L = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
      
    }
    
    pixels_field = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
    pixels_field_J = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
    pixels_field_K = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
    pixels_field_L = array(0, dim = c(n_samples, no_T, dim(pixels_file@coords)[1]))
    
    # reduce the size of samp by removing the reduntant information #
    #for (i in 1:n_samples)
    #{
    #  samp[[i]]$latent = samp[[i]]$latent[-c(grep('Predictor', rownames(samp[[i]]$latent), fixed = T)),]
    #}
    
    # - choose plotting region to be the same as the study area polygon
    #proj = inla.spde.make.A(mesh, loc = coordinates(pixels_file))
    ##### MAY NEED TO BE NAMES NOT ROWNAMES CHECK IF FAILS!
    spatial_field_ind = which(startsWith(prefix = 'simple.field:',x = names(samp[[1]]$latent)))
    spatial_field_ind_L = which(startsWith(prefix = 'simple.field.L:',x = names(samp[[1]]$latent)))
    MONTH_INLA_ind = which(startsWith(prefix = 'MONTH_INLA',x = names(samp[[1]]$latent)))
    
    print('Starting results compilation')
    
    for(i in 1:n_samples)
    {
      counter = 1
      
      pars[i,] = as.numeric(samp[[i]]$latent[names_pars])
      hyper_pars[i,] = as.numeric(samp[[i]]$hyperpar)
      MONTH_INLA_J[i,] = as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[1:no_T]
      MONTH_INLA_K[i,] = as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[no_T+c(1:no_T)]
      MONTH_INLA_L[i,] = as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[(no_T*2)+c(1:no_T)]
      
      print('stored parameter estimates')
      
      for(j in 1:no_T) # loop over the months
      {
        ind = (mesh$n * (j-1) + 1) : (mesh$n * j)
        
        mesh_field_J[i,j,] = #samp[[i]]$latent['Intercept'] +
          as.numeric(samp[[i]]$latent[spatial_field_ind]) +
          as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j] +#[ind]  # Spatiotemporal field 
          as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)1',x = names(samp[[i]]$latent)))])
        
        mesh_field_K[i,j,] = #samp[[i]]$latent['Intercept'] +
          as.numeric(samp[[i]]$latent[spatial_field_ind]) +
          as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+6] +#[ind]  # Spatiotemporal field
          as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)2',x = names(samp[[i]]$latent)))])
        
        mesh_field_L[i,j,] = #samp[[i]]$latent['Intercept'] +
          as.numeric(samp[[i]]$latent[spatial_field_ind]) +
          as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+12] +#[ind]  # Spatiotemporal field
          as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)3',x = names(samp[[i]]$latent)))]) +
          as.numeric(samp[[i]]$latent[spatial_field_ind_L])
        
        print('temp 1')
        
        mesh_field[i,j,] = mesh_field_J[i,j,] + mesh_field_K[i,j,] + mesh_field_L[i,j,]
        
        # project onto pixels
        pixels_field_J[i,j,] = as.numeric(proj %*% mesh_field_J[i,j,])
        pixels_field_K[i,j,] = as.numeric(proj %*% mesh_field_K[i,j,])
        pixels_field_L[i,j,] = as.numeric(proj %*% mesh_field_L[i,j,])
        pixels_field[i,j,] = as.numeric(proj %*% mesh_field[i,j,])
        
        print('temp 2')
        
        if(model=='1')
        {
          #ind = (mesh$n * (j-1) + 1) : (mesh$n * j)
          mesh_covariates[i,j,] = samp[[i]]$latent['Intercept'] +
            samp[[i]]$latent['SSTmonthavg'] * mesh_data[ind,'SSTmonthavg'] +
            samp[[i]]$latent['SSTminusmonth'] * mesh_data[ind,'SSTminusmonth'] +
            samp[[i]]$latent['chloromonthavg'] * mesh_data[ind,'chloromonthavg'] +
            samp[[i]]$latent['chlorominusmonth'] * mesh_data[ind,'chlorominusmonth'] +
            samp[[i]]$latent['vessel'] * mesh_data[ind,'vessel'] 
          #samp[[i]]$latent['depth'] * mesh_data[ind,'depth'] 
          
          pixels_covariates[i,j,] = as.numeric(proj %*% mesh_covariates[i,j,])
          
        }
        if(model=='2')
        {
          #ind = (mesh$n * (j-1) + 1) : (mesh$n * j)
          mesh_covariates[i,j,] = samp[[i]]$latent['Intercept'] +
            samp[[i]]$latent['SSTmonthavg'] * mesh_data[ind,'SSTmonthavg'] +
            samp[[i]]$latent['SSTspatialavg'] * mesh_data[ind,'SSTspatialavg'] +
            samp[[i]]$latent['SSTminusspacetime'] * mesh_data[ind,'SSTminusspacetime'] +
            samp[[i]]$latent['chloromonthavg'] * mesh_data[ind,'chloromonthavg'] +
            samp[[i]]$latent['chlorospatialavg'] * mesh_data[ind,'chlorospatialavg'] +
            samp[[i]]$latent['chlorominusspacetime'] * mesh_data[ind,'chlorominusspacetime'] +
            samp[[i]]$latent['vessel'] * mesh_data[ind,'vessel'] 
          #+samp[[i]]$latent['depth'] * mesh_data[ind,'depth'] 
          
          pixels_covariates[i,j,] = as.numeric(proj %*% mesh_covariates[i,j,])
          
        }
        if(model=='3')
        {
          #ind = (mesh$n * (j-1) + 1) : (mesh$n * j)
          mesh_covariates[i,j,] = #samp[[i]]$latent['Intercept'] +
            samp[[i]]$latent['SSTmonthavg'] * mesh_data[ind,'SSTmonthavg'] +
            samp[[i]]$latent['SSTspatialavg'] * mesh_data[ind,'SSTspatialavg'] +
            samp[[i]]$latent['SSTminusspacetime'] * mesh_data[ind,'SSTminusspacetime'] +
            samp[[i]]$latent['chloromonthavg'] * mesh_data[ind,'chloromonthavg'] +
            samp[[i]]$latent['chlorospatialavg'] * mesh_data[ind,'chlorospatialavg'] +
            samp[[i]]$latent['chlorominusspacetime'] * mesh_data[ind,'chlorominusspacetime'] #+
          #samp[[i]]$latent['depth'] * mesh_data[ind,'depth'] 
          #samp[[i]]$latent['vessel'] * mesh_data[ind,'vessel'] 
          
          print('temp 3')
          pixels_covariates[i,j,] = as.numeric(proj %*% mesh_covariates[i,j,])
          
        }
        # predict how many SRKW sightings we expect from each pod
        predicted_sightings[i,j,1] = sum(rpois(mesh$n, exp((mesh_covariates[i,j,] + mesh_field_J[i,j,]))*E_J[i,,j]))
        predicted_sightings[i,j,2] = sum(rpois(mesh$n, exp((mesh_covariates[i,j,] + mesh_field_K[i,j,]))*E_K[i,,j]))
        predicted_sightings[i,j,3] = sum(rpois(mesh$n, exp((mesh_covariates[i,j,] + mesh_field_L[i,j,]))*E_L[i,,j]))
        
        print('temp 4')
        # save files to use for computing pod-specific probabilities
        pod_probabilities_J[i,j,] = exp( (as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j] +
                                            as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)1',x = names(samp[[i]]$latent)))]) ) -
                                           log(  (    exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j] +
                                                             as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)1',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+6] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)2',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+12] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)3',x = names(samp[[i]]$latent)))]) +
                                                               as.numeric((proj %*% as.numeric(samp[[i]]$latent[spatial_field_ind_L]))) ) ) ) )
        
        
        pod_probabilities_K[i,j,] = exp( (as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+6] +
                                            as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)2',x = names(samp[[i]]$latent)))]) ) -
                                           log(  (    exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j] +
                                                             as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)1',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+6] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)2',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+12] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)3',x = names(samp[[i]]$latent)))]) +
                                                               as.numeric((proj %*% as.numeric(samp[[i]]$latent[spatial_field_ind_L]))) ) ) ) )
        
        pod_probabilities_L[i,j,] = exp( (as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+12] +
                                            as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)3',x = names(samp[[i]]$latent)))]) +
                                            as.numeric((proj %*% as.numeric(samp[[i]]$latent[spatial_field_ind_L]) ))) -
                                           log(  (    exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j] +
                                                             as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)1',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+6] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)2',x = names(samp[[i]]$latent)))])) +
                                                        exp( as.numeric(samp[[i]]$latent[MONTH_INLA_ind])[j+12] +
                                                               as.numeric(samp[[i]]$latent[which(startsWith(prefix = 'factor(POD_INLA)3',x = names(samp[[i]]$latent)))]) +
                                                               as.numeric((proj %*% as.numeric(samp[[i]]$latent[spatial_field_ind_L]))) ) ) ) )
      }
      samp[[i]] <- 1 # reduce size in memory
    }
    
    print('computed linear predictors and pod probabilities')
    
    if(model != '0')
    {
      # Compute variances and covariances between covariates and GMRF
      covariates_mean = mean(pixels_covariates)
      pixels_covariates_spatial_variance = apply((pixels_covariates-covariates_mean)^2, 3, mean)
      covariates_variance = mean((pixels_covariates-covariates_mean)^2)
      
      field_mean = mean(pixels_field)
      pixels_field_spatial_variance = apply((pixels_field-field_mean)^2, 3, mean)
      field_variance = mean((pixels_field-field_mean)^2)
      
      # save the GMRF field for later plotting
      Post_mean_GMRF = t(apply(pixels_field, c(2,3), function(x){mean(x, na.rm=T)}))
      Post_SD_GMRF = t(apply(pixels_field, c(2,3), function(x){sd(x, na.rm=T)}))
      
      # covariance + correlation between covariates and field
      field_covariates_covariance = cov(as.numeric(pixels_field), as.numeric(pixels_covariates))
      field_covariates_correlation = cor(as.numeric(pixels_field), as.numeric(pixels_covariates))
      
      pixels_field = pixels_field + 3*pixels_covariates # add the covariates and the field to form LP
      pixels_field_J = pixels_field_J + pixels_covariates
      pixels_field_K = pixels_field_K + pixels_covariates
      pixels_field_L = pixels_field_L + pixels_covariates
      
      linearpred_mean = mean(pixels_field)
      linearpred_deciles = quantile(pixels_field, probs = seq(0.5, 1, 0.1), na.rm = T)
      linearpred_deciles_monthly = apply(pixels_field, 2, quantile, probs = seq(0.5, 1, 0.1), na.rm=T)
      linearpred_deciles_low = quantile(pixels_field, probs = seq(0.1, 0.5, 0.1), na.rm = T)
      linearpred_deciles_low_monthly = apply(pixels_field, 2, quantile, probs = seq(0.1, 0.5, 0.1), na.rm=T)
      # the monthly one are n_probs x no_T in dimension
      
      linearpred_deciles_pod = quantile(c(pixels_field_J,pixels_field_K,pixels_field_L), probs = seq(0.5, 1, 0.1), na.rm = T)
      linearpred_deciles_pod_monthly = apply(abind(pixels_field_J,pixels_field_K,pixels_field_L,along=1), 2, quantile, probs = seq(0.5, 1, 0.1), na.rm=T)
      linearpred_deciles_pod_low = quantile(c(pixels_field_J,pixels_field_K,pixels_field_L), probs = seq(0.1, 0.5, 0.1), na.rm = T)
      linearpred_deciles_pod_low_monthly = apply(abind(pixels_field_J,pixels_field_K,pixels_field_L,along=1), 2, quantile, probs = seq(0.1, 0.5, 0.1), na.rm=T)
      # this is n_probs x no_T in dimension
      
      pixels_linearpred_spatial_variance = apply((pixels_field-linearpred_mean)^2, 3, mean)
      linearpred_variance = mean((pixels_field - linearpred_mean)^2)
      
      pixels_covariates_field_spatial_covariance = apply((pixels_covariates-covariates_mean)*(pixels_field-field_mean), 3, mean)
      
      
      # Summary statistics on posterior predictive estimated number of sightings
      predicted_sightings_LCL = t(apply(predicted_sightings, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))}))  
      predicted_sightings_UCL = t(apply(predicted_sightings, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))  
      predicted_sightings_Median = t(apply(predicted_sightings, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.5))}))  
    }
    
    # compute median, quantiles etc of results
    Post_mean = t(apply(pixels_field, c(2,3), function(x){mean(x, na.rm=T)}))
    Post_LCL = apply((apply(pixels_field, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
    Post_UCL = apply((apply(pixels_field, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
    Post_SD = t(apply(pixels_field, c(2,3), function(x){sd(x, na.rm=T)}))
    Post_LCL2 = t(apply(pixels_field, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))})) 
    Post_UCL2 = t(apply(pixels_field, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))
    Post_Bigger_D5 = t(apply(pixels_field, c(2,3), function(x){mean(x > linearpred_deciles[1], na.rm=T)}))
    Post_Bigger_D6 = t(apply(pixels_field, c(2,3), function(x){mean(x > linearpred_deciles[2], na.rm=T)}))
    Post_Bigger_D7 = t(apply(pixels_field, c(2,3), function(x){mean(x > linearpred_deciles[3], na.rm=T)}))
    Post_Bigger_D8 = t(apply(pixels_field, c(2,3), function(x){mean(x > linearpred_deciles[4], na.rm=T)}))
    Post_Bigger_D9 = t(apply(pixels_field, c(2,3), function(x){mean(x > linearpred_deciles[5], na.rm=T)}))
    Post_Smaller_D1 = t(apply(pixels_field, c(2,3), function(x){mean(x < linearpred_deciles_low[1], na.rm=T)}))
    Post_Smaller_D2 = t(apply(pixels_field, c(2,3), function(x){mean(x < linearpred_deciles_low[2], na.rm=T)}))
    Post_Smaller_D3 = t(apply(pixels_field, c(2,3), function(x){mean(x < linearpred_deciles_low[3], na.rm=T)}))
    Post_Smaller_D4 = t(apply(pixels_field, c(2,3), function(x){mean(x < linearpred_deciles_low[4], na.rm=T)}))
    Post_Smaller_D5 = t(apply(pixels_field, c(2,3), function(x){mean(x < linearpred_deciles_low[5], na.rm=T)}))
    
    # compute the exceedance maps per month
    for(i in 1:no_T)
    {
      if(i==1)
      {
        print('Creating the first monthly exceedance maps')
        Post_Bigger_D5_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[1,i], na.rm=T)})
        Post_Bigger_D6_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[2,i], na.rm=T)})
        Post_Bigger_D7_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[3,i], na.rm=T)})
        Post_Bigger_D8_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[4,i], na.rm=T)})
        Post_Bigger_D9_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[5,i], na.rm=T)})
        Post_Smaller_D1_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[1,i], na.rm=T)})
        Post_Smaller_D2_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[2,i], na.rm=T)})
        Post_Smaller_D3_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[3,i], na.rm=T)})
        Post_Smaller_D4_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[4,i], na.rm=T)})
        Post_Smaller_D5_Monthly = apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[5,i], na.rm=T)})
        # returns a vector of length dim(pixels)[1]
      }
      if(i!=1)
      {
        print('binding the monthly exceedance maps')
        Post_Bigger_D5_Monthly = cbind(Post_Bigger_D5_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[1,i], na.rm=T)}))
        Post_Bigger_D6_Monthly = cbind(Post_Bigger_D6_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[2,i], na.rm=T)}))
        Post_Bigger_D7_Monthly = cbind(Post_Bigger_D7_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[3,i], na.rm=T)}))
        Post_Bigger_D8_Monthly = cbind(Post_Bigger_D8_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[4,i], na.rm=T)}))
        Post_Bigger_D9_Monthly = cbind(Post_Bigger_D9_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x > linearpred_deciles_monthly[5,i], na.rm=T)}))
        Post_Smaller_D1_Monthly = cbind(Post_Smaller_D1_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[1,i], na.rm=T)}))
        Post_Smaller_D2_Monthly = cbind(Post_Smaller_D2_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[2,i], na.rm=T)}))
        Post_Smaller_D3_Monthly = cbind(Post_Smaller_D3_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[3,i], na.rm=T)}))
        Post_Smaller_D4_Monthly = cbind(Post_Smaller_D4_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[4,i], na.rm=T)}))
        Post_Smaller_D5_Monthly = cbind(Post_Smaller_D5_Monthly, apply(pixels_field[,i,], c(2), function(x){mean(x < linearpred_deciles_low_monthly[5,i], na.rm=T)}))
        # returns a matrix of dim dim(pixels)[1] x no_T
      }
    }
    
    funct_pixels_field = funct(pixels_field)
    linearpred_funct_deciles = quantile(funct_pixels_field, probs = seq(0.5, 1, 0.1), na.rm = T)
    
    # compute median, quantiles etc of results
    Post_funct_median = t(apply(funct_pixels_field, c(2,3), function(x){median(x, na.rm=T)}))
    #Post_funct_LCL = apply((apply(funct_pixels_field, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
    #Post_funct_UCL = apply((apply(funct_pixels_field, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
    Post_funct_MAD = t(apply(funct_pixels_field, c(2,3), function(x){mean(abs(x-median(x)), na.rm=T)}))
    #Post_funct_LCL2 = t(apply(funct_pixels_field, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))})) 
    #Post_funct_UCL2 = t(apply(funct_pixels_field, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))
    #Post_funct_Bigger_D5 = t(apply(funct_pixels_field, c(2,3), function(x){mean(x > linearpred_funct_deciles[1], na.rm=T)}))
    
    ################################
    # STORE THE POD-SPECIFIC RESULTS
    ################################
    # compute median, quantiles etc of results
    print('save results for each pod')
    
    Post_mean_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x, na.rm=T)}))
    Post_LCL_J = apply((apply(pixels_field_J, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
    Post_UCL_J = apply((apply(pixels_field_J, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
    Post_SD_J = t(apply(pixels_field_J, c(2,3), function(x){sd(x, na.rm=T)}))
    Post_LCL2_J = t(apply(pixels_field_J, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))})) 
    Post_UCL2_J = t(apply(pixels_field_J, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))
    Post_Bigger_D5_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x > linearpred_deciles_pod[1], na.rm=T)}))
    Post_Bigger_D6_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x > linearpred_deciles_pod[2], na.rm=T)}))
    Post_Bigger_D7_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x > linearpred_deciles_pod[3], na.rm=T)}))
    Post_Bigger_D8_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x > linearpred_deciles_pod[4], na.rm=T)}))
    Post_Bigger_D9_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x > linearpred_deciles_pod[5], na.rm=T)}))
    Post_Smaller_D1_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[1], na.rm=T)}))
    Post_Smaller_D2_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[2], na.rm=T)}))
    Post_Smaller_D3_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[3], na.rm=T)}))
    Post_Smaller_D4_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[4], na.rm=T)}))
    Post_Smaller_D5_J = t(apply(pixels_field_J, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[5], na.rm=T)}))
    Post_Prob_J = apply(pod_probabilities_J, c(2,3), function(x){median(x, na.rm=T)})
    
    Post_mean_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x, na.rm=T)}))
    Post_LCL_K = apply((apply(pixels_field_K, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
    Post_UCL_K = apply((apply(pixels_field_K, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
    Post_SD_K = t(apply(pixels_field_K, c(2,3), function(x){sd(x, na.rm=T)}))
    Post_LCL2_K = t(apply(pixels_field_K, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))})) 
    Post_UCL2_K = t(apply(pixels_field_K, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))
    Post_Bigger_D5_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x > linearpred_deciles_pod[1], na.rm=T)}))
    Post_Bigger_D6_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x > linearpred_deciles_pod[2], na.rm=T)}))
    Post_Bigger_D7_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x > linearpred_deciles_pod[3], na.rm=T)}))
    Post_Bigger_D8_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x > linearpred_deciles_pod[4], na.rm=T)}))
    Post_Bigger_D9_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x > linearpred_deciles_pod[5], na.rm=T)}))
    Post_Smaller_D1_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[1], na.rm=T)}))
    Post_Smaller_D2_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[2], na.rm=T)}))
    Post_Smaller_D3_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[3], na.rm=T)}))
    Post_Smaller_D4_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[4], na.rm=T)}))
    Post_Smaller_D5_K = t(apply(pixels_field_K, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[5], na.rm=T)}))
    Post_Prob_K = apply(pod_probabilities_K, c(2,3), function(x){median(x, na.rm=T)})
    
    Post_mean_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x, na.rm=T)}))
    Post_LCL_L = apply((apply(pixels_field_L, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.025))
    Post_UCL_L = apply((apply(pixels_field_L, c(1,2), function(x) {mean(x, na.rm=TRUE)})), 2, quantile, probs = c(0.975))
    Post_SD_L = t(apply(pixels_field_L, c(2,3), function(x){sd(x, na.rm=T)}))
    Post_LCL2_L = t(apply(pixels_field_L, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.025))})) 
    Post_UCL2_L = t(apply(pixels_field_L, c(2,3), function(x){quantile(x, na.rm=T, probs = c(0.975))}))
    Post_Bigger_D5_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x > linearpred_deciles_pod[1], na.rm=T)}))
    Post_Bigger_D6_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x > linearpred_deciles_pod[2], na.rm=T)}))
    Post_Bigger_D7_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x > linearpred_deciles_pod[3], na.rm=T)}))
    Post_Bigger_D8_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x > linearpred_deciles_pod[4], na.rm=T)}))
    Post_Bigger_D9_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x > linearpred_deciles_pod[5], na.rm=T)}))
    Post_Smaller_D1_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[1], na.rm=T)}))
    Post_Smaller_D2_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[2], na.rm=T)}))
    Post_Smaller_D3_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[3], na.rm=T)}))
    Post_Smaller_D4_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[4], na.rm=T)}))
    Post_Smaller_D5_L = t(apply(pixels_field_L, c(2,3), function(x){mean(x < linearpred_deciles_pod_low[5], na.rm=T)}))
    Post_Prob_L = apply(pod_probabilities_L, c(2,3), function(x){median(x, na.rm=T)})
    
    # compute the exceedance maps per month
    for(i in 1:no_T)
    {
      if(i==1)
      {
        print('Creating the first monthly exceedance maps')
        Post_Bigger_D5_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)})
        Post_Bigger_D6_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)})
        Post_Bigger_D7_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)})
        Post_Bigger_D8_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)})
        Post_Bigger_D9_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)})
        Post_Smaller_D1_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)})
        Post_Smaller_D2_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)})
        Post_Smaller_D3_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)})
        Post_Smaller_D4_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)})
        Post_Smaller_D5_Monthly_J = apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)})
        
        Post_Bigger_D5_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)})
        Post_Bigger_D6_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)})
        Post_Bigger_D7_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)})
        Post_Bigger_D8_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)})
        Post_Bigger_D9_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)})
        Post_Smaller_D1_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)})
        Post_Smaller_D2_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)})
        Post_Smaller_D3_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)})
        Post_Smaller_D4_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)})
        Post_Smaller_D5_Monthly_K = apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)})
        
        Post_Bigger_D5_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)})
        Post_Bigger_D6_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)})
        Post_Bigger_D7_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)})
        Post_Bigger_D8_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)})
        Post_Bigger_D9_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)})
        Post_Smaller_D1_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)})
        Post_Smaller_D2_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)})
        Post_Smaller_D3_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)})
        Post_Smaller_D4_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)})
        Post_Smaller_D5_Monthly_L = apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)})
        # returns a vector of length dim(pixels)[1]
      }
      if(i!=1)
      {
        print('binding the monthly exceedance maps')
        Post_Bigger_D5_Monthly_J = cbind(Post_Bigger_D5_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)}))
        Post_Bigger_D6_Monthly_J = cbind(Post_Bigger_D6_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)}))
        Post_Bigger_D7_Monthly_J = cbind(Post_Bigger_D7_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)}))
        Post_Bigger_D8_Monthly_J = cbind(Post_Bigger_D8_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)}))
        Post_Bigger_D9_Monthly_J = cbind(Post_Bigger_D9_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)}))
        Post_Smaller_D1_Monthly_J = cbind(Post_Smaller_D1_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)}))
        Post_Smaller_D2_Monthly_J = cbind(Post_Smaller_D2_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)}))
        Post_Smaller_D3_Monthly_J = cbind(Post_Smaller_D3_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)}))
        Post_Smaller_D4_Monthly_J = cbind(Post_Smaller_D4_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)}))
        Post_Smaller_D5_Monthly_J = cbind(Post_Smaller_D5_Monthly_J, apply(pixels_field_J[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)}))
        
        Post_Bigger_D5_Monthly_K = cbind(Post_Bigger_D5_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)}))
        Post_Bigger_D6_Monthly_K = cbind(Post_Bigger_D6_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)}))
        Post_Bigger_D7_Monthly_K = cbind(Post_Bigger_D7_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)}))
        Post_Bigger_D8_Monthly_K = cbind(Post_Bigger_D8_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)}))
        Post_Bigger_D9_Monthly_K = cbind(Post_Bigger_D9_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)}))
        Post_Smaller_D1_Monthly_K = cbind(Post_Smaller_D1_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)}))
        Post_Smaller_D2_Monthly_K = cbind(Post_Smaller_D2_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)}))
        Post_Smaller_D3_Monthly_K = cbind(Post_Smaller_D3_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)}))
        Post_Smaller_D4_Monthly_K = cbind(Post_Smaller_D4_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)}))
        Post_Smaller_D5_Monthly_K = cbind(Post_Smaller_D5_Monthly_K, apply(pixels_field_K[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)}))
        
        Post_Bigger_D5_Monthly_L = cbind(Post_Bigger_D5_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[1,i], na.rm=T)}))
        Post_Bigger_D6_Monthly_L = cbind(Post_Bigger_D6_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[2,i], na.rm=T)}))
        Post_Bigger_D7_Monthly_L = cbind(Post_Bigger_D7_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[3,i], na.rm=T)}))
        Post_Bigger_D8_Monthly_L = cbind(Post_Bigger_D8_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[4,i], na.rm=T)}))
        Post_Bigger_D9_Monthly_L = cbind(Post_Bigger_D9_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x > linearpred_deciles_pod_monthly[5,i], na.rm=T)}))
        Post_Smaller_D1_Monthly_L = cbind(Post_Smaller_D1_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[1,i], na.rm=T)}))
        Post_Smaller_D2_Monthly_L = cbind(Post_Smaller_D2_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[2,i], na.rm=T)}))
        Post_Smaller_D3_Monthly_L = cbind(Post_Smaller_D3_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[3,i], na.rm=T)}))
        Post_Smaller_D4_Monthly_L = cbind(Post_Smaller_D4_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[4,i], na.rm=T)}))
        Post_Smaller_D5_Monthly_L = cbind(Post_Smaller_D5_Monthly_L, apply(pixels_field_L[,i,], c(2), function(x){mean(x < linearpred_deciles_pod_low_monthly[5,i], na.rm=T)}))
        # returns a matrix of dim dim(pixels)[1] x no_T
      }
    }
    
    if(model == '0')
    {
      pixels2 = SpatialPixelsDataFrame(pixels_file, data = data.frame(Post_mean, Post_LCL, Post_UCL, Post_SD, Post_LCL2, Post_UCL2,
                                                                      Post_funct_mean, Post_funct_LCL, Post_funct_UCL, Post_funct_SD, Post_funct_LCL2, Post_funct_UCL2))
    }
    if(model != '0')
    {
      # pixels2 = SpatialPixelsDataFrame(pixels_file, data = data.frame(Post_mean, Post_LCL, Post_UCL, Post_SD, Post_LCL2, Post_UCL2,
      #                                                            Post_funct_mean, Post_funct_LCL, Post_funct_UCL, Post_funct_SD, Post_funct_LCL2, Post_funct_UCL2,
      #                                                            pixels_covariates_spatial_variance, pixels_field_spatial_variance,
      #                                                            pixels_linearpred_spatial_variance, covariates_variance, field_variance,
      #                                                            linearpred_variance, field_covariates_covariance, field_covariates_correlation,
      #                                                            Post_mean_GMRF, Post_SD_GMRF))
      pixels2 = list(pixels_file,Post_mean, Post_LCL, Post_UCL, Post_SD, Post_LCL2, Post_UCL2,
                     Post_funct_median, Post_funct_MAD,
                     pixels_covariates_spatial_variance, pixels_field_spatial_variance,
                     pixels_linearpred_spatial_variance, pixels_covariates_field_spatial_covariance, 
                     covariates_variance, field_variance,
                     linearpred_variance, field_covariates_covariance, field_covariates_correlation,
                     Post_mean_GMRF, Post_SD_GMRF,
                     Post_Bigger_D5,Post_Bigger_D6,Post_Bigger_D7,Post_Bigger_D8,Post_Bigger_D9,
                     Post_Smaller_D1,Post_Smaller_D2,Post_Smaller_D3,Post_Smaller_D4,Post_Smaller_D5,
                     Post_Bigger_D5_Monthly,Post_Bigger_D6_Monthly,Post_Bigger_D7_Monthly,Post_Bigger_D8_Monthly,
                     Post_Bigger_D9_Monthly,Post_Smaller_D1_Monthly,Post_Smaller_D2_Monthly,Post_Smaller_D3_Monthly,
                     Post_Smaller_D4_Monthly,Post_Smaller_D5_Monthly,
                     predicted_sightings_LCL,predicted_sightings_UCL,predicted_sightings_Median)
      names(pixels2) = c('pixels_file','Post_mean', 'Post_LCL', 'Post_UCL', 'Post_SD', 'Post_LCL2', 'Post_UCL2',
                         'Post_funct_median', 'Post_funct_MAD',
                         'pixels_covariates_spatial_variance', 'pixels_field_spatial_variance',
                         'pixels_linearpred_spatial_variance','pixels_covariates_field_spatial_covariance', 
                         'covariates_variance', 'field_variance',
                         'linearpred_variance', 'field_covariates_covariance', 'field_covariates_correlation',
                         'Post_mean_GMRF', 'Post_SD_GMRF', 'Post_Bigger_D5', 'Post_Bigger_D6', 'Post_Bigger_D7',
                         'Post_Bigger_D8', 'Post_Bigger_D9','Post_Smaller_D1', 'Post_Smaller_D2', 'Post_Smaller_D3',
                         'Post_Smaller_D4', 'Post_Smaller_D5','Post_Bigger_D5_Monthly','Post_Bigger_D6_Monthly',
                         'Post_Bigger_D7_Monthly','Post_Bigger_D8_Monthly',
                         'Post_Bigger_D9_Monthly','Post_Smaller_D1_Monthly','Post_Smaller_D2_Monthly',
                         'Post_Smaller_D3_Monthly','Post_Smaller_D4_Monthly','Post_Smaller_D5_Monthly',
                         'predicted_sightings_LCL','predicted_sightings_UCL',
                         'predicted_sightings_Median')
      
      #########
      # REPEAT FOR PODS
      #########
      pixels2_J = list(pixels_file, Post_mean_J, Post_LCL_J, Post_UCL_J, Post_SD_J, Post_LCL2_J, Post_UCL2_J,
                       Post_Bigger_D5_J,Post_Bigger_D6_J,Post_Bigger_D7_J,Post_Bigger_D8_J,Post_Bigger_D9_J,
                       Post_Smaller_D1_J,Post_Smaller_D2_J,Post_Smaller_D3_J,Post_Smaller_D4_J,Post_Smaller_D5_J,
                       Post_Bigger_D5_Monthly_J,Post_Bigger_D6_Monthly_J,Post_Bigger_D7_Monthly_J,Post_Bigger_D8_Monthly_J,
                       Post_Bigger_D9_Monthly_J,Post_Smaller_D1_Monthly_J,Post_Smaller_D2_Monthly_J,Post_Smaller_D3_Monthly_J,
                       Post_Smaller_D4_Monthly_J,Post_Smaller_D5_Monthly_J,Post_Prob_J)
      names(pixels2_J) = c('pixels_file','Post_mean', 'Post_LCL', 'Post_UCL', 'Post_SD', 'Post_LCL2', 'Post_UCL2',
                           'Post_Bigger_D5', 'Post_Bigger_D6', 'Post_Bigger_D7',
                           'Post_Bigger_D8', 'Post_Bigger_D9', 'Post_Smaller_D1',
                           'Post_Smaller_D2', 'Post_Smaller_D3', 'Post_Smaller_D4',
                           'Post_Smaller_D5','Post_Bigger_D5_Monthly','Post_Bigger_D6_Monthly',
                           'Post_Bigger_D7_Monthly','Post_Bigger_D8_Monthly',
                           'Post_Bigger_D9_Monthly','Post_Smaller_D1_Monthly','Post_Smaller_D2_Monthly',
                           'Post_Smaller_D3_Monthly','Post_Smaller_D4_Monthly','Post_Smaller_D5_Monthly',
                           'Post_Prob_J')
      pixels2_K = list(pixels_file, Post_mean_K, Post_LCL_K, Post_UCL_K, Post_SD_K, Post_LCL2_K, Post_UCL2_K,
                       Post_Bigger_D5_K,Post_Bigger_D6_K,Post_Bigger_D7_K,Post_Bigger_D8_K,Post_Bigger_D9_K,
                       Post_Smaller_D1_K,Post_Smaller_D2_K,Post_Smaller_D3_K,Post_Smaller_D4_K,Post_Smaller_D5_K,
                       Post_Bigger_D5_Monthly_K,Post_Bigger_D6_Monthly_K,Post_Bigger_D7_Monthly_K,Post_Bigger_D8_Monthly_K,
                       Post_Bigger_D9_Monthly_K,Post_Smaller_D1_Monthly_K,Post_Smaller_D2_Monthly_K,Post_Smaller_D3_Monthly_K,
                       Post_Smaller_D4_Monthly_K,Post_Smaller_D5_Monthly_K,Post_Prob_K)
      names(pixels2_K) = c('pixels_file','Post_mean', 'Post_LCL', 'Post_UCL', 'Post_SD', 'Post_LCL2', 'Post_UCL2',
                           'Post_Bigger_D5', 'Post_Bigger_D6', 'Post_Bigger_D7',
                           'Post_Bigger_D8', 'Post_Bigger_D9', 'Post_Smaller_D1',
                           'Post_Smaller_D2', 'Post_Smaller_D3', 'Post_Smaller_D4',
                           'Post_Smaller_D5','Post_Bigger_D5_Monthly','Post_Bigger_D6_Monthly',
                           'Post_Bigger_D7_Monthly','Post_Bigger_D8_Monthly',
                           'Post_Bigger_D9_Monthly','Post_Smaller_D1_Monthly','Post_Smaller_D2_Monthly',
                           'Post_Smaller_D3_Monthly','Post_Smaller_D4_Monthly','Post_Smaller_D5_Monthly',
                           'Post_Prob_K')
      pixels2_L = list(pixels_file, Post_mean_L, Post_LCL_L, Post_UCL_L, Post_SD_L, Post_LCL2_L, Post_UCL2_L,
                       Post_Bigger_D5_L,Post_Bigger_D6_L,Post_Bigger_D7_L,Post_Bigger_D8_L,Post_Bigger_D9_L,
                       Post_Smaller_D1_L,Post_Smaller_D2_L,Post_Smaller_D3_L,Post_Smaller_D4_L,Post_Smaller_D5_L,
                       Post_Bigger_D5_Monthly_L,Post_Bigger_D6_Monthly_L,Post_Bigger_D7_Monthly_L,Post_Bigger_D8_Monthly_L,
                       Post_Bigger_D9_Monthly_L,Post_Smaller_D1_Monthly_L,Post_Smaller_D2_Monthly_L,Post_Smaller_D3_Monthly_L,
                       Post_Smaller_D4_Monthly_L,Post_Smaller_D5_Monthly_L,
                       Post_Prob_L)
      names(pixels2_L) = c('pixels_file','Post_mean', 'Post_LCL', 'Post_UCL', 'Post_SD', 'Post_LCL2', 'Post_UCL2',
                           'Post_Bigger_D5', 'Post_Bigger_D6', 'Post_Bigger_D7',
                           'Post_Bigger_D8', 'Post_Bigger_D9','Post_Smaller_D1',
                           'Post_Smaller_D2', 'Post_Smaller_D3', 'Post_Smaller_D4',
                           'Post_Smaller_D5','Post_Bigger_D5_Monthly','Post_Bigger_D6_Monthly',
                           'Post_Bigger_D7_Monthly','Post_Bigger_D8_Monthly',
                           'Post_Bigger_D9_Monthly','Post_Smaller_D1_Monthly','Post_Smaller_D2_Monthly',
                           'Post_Smaller_D3_Monthly','Post_Smaller_D4_Monthly','Post_Smaller_D5_Monthly',
                           'Post_Prob_L')
      print('Finally save parameter estimates')
      pars_results = data.frame( Mean = apply(pars, 2, mean, na.rm=T),
                                 SD = apply(pars, 2, sd, na.rm=T),
                                 LCL_95 = apply(pars, 2, quantile, probs = 0.025, na.rm=T),
                                 UCL_95 = apply(pars, 2, quantile, probs = 0.975, na.rm=T),
                                 Median = apply(pars, 2, quantile, probs = 0.5, na.rm=T))
      hyper_pars_results = data.frame( Mean = apply(hyper_pars, 2, mean, na.rm=T),
                                 SD = apply(hyper_pars, 2, sd, na.rm=T),
                                 LCL_95 = apply(hyper_pars, 2, quantile, probs = 0.025, na.rm=T),
                                 UCL_95 = apply(hyper_pars, 2, quantile, probs = 0.975, na.rm=T),
                                 Median = apply(hyper_pars, 2, quantile, probs = 0.5, na.rm=T))
      MONTH_INLA_J_results = data.frame( Mean = apply(MONTH_INLA_J, 2, mean, na.rm=T),
                                       SD = apply(MONTH_INLA_J, 2, sd, na.rm=T),
                                       LCL_95 = apply(MONTH_INLA_J, 2, quantile, probs = 0.025, na.rm=T),
                                       UCL_95 = apply(MONTH_INLA_J, 2, quantile, probs = 0.975, na.rm=T),
                                       Median = apply(MONTH_INLA_J, 2, quantile, probs = 0.5, na.rm=T))
      MONTH_INLA_K_results = data.frame( Mean = apply(MONTH_INLA_K, 2, mean, na.rm=T),
                                       SD = apply(MONTH_INLA_K, 2, sd, na.rm=T),
                                       LCL_95 = apply(MONTH_INLA_K, 2, quantile, probs = 0.025, na.rm=T),
                                       UCL_95 = apply(MONTH_INLA_K, 2, quantile, probs = 0.975, na.rm=T),
                                       Median = apply(MONTH_INLA_K, 2, quantile, probs = 0.5, na.rm=T))
      MONTH_INLA_L_results = data.frame( Mean = apply(MONTH_INLA_L, 2, mean, na.rm=T),
                                       SD = apply(MONTH_INLA_L, 2, sd, na.rm=T),
                                       LCL_95 = apply(MONTH_INLA_L, 2, quantile, probs = 0.025, na.rm=T),
                                       UCL_95 = apply(MONTH_INLA_L, 2, quantile, probs = 0.975, na.rm=T),
                                       Median = apply(MONTH_INLA_L, 2, quantile, probs = 0.5, na.rm=T))
      parameters = list(pars = pars_results, hyper_pars = hyper_pars_results, MONTH_INLA_J = MONTH_INLA_J_results,
                        MONTH_INLA_K = MONTH_INLA_K_results, MONTH_INLA_L = MONTH_INLA_L_results)
      print('success')
    }
    return(list(combined_results = pixels2,
                J_results = pixels2_J,
                K_results = pixels2_K,
                L_results = pixels2_L,
                parameters = parameters))
  }
  results_simple2 = local.plot.field(mesh, samp, E_J = E_J, E_K = E_K, E_L = E_L,
                                     COAST_simp, n_samples = no_samples, 
                                     nx = 300, ny = 300, 
                                     no_T=no_T, funct = exp, mesh_data = covariates_pp,
                                     model = '3',
                                     pixels_file = pixels_file)
  #saveRDS(results_simple2, 'results_spacetime_covars3.rds')
  saveRDS(results_simple2, 'results_space_covars3_MCMC_podSPDE_monthly_dup_removed.rds')
  
  toc = proc.time() - tic
  print('Operation took')
  print(toc[3])
  print(' seconds')
