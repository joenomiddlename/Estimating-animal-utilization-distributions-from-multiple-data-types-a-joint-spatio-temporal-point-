# Function for creating sampling intensity / effort layer

sampling_intensity_script = function(barrier_mesh, barrier_triangles, port_coords, 
                                     port_number_trips, port_countries,
                                     port_trip_length, port_trip_range, barrier_range,
                                     number_time_periods, relative_effort_period,
                                     annual_effort_relative, number_years, number_days_per_period,
                                     cancelled_days_per_period, lost_days_from_follow,
                                     polygons, no_MC_samples, port_sim='yes',
                                     port_lines = list(ind = c(1e9))){
  # barrier mesh is an inla.barrier.mesh (with CRS) with corresponding barrier_triangles
  # Provide a very high res mesh to reduce computational issues with barrier model.
  # port_coords is a SpatialPoints object with CRS
  # w_barrier are the weights (areas) for the dual mesh
  # port_countries is a vector of length == dim(port_coords)[1] giving the country of each port.
  # port_number_trips is a vector of length == dim(port_coords)[1] detailing max number of 
  #                   total trips per day from each port
  # port_trip_length is a function call producing a vector of same size as above detailing 
  #                  the time in hours of each trip per port.
  # port_trip_range is a vector giving the maximum distance in km from each port
  # barrier_range is a scalar giving the range inside land. Note must be > 2xlength of
  # longest mesh triangles.
  # number_time_periods is a scalar giving the number of periods (e.g. months/years) to 
  #                     simulate a unique sampling intensity layer for
  # relative_effort_period is a function call producing a vector of length == number_time_periods.
  #                  representing the relative proportion of total numbers
  #                  of operational WW boats per time period. Usually a prior distn. 
  # annual_effort_relative is a named list of length 2, with function calls. The functions simulate 
  #                  the increase/decrease in the number of boats relative to 2011 for the US and 
  #                  CA respectively. Names must be US and CA. Both produce matrix of 
  #                  dim ==c( nports_region, n_years). Usually a prior distn. 
  # number_years is a scalar giving the number of years of study
  # number_days_per_period is a self explanatory vector of length == number_time_periods
  # cancelled_days_per_period is a function call producing a matrix of dim == c(number_time_periods, nyears) 
  #                       giving the total number of days in each period that the WW companies 
  #                       were cancelled (e.g. due to weather constraints). Usually a prior distn.
  # lost_days_from_follow is a list containing 3 matrices of dim == c(number_time_periods, nyears). 
  #                       These give the total number of days in each period and year that the WW companies 
  #                       lost due to follow. The elements of the list are for J, K and L pods respectively.
  # Pixels is a SpatialPixelsDataframe created from the pixels() function from INLAbru
  # polygons is a SpatialPolygons object with the coastline of interest.
  # no_MC_samples is a scalar giving the number of Monte Carlo fields to produce
  # fitting_mesh is the mesh used for the INLA model fitting
  # port_lines is a list containing 3 arguments 'Points' (a list with elements equal SpatialPoints), 
  # 'Range' (a vector of max ranges to deviate from the lines) and 'ind' (a vector stating which ports
  # need special attentio)
  
  # returns a list with first argument equal to a list of the no_MC_samples copies of 
  # the number of boat hours per period, per port and 
  # second argument equal toa list of the normalized (volume 1) sampling intensity layers per port
  # 3rd argument is a list of plottable fields using gg
  
  # Helper functions 
  # Haakon Bakka's code
  # define a function defining how to plot any spatial fields
  colsc_corr <- function(...) {
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                         limits = c(0.1,1))
  }
  plot.field2 = function(field, poly, mesh, pixels, ...){
    stopifnot(length(field) == mesh$n)
    
    # - choose plotting region to be the same as the study area polygon
    #pixels_plotting = pixels(mesh, mask = poly, nx = 300, ny = 300)
    proj = inla.spde.make.A(mesh, loc = coordinates(pixels) ) 
    #proj = inla.mesh.projector(mesh, xlim = xlim, 
    #                           ylim = ylim, dims=c(300, 300))
    # - Can project from the mesh onto a 300x300 grid 
    #   for plots
    field.proj = as.numeric(proj %*% field)
    
    field_df = SpatialPixelsDataFrame(pixels, 
                                      data = data.frame(field = field.proj))
    # - Do the projection
    #image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
    #           xlim = xlim, ylim = ylim, ...)  
    plot = ggplot() + gg(field_df) + gg(poly) + colsc_corr()
    return(list(plot = plot, pixels_df = field_df))
  }
  
  
  
  # Compute the correlation between the field at any point and this reference location
  local.find.correlation2 = function(Q, location, mesh) {
    sd = sqrt(diag(inla.qinv(Q)))
    # - the marginal standard deviations
    
    A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
    # - create a fake A matrix, to extract the closest mesh node index
    id.node = which.max(A.tmp[1, ])
    # - index of the closest node
    
    print(paste('The location used was c(', 
                round(mesh$loc[id.node, 1], 4), ', ', 
                round(mesh$loc[id.node, 2], 4), ')' ))
    # - location of the closest node
    # - should be close to the location input
    # - sometimes used to plot a black dot
    
    ## Solve a matrix system to find the column of the covariance matrix
    Inode = rep(0, dim(Q)[1]); Inode[id.node] = 1
    covar.column = solve(Q, Inode)
    corr = drop(matrix(covar.column)) / (sd*sd[id.node])
    return(corr)
  }
  
  # Compute the correlation between the field at any point and this set of reference locations
  local.find.correlation3 = function(Q, locations, mesh) {
    sd = sqrt(diag(inla.qinv(Q)))
    # - the marginal standard deviations
    n_locations = dim(locations)[1]
    
    A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(locations, nrow = dim(locations)[1], ncol = 2))
    # - create a fake A matrix, to extract the closest mesh node index
    id.nodes = apply(A.tmp, 1, which.max)# which.max(A.tmp[1, ])
    id.nodes = unique(id.nodes)
    # - index of the closest node
    
    ## Solve a matrix system to find the column of the covariance matrix
    Inode = matrix(0, ncol = length(id.nodes), nrow = dim(Q)[1])
    Inode[cbind(id.nodes,c(1:length(id.nodes)))] = 1
    covar.columns = solve(Q, Inode)
    
    #browser()
    
    for(i in 1:dim(covar.columns)[2])
    {
      if(i == 1)
      {
        corr = drop(matrix(covar.columns[,i])) / (sd*sd[id.nodes[i]])
      }
      if(i != 1)
      {
        corr = pmax(corr, drop(matrix(covar.columns[,i])) / (sd*sd[id.nodes[i]]))
        #corr = corr + drop(matrix(covar.columns[,i])) / (sd*sd[id.nodes[i]])
      }
    }
    
    # for(i in 1:n_locations)
    # {
    #   A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(locations[i,], nrow = 1, ncol = 2))
    #   # - create a fake A matrix, to extract the closest mesh node index
    #   id.node = which.max(A.tmp[1, ])
    #   # - index of the closest node
    #   
    #   ## Solve a matrix system to find the column of the covariance matrix
    #   Inode = rep(0, dim(Q)[1]); Inode[id.node] = 1
    #   covar.column = solve(Q, Inode)
    #   
    #   if(i == 1)
    #   {
    #     corr = drop(matrix(covar.column)) / (sd*sd[id.node])
    #   }
    #   if(i != 1)
    #   {
    #     corr = corr + drop(matrix(covar.column)) / (sd*sd[id.node])
    #   }
    # }
    
    return(corr)
  }
  
  West_decay_f = function(field, x_coords)
  {
    # A function to linearly decay the search effort for all points west of Sooke
    # Based on conversations with experts (Dom and Jason).
    if(length(field) != length(x_coords))
    {
      stop('The x coordinates and field should have the same length')
    }
    Sooke_x = 1167.678
    Min_x = 1100 # West of here we kill the search effort
    
    ratio = ifelse( x_coords < Sooke_x,
                   pmax((x_coords - Min_x)/(Sooke_x - Min_x),0),
                    1)
    
    return(pmin(field, ratio))
  }
  
  # Check the arguments
  # if(is.null(barrier_mesh))
  # {
  #   stop("barrier mesh needs a CRS")
  # }
  # if(class(port_coords) != 'SpatialPoints' | is.null(proj4string(port_coords)) )
  # {
  #   stop('port_coords needs to be a SpatialPoints object with CRS')
  # }
  # if(length(port_number_trips) != dim(port_coords@coords)[1])
  # {
  #   stop('port_number_trips needs length = port_coords')
  # }
  # if(length(port_trip_length) != dim(port_coords@coords)[1])
  # {
  #   stop('port_trip_length needs length = port_coords')
  # }
  # if(length(port_trip_range) != dim(port_coords@coords)[1])
  # {
  #   stop('port_trip_range needs length = port_coords')
  # }
  # if(length(total_fleet_size) != number_time_periods)
  # {
  #   stop('total_fleet_size needs to be length number_time_periods')
  # }
  # 
  # 
  # Step 1 for each port, create a normalized (volume 1) sampling intensity layer.
  number_ports = dim(port_coords@coords)[1]
  
  # define a mapping from the high res mesh to the mesh used in the model
  #mesh_mapping = inla.spde.make.A(mesh=barrier_mesh, loc = fitting_mesh$loc[,1:2])
  
  # which ports need special treatment 
  ind_special = port_lines$ind
  count_special = 1
  
  pixels_plotting =  pixels(barrier_mesh, mask = polygons, nx = 300, ny = 300)
  
  if(port_sim == 'yes')
  {
    for(i in 1:number_ports)
    {
      if(!(i %in% ind_special))
      {
        # define barrier model with land range of barrier_range km.
        barrier.model = inla.barrier.pcmatern(barrier_mesh, barrier_triangles, prior.range = c(port_trip_range[i], 0.5),
                                              prior.sigma = c(1,0.5), range.fraction = barrier_range/port_trip_range[i])
        
        # plot barrier correlation
        Q_barrier = inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(port_trip_range[i]))) # (log sd and log range) 
        
        corr_barrier = local.find.correlation2(Q_barrier, loc = port_coords@coords[i,], barrier_mesh)
        
        corr_plot_barrier = plot.field2(corr_barrier, polygons, barrier_mesh, pixels_plotting)
        corr_plot_barrier$plot = corr_plot_barrier$plot + gg(port_coords[i,]) +
          ggtitle('barrier model correlation field')
        
        # scale by the number of boathours from port each day
        #corr_plot_barrier = corr_plot_barrier * boathours[i] 
        
        # project onto the fitting mesh
        #corr_barrier =  as.numeric(mesh_mapping %*% corr_barrier)
        
        #volume = sum(corr_barrier) #* ( w_barrier/sum(w_barrier) ))
        #field_at_mesh = corr_barrier/volume # make sure the volume integral is 1
        
        #field_at_mesh = (field_at_mesh/volume) #* boathours[i] # scale by the number of boathours from port each day
        
        pixels_field = corr_plot_barrier$pixels_df
        
        if(i == 1)
        {
          plotting_fields = vector('list',number_ports)
          
          plotting_fields[[i]] = corr_plot_barrier$plot #SpatialPixelsDataFrame(pixels, data = pixels_field)
          
          normalized_effort_layers = vector('list',number_ports)
          
          normalized_effort_layers[[i]] = corr_plot_barrier$pixels_df#data.frame(field_at_mesh)
        }
        if(i != 1)
        {
          plotting_fields[[i]] = corr_plot_barrier$plot
          
          normalized_effort_layers[[i]] = corr_plot_barrier$pixels_df#field_at_mesh
        }
      }
      if(i %in% ind_special)
      {
        range_special = port_lines$Range[count_special]
        Points_special = port_lines$Points[[count_special]]@coords
        # define barrier model with land range of barrier_range km.
        barrier.model = inla.barrier.pcmatern(barrier_mesh, barrier_triangles, prior.range = c(range_special, 0.5),
                                              prior.sigma = c(1,0.5), range.fraction = barrier_range/range_special)
        
        # plot barrier correlation
        Q_barrier = inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range_special))) # (log sd and log range) 
        
        print('Creating a search effort layer from polygon... Could be slow')
        corr_barrier =  pmin(local.find.correlation3(Q_barrier, locations = Points_special, barrier_mesh),1)
        corr_barrier = West_decay_f(corr_barrier, barrier_mesh$loc[,1])
        print('Done!')
        # # loop through the points - adding the search effort layers and taking max(sum, 1)
        # for(j in 1:dim(Points_special)[1])
        # {
        #   if(j == 1)
        #   {
        #     browser()
        #     corr_barrier = local.find.correlation2(Q_barrier, loc = Points_special[j,], barrier_mesh)
        #   }
        #   if(j != 1)
        #   {
        #     browser()
        #     corr_barrier = pmin(corr_barrier + local.find.correlation2(Q_barrier, loc = Points_special[j,], barrier_mesh),
        #                         1)
        #   }
        #   print(paste(j,'out of',dim(Points_special)[1],'tracklines complete for port',i, sep = ' '))
        # }
        corr_plot_barrier = plot.field2(corr_barrier, polygons, barrier_mesh, pixels_plotting)
        corr_plot_barrier$plot = corr_plot_barrier$plot + gg(port_coords[i,]) +
          ggtitle('barrier model correlation field')
        print(corr_plot_barrier$plot)
        
        pixels_field = corr_plot_barrier$pixels_df
        
        if(i == 1)
        {
          plotting_fields = vector('list',number_ports)
          
          plotting_fields[[i]] = corr_plot_barrier$plot #SpatialPixelsDataFrame(pixels, data = pixels_field)
          
          normalized_effort_layers = vector('list',number_ports)
          
          normalized_effort_layers[[i]] = corr_plot_barrier$pixels_df#data.frame(field_at_mesh)
        }
        if(i != 1)
        {
          plotting_fields[[i]] = corr_plot_barrier$plot
          
          normalized_effort_layers[[i]] = corr_plot_barrier$pixels_df#field_at_mesh
        }
        count_special = count_special + 1
      }
      
      print(paste('finished making field number',i,'out of',number_ports,'ports'))
    } 
  }
  if(port_sim != 'yes')
  {
    normalized_effort_layers = vector('list',number_ports) # return empty lists
    plotting_fields = vector('list',number_ports)
  }
  
  # create object for storing the boat hours for each pod
  total_boat_hours_per_period_per_port2_J = vector('list',no_MC_samples)
  total_boat_hours_per_period_per_port2_K = vector('list',no_MC_samples)
  total_boat_hours_per_period_per_port2_L = vector('list',no_MC_samples)
  
  # Step 2 Create multiple effort layers by sampling various quantities from their prior distributions
  for(i in 1:no_MC_samples)
  {
    # sample the max number of trips per day from each port
    port_number_trips2 = matrix(rep(port_number_trips, times = number_years), 
                                nrow = number_ports, ncol = number_years,
                                byrow = FALSE) # set equal to 2011 levels
    # NEW INFORMATION - 2014 saw an increase of max 2 trips per day for Steveston
    port_number_trips2[14,6:8] = port_number_trips2[14,6:8] + 2
    
    # sample the increase/decrease across the ports each year from the prior distributions
    port_number_trips2[port_countries == 'US',] = port_number_trips2[port_countries == 'US',] + 
      eval(annual_effort_relative$US)
    port_number_trips2[port_countries == 'CA',] = port_number_trips2[port_countries == 'CA',] + 
      eval(annual_effort_relative$CA)
  
    #browser()
    
    # sample the length/duration of each trip from the port
    port_trip_length2 = eval(port_trip_length)
    
    # sample the amount of time in hours spent stationary
    #time_stationary_per_trip2 = eval(time_stationary_per_trip)
    
    # decrease the search time by the stationary time
    #port_trip_length2 = port_trip_length2 - time_stationary_per_trip2
    
    # sample the relative effort (as a fraction) across the months
    relative_effort_period2 = eval(relative_effort_period)
    
    # sample the number of cancelled days per period, per year based on weather
    cancelled_days_per_period2 = eval(cancelled_days_per_period)
   
    #browser()
    
    # calculate the number of days for each period
    #number_days_per_period
    
    # compute number of boat hours per port as follows:
    # 1 - number of years * number of days
    number_days_per_period2 = matrix(rep(number_days_per_period, times = number_years),
                                     nrow = number_time_periods, ncol = number_years,
                                     byrow = FALSE)
    number_days_per_period2 = number_days_per_period2 - cancelled_days_per_period2
    
    # create the pod-specific number of days - remembering to subtract the lost days of effort
    # due to following (under the monotonic discovery assumption)
    number_days_per_period_J = number_days_per_period2 - lost_days_from_follow[[1]]
    number_days_per_period_K = number_days_per_period2 - lost_days_from_follow[[2]]
    number_days_per_period_L = number_days_per_period2 - lost_days_from_follow[[3]]
    #browser()
    
    # 2 - multiply the number of days per period by the number of boat hours per day, per port, per pod
    total_boat_hours_per_period_per_port_J = matrix(0, nrow = number_ports, ncol = number_time_periods) 
    total_boat_hours_per_period_per_port_K = matrix(0, nrow = number_ports, ncol = number_time_periods)
    total_boat_hours_per_period_per_port_L = matrix(0, nrow = number_ports, ncol = number_time_periods)
    
    # reshape into an nmonth x nyear matrix 
    relative_effort_period2 = matrix(rep(relative_effort_period2, times = number_years), 
                                     nrow = number_time_periods,
                                     ncol = number_years, byrow = FALSE)
    
    #browser()
    
    for(j in 1:number_ports)
    {
      for(k in 1:number_years)
      {
        for(l in 1:number_time_periods)
        {
          total_boat_hours_per_period_per_port_J[j,l] = total_boat_hours_per_period_per_port_J[j,l] +
            number_days_per_period_J[l,k] * (port_number_trips2[j,k] * relative_effort_period2[l] )
          
          total_boat_hours_per_period_per_port_K[j,l] = total_boat_hours_per_period_per_port_K[j,l] +
            number_days_per_period_K[l,k] * (port_number_trips2[j,k] * relative_effort_period2[l] ) 
          
          total_boat_hours_per_period_per_port_L[j,l] = total_boat_hours_per_period_per_port_L[j,l] +
            number_days_per_period_L[l,k] * (port_number_trips2[j,k] * relative_effort_period2[l] ) 
        }
      }
      # scale the normalized effort layer for port j by the total boat hours for months l
    }
    
    #browser()
    
    total_boat_hours_per_period_per_port2_J[[i]] = total_boat_hours_per_period_per_port_J
    total_boat_hours_per_period_per_port2_K[[i]] = total_boat_hours_per_period_per_port_K
    total_boat_hours_per_period_per_port2_L[[i]] = total_boat_hours_per_period_per_port_L
      
    # combine each of the scaled effort layers into one effort layer per day
    # Seperate US and Canadian ports
    
    #total_daily_effort = rowSums(normalized_effort_layers)
    #plotting_fields@data[,(i+1)] = total_daily_effort # see the final layer for plotting
    
    # Step 3 scale that daily field by the relative fleet size and the total number of valid days, per time period
    
    #relative_fleet_size = total_fleet_size/max(total_fleet_size)
    #valid_days_per_period
    
    
    
    
  }
  
  
  return(list(total_boat_hours_per_period_per_port_J = total_boat_hours_per_period_per_port2_J,
              total_boat_hours_per_period_per_port_K = total_boat_hours_per_period_per_port2_K,
              total_boat_hours_per_period_per_port_L = total_boat_hours_per_period_per_port2_L,
              normalized_effort_layers = normalized_effort_layers,
              plotting_fields = plotting_fields))
}