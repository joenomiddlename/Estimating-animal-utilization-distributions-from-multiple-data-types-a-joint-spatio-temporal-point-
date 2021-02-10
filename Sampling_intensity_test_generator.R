# Sampling intensity MC sample compiler function
#.libPaths('/zfs/users/joe.watson/joe.watson/rlibs')

library(sp)
library(rgeos)
#library(rgdal)
library(foreach)
library(parallel)
library(doParallel)
#library(doSNOW)
#cl = makeCluster(10)
#registerDoParallel(cl)
registerDoParallel(3)
setwd("~/ownCloud/Whale Project/Final Project Scripts")
SI_WW = readRDS('SI_WW2_dup_removed.rds')
SI_WW2 = SI_WW#readRDS('SI_WW2.rds')
dmesh_barrier_transformed = readRDS('dmesh_barrier_transformed.rds')
dmesh_transformed = readRDS('dmesh_transformed.rds')

grid_length = length(SI_WW$normalized_effort_layers[[1]]@data$field)

#pixels_SI = pixels(mesh_for_SI, COAST_transformed, nx=300,ny=300)

# SI_WW_complete_monthly_J_mean = matrix(0, nrow = length(dmesh_transformed), ncol = 3) 
# SI_WW_complete_monthly_K_mean = matrix(0, nrow = length(dmesh_transformed), ncol = 3)
# SI_WW_complete_monthly_L_mean = matrix(0, nrow = length(dmesh_transformed), ncol = 3)
# 
# SI_WW_complete_monthly_J_var = matrix(0, nrow = length(dmesh_transformed), ncol = 3) 
# SI_WW_complete_monthly_K_var = matrix(0, nrow = length(dmesh_transformed), ncol = 3)
# SI_WW_complete_monthly_L_var = matrix(0, nrow = length(dmesh_transformed), ncol = 3)
# 
# SI_WW_complete_J_mean = rep(0, length(dmesh_transformed)) 
# SI_WW_complete_K_mean = rep(0, length(dmesh_transformed)) 
# SI_WW_complete_L_mean = rep(0, length(dmesh_transformed)) 
# 
# SI_WW_complete_J_var = rep(0, length(dmesh_transformed))  
# SI_WW_complete_K_var = rep(0, length(dmesh_transformed)) 
# SI_WW_complete_L_var = rep(0, length(dmesh_transformed)) 
# 
# SI_WW_complete_monthly_J_mean_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3) 
# SI_WW_complete_monthly_K_mean_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3)
# SI_WW_complete_monthly_L_mean_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3)
# 
# SI_WW_complete_monthly_J_var_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3) 
# SI_WW_complete_monthly_K_var_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3)
# SI_WW_complete_monthly_L_var_barrier = matrix(0, nrow = length(dmesh_barrier_transformed), ncol = 3)
# 
# SI_WW_complete_J_mean_barrier = rep(0, length(dmesh_barrier_transformed)) 
# SI_WW_complete_K_mean_barrier = rep(0, length(dmesh_barrier_transformed)) 
# SI_WW_complete_L_mean_barrier = rep(0, length(dmesh_barrier_transformed)) 
# 
# SI_WW_complete_J_var_barrier = rep(0, length(dmesh_barrier_transformed))  
# SI_WW_complete_K_var_barrier = rep(0, length(dmesh_barrier_transformed)) 
# SI_WW_complete_L_var_barrier = rep(0, length(dmesh_barrier_transformed)) 


# remember that the port fields are stored in SI_WW and the new boat hours in SI_WW2
# remember that the port fields are stored on a regular grid. Need to project onto meshes
sum(SI_WW$normalized_effort_layers[[1]]@data$field) # don't sum to 1 despite the name
effort_layers = SI_WW$normalized_effort_layers

# standardise to sum to 1
for( i in 1:length(effort_layers))
{
  effort_layers[[i]]$field = (effort_layers[[i]]$field / 
                      sum(effort_layers[[i]]$field))
}


SI_temp = SI_WW$normalized_effort_layers[[1]]
SI_temp2 = SI_WW$normalized_effort_layers[[1]]

temp_df1 = data.frame(matrix(0, nrow = length(effort_layers[[1]]$field), ncol = 3))
temp_df2 = data.frame(matrix(0, nrow = length(effort_layers[[1]]$field), ncol = 3*6))

# define our combiner funciton for parallel R
# comb <- function(x, ...) {
#   lapply(seq_along(x),
#          function(i) rbind(x[[i]], lapply(list(...), function(y) y[[i]])))
# }

tic = proc.time()
i=1
# oper <- foreach(i=1:500, .combine='comb', .multicombine=TRUE, .inorder = F, .packages = c('sp','rgeos'),
#                 .init=list(list(), list())) %dopar% {
#                   
                  #.libPaths('/zfs/users/joe.watson/joe.watson/rlibs')
                  
                  #library(sp)
                  #library(rgeos)
                  
                  # compute the total boat hours (sum the months)
                  boat_hours_monthly_port_J = SI_WW2$total_boat_hours_per_period_per_port_J[[i]] 
                  boat_hours_monthly_port_K = SI_WW2$total_boat_hours_per_period_per_port_K[[i]] 
                  boat_hours_monthly_port_L = SI_WW2$total_boat_hours_per_period_per_port_L[[i]] 
                  
                  SI_temp2@data = temp_df1
                  
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
                    
                    #SI_temp@data[,((count+3):(count+5))] = SI_temp@data[,(count:(count+2))]^2 
                    
                    SI_temp2@data[,(1:3)] = SI_temp2@data[,(1:3)] + SI_temp@data[,(count:(count+2))]
                    
                    count = count + 3 #6
                  }
                  
                  #SI_temp2@data[,c(4,5,6)] = SI_temp2@data[,c(1,2,3)]^2
                  
                  # SI_temp2@data$K2 = SI_temp2@data$K^2
                  # SI_temp2@data$L2 = SI_temp2@data$L^2
                  
                  # aggregate in a spatially weighted manner
                  # combine the SI and SI2 together to save time
                  SI_temp@data = cbind(SI_temp@data, SI_temp2@data)
                  
                  # Note that areaWeighted = TRUE always computes weighted average. Need to multiply by area of dmesh
                  temp_agg = aggregate( SI_temp, dmesh_transformed, mean, areaWeighted=TRUE)@data
                  temp_agg_barrier = aggregate( SI_temp, dmesh_barrier_transformed, mean, areaWeighted=TRUE)@data
                  
                  
                  
#                  list(temp_agg, temp_agg_barrier)
#                }
toc = proc.time()
#print('500 iterations took in seconds')
print((toc-tic)[3])

seq_temp = seq(from= 1, by = 3, length.out = 6)

# temp_agg = oper[[1]][[1]]
# temp_agg_barrier = oper[[2]][[1]]
# for(i in 2:length(oper[[1]]))
# {
#   temp_agg = temp_agg + oper[[1]][[i]]
#   temp_agg_barrier = temp_agg_barrier + oper[[2]][[i]]
# }
mesh_n = dim(temp_agg)[1]
barrier_n = dim(temp_agg_barrier)[1]

SI_WW_complete_monthly_J_mean = temp_agg[,seq_temp]
SI_WW_complete_monthly_J_mean_barrier = temp_agg_barrier[,seq_temp]
  
SI_WW_complete_monthly_K_mean = temp_agg[,seq_temp+1]
SI_WW_complete_monthly_K_mean_barrier = temp_agg_barrier[,seq_temp+1]
  
SI_WW_complete_monthly_L_mean = temp_agg[,seq_temp+2]
SI_WW_complete_monthly_L_mean_barrier = temp_agg_barrier[,seq_temp+2]
  
SI_WW_complete_J_mean = temp_agg[,19]
SI_WW_complete_J_mean_barrier = temp_agg_barrier[,19]
  
SI_WW_complete_K_mean = temp_agg[,20]
SI_WW_complete_K_mean_barrier = temp_agg_barrier[,20]
  
SI_WW_complete_L_mean = temp_agg[,21]
SI_WW_complete_L_mean_barrier = temp_agg_barrier[,21]
  
#colMeans(temp_agg_J_month) #temp_agg[,seq_temp] / 1000
#SI_WW_complete_monthly_J_var = apply(temp_agg_J_month, c(1),FUN = function(y){return(diag(var(t(y))))}) #FUN = function(x){apply(x, 2, FUN = function(y){return(diag(var(y)))})}))#diag(var(temp_agg_J_month))#temp_agg[,seq_temp+3] / 1000
rm(temp_agg_J_month)
rm(temp_agg_K_month)
rm(temp_agg_L_month)

# repeat for the barrier mesh
rm(temp_agg_J_month_barrier)
rm(temp_agg_K_month_barrier)
rm(temp_agg_L_month_barrier)

rm(temp_agg_J)
rm(temp_agg_K)
rm(temp_agg_L)

rm(temp_agg_J_barrier)
rm(temp_agg_K_barrier)
rm(temp_agg_L_barrier)

rm(SI_temp,SI_temp2, dmesh_barrier_transformed, dmesh_transformed, effort_layers, SI_WW, SI_WW2,
   temp_df1, temp_df2, seq_temp, i, grid_length, cl, comb)

save.image('WW_SI_compiled_remove_dup.RData')

