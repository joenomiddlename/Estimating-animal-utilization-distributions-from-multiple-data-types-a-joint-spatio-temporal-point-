# 2009_2016_Model_fitting
library(INLA)
library(inlabru)
library(ggplot2)
library(xtable)

setwd("~/ownCloud/Whale Project/Final Project Scripts")
# load pardiso library
# Mac
inla.setOption("pardiso.license", "~/licenses/pardiso.lic")
# Linux

# Load all the precompiled data and fit the competing models
temp = readRDS('Model_fitting_files_scaled_dup_removed.rds')
temp2 = readRDS('Model_fitting_files_scaled_new_constr_dup_removed.rds')
temp3 = readRDS('Model_fitting_files_scaled_depthspde_dup_removed.rds')
temp4 = readRDS('model_fitting_files_podSPDE_dup_removed2.rds')
temp5 = readRDS('Model_fitting_files_scaled_depthspde_podspde_dup_removed.rds')

list2env(temp, .GlobalEnv)
list2env(temp2, .GlobalEnv)
list2env(temp3, .GlobalEnv)
list2env(temp4, .GlobalEnv)
list2env(temp5, .GlobalEnv)
rm(temp, temp2, temp3, temp4, temp5)

mesh_barrier = readRDS('mesh_barrier.rds')
polygon.triangles = readRDS('polygon_triangles.rds')
mesh = readRDS('mesh.rds')

covariates_pp_new <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/covariates_pp_dup_removed.rds")

depth_mesh = inla.mesh.1d(loc = seq(from = min(covariates_pp_new$depth,na.rm=T)-1,
                                    to = max(covariates_pp_new$depth,na.rm=T)+1, length.out = 20),
                          degree = 2)

# Set the priors
barrier.model = inla.barrier.pcmatern(inla.spTransform(mesh_barrier, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")), 
                                      barrier.triangles = polygon.triangles, 
                                      prior.range = c(15, 0.01), 
                                      prior.sigma = c(3, 0.1), 
                                      range.fraction = 0.2) # remember we're not in lon/lat, 3e5 is domain size
# reasonable lower range (median) is half of domain prior.range = c(15e4, .5)
# reasonable upper scale is sd of 2

simple.model = inla.spde2.pcmatern(inla.spTransform(mesh, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")), 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1),
                                   constr = T) 

simple.model.barriermesh = inla.spde2.pcmatern(inla.spTransform(mesh_barrier, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")), 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1),
                                   constr = T) 

Big_Stack = inla.stack(PP_stack_J, PP_stack_K, PP_stack_L)
#Big_Stack$effects$data$Search_Effort = Big_Stack$effects$data$Search_Effort/max(Big_Stack$effects$data$Search_Effort,na.rm = T)
Big_Stack_Barrier = inla.stack(PP_stack_barrier_J, PP_stack_barrier_K, PP_stack_barrier_L)
#Big_Stack_Barrier$effects$data$Search_Effort = Big_Stack_Barrier$effects$data$Search_Effort/max(Big_Stack_Barrier$effects$data$Search_Effort,na.rm = T)
Big_Stack_constr = inla.stack(PP_stack_constr_J, PP_stack_constr_K, PP_stack_constr_L)
#Big_Stack$effects$data$Search_Effort = Big_Stack$effects$data$Search_Effort/max(Big_Stack$effects$data$Search_Effort,na.rm = T)
Big_Stack_constr_Barrier = inla.stack(PP_stack_constr_barrier_J, PP_stack_constr_barrier_K, PP_stack_constr_barrier_L)
#now for spacetime constrained model
Big_Stack_constr_spacetime = inla.stack(PP_stack_constr_spacetime_J, PP_stack_constr_spacetime_K, PP_stack_constr_spacetime_L)
#Big_Stack$effects$data$Search_Effort = Big_Stack$effects$data$Search_Effort/max(Big_Stack$effects$data$Search_Effort,na.rm = T)
Big_Stack_constr_Barrier_spacetime = inla.stack(PP_stack_constr_spacetime_barrier_J, PP_stack_constr_spacetime_barrier_K, PP_stack_constr_spacetime_barrier_L)
#
Big_Stack_podSPDE = inla.stack(PP_stack_J_podSPDE, PP_stack_K_podSPDE, PP_stack_L_podSPDE)
#
Big_Stack_podSPDE_barrier = inla.stack(PP_stack_barrier_J_podSPDE, PP_stack_barrier_K_podSPDE, PP_stack_barrier_L_podSPDE)
rm(PP_stack_J_podSPDE, PP_stack_K_podSPDE, PP_stack_L_podSPDE, PP_stack_barrier_J_podSPDE, PP_stack_barrier_K_podSPDE, PP_stack_barrier_L_podSPDE)

# Restrict the field to be orthogonal to covariates
# load the indices of the mesh vertices lying in water
#mesh_vertices_in_water <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/mesh_vertices_in_water.rds")
#mesh_vertices_in_water_rep = rep(mesh_vertices_in_water, 6) + rep(c(mesh$n*c(0:5)), each = length(mesh_vertices_in_water))
mesh_vertices_in_water = 1:mesh$n
mesh_vertices_in_water_rep = mesh_vertices_in_water#1:(6*mesh$n)

n.covariates = 5#8
n.data = mesh$n#*6#length(mesh_vertices_in_water_rep)

X_water = cbind(rep(1,n.data), covariates_pp_new$SSTmonthavg[mesh_vertices_in_water_rep],
          covariates_pp_new$SSTspatialavg[mesh_vertices_in_water_rep],
          covariates_pp_new$chloromonthavg[mesh_vertices_in_water_rep],
          covariates_pp_new$chlorospatialavg[mesh_vertices_in_water_rep])#,
          #covariates_pp_new$SSTminusspacetime[mesh_vertices_in_water_rep],
          #covariates_pp_new$chlorominusspacetime[mesh_vertices_in_water_rep])#,covariates_pp_new$depth)
# Combine observation locations with water locations
X = rbind(X_water, cbind(rep(1,dim(X_sightings)[1]), 
                         as.matrix(X_sightings[,c('SSTmonthavg','SSTspatialavg','chloromonthavg',
                                  'chlorospatialavg')])))#,'SSTminusspacetime',
                                  #'chlorominusspacetime')])) )

Q = qr.Q(qr(X))
#Q = qr.Q(qr(X_water))

A_constr = inla.spde.make.A(mesh, loc = cbind(rep(mesh$loc[mesh_vertices_in_water,1],times = 1),rep(mesh$loc[mesh_vertices_in_water,2],times = 1)))#, group = rep(1:6, each = mesh$n), n.group = 6)
# join the water locations with the sightings projector matrix
A_constr = rbind(A_constr, A_constr_sightings)

# obtain effort
e_constr = PP_stack_constr_J$data$data$e[1:(mesh$n*6)][mesh_vertices_in_water_rep]
# scale projector matrix by effort
A_constr_effort = A_constr
for(i in 1:length(e_constr))
{
  A_constr_effort[i,] = A_constr_effort[i,] * e_constr[i]
}


# ensure the spde is orthogonal to the column space spanned by the covariates at the water locations
# This is an integral constraint
simple.model.constr = inla.spde2.pcmatern(inla.spTransform(mesh, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")), 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1),
                                   #n.iid.group = 6,
                                   #constr = T,
                                   extraconstr.int = list(A = as.matrix(t(Q)%*%A_constr), e = rep(0, n.covariates))) 

# Next ensure at the observation locations. Regular constraint.
simple.model.constr.obs = inla.spde2.pcmatern(inla.spTransform(mesh, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")), 
                                          prior.range = c(15, 0.01), 
                                          prior.sigma = c(3, 0.1),
                                          #n.iid.group = 6,
                                          #constr = T,
                                          extraconstr = list(A = as.matrix(t(Q)%*%A_constr), e = rep(0, n.covariates))) 

#constr_field_index = inla.spde.make.index('field_constr', n.spde=simple.model.constr$n.spde)

#####

# Models with vessel as a predictor

####

# inla.setOption("pardiso.license", "~/ownCloud/Whale Project/Final Project Scripts/licenses/pardiso.lic")
# simple_fit_covars0_monthly = inla(y ~ -1 + factor(POD_INLA) +
#                                     SSTmonthavg + SSTminusmonth + chloromonthavg + chlorominusmonth +
#                                     vessel + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),
#                                   data = inla.stack.data(Big_Stack),
#                                   family = 'poisson',
#                                   E = e,
#                                   control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                   control.mode = list(theta = c(2.0), restart=TRUE ), 
#                                   control.inla = list(int.strategy = 'eb'),
#                                   verbose = T)
# # DIC 853.31
# # DIC after revision 1099.94
# # DIC new e 
# summary(simple_fit_covars0_monthly)
# saveRDS(simple_fit_covars0_monthly, 'simple_fit_covars0_monthly.rds')
# rm(simple_fit_covars0_monthly)
# simple_fit_covars0_monthly$summary.fixed
# plot(x = 1:6, y = simple_fit_covars0_monthly$summary.random$MONTH_INLA[,2], pch = 'l')
# lines(x = 1:6, y = simple_fit_covars0_monthly$summary.random$MONTH_INLA[,2], pch = 'l')
# 
# simple_fit_covars01_monthly = inla(y ~ -1 + factor(POD_INLA) +
#                                      SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
#                                      SSTminusspacetime + chlorominusspacetime + vessel +  
#                                      f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),
#                                    data = inla.stack.data(Big_Stack),
#                                    family = 'poisson',
#                                    E = e,
#                                    control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                    control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                    control.mode = list(theta = c(1.36), restart=TRUE ), 
#                                    control.inla = list(int.strategy = 'eb'),
#                                    verbose = T)
# # DIC 751.94 - better
# # DIC 1008.13 after revision better
# # 21284.3 new e
# summary(simple_fit_covars01_monthly)
# saveRDS(simple_fit_covars01_monthly, 'simple_fit_covars01_monthly.rds')
# 
# simple_fit_covars01_monthly$summary.fixed
# plot(x = 1:6, y = simple_fit_covars01_monthly$summary.random$MONTH_INLA[,2], pch = 'l')
# lines(x = 1:6, y = simple_fit_covars01_monthly$summary.random$MONTH_INLA[,2], pch = 'l')
# 
# rm(simple_fit_covars01_monthly)
# # Depth variable has no effect and increases DIC of both models - leave out 
# 
# inla.setOption("pardiso.license", "~/licenses/pardiso.lic")
# simple_fit_covars1_monthly = inla(y ~ -1 + f(simple.field, model = simple.model) + factor(POD_INLA) +
#                                     SSTmonthavg + SSTminusmonth + chloromonthavg + chlorominusmonth +
#                                     vessel + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),
#                                   data = inla.stack.data(Big_Stack),
#                                   family = 'poisson',
#                                   E = e,
#                                   control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                   control.mode = list(theta = c(3.3714, 1.8268, 1.7921), restart=TRUE ), 
#                                   control.inla = list(int.strategy = 'eb'),
#                                   verbose = T)
# #DIC -5579.86 with RW 
# #DIC -5335.01 after SI revisions
# #DIC -394.97 new
# summary(simple_fit_covars1_monthly)
# saveRDS(simple_fit_covars1_monthly, 'simple_fit_covars1_monthly.rds')
# rm(simple_fit_covars1_monthly)
# # Model 1 SSTminusmonth negative for model without spatial term, chloromonthavg and chlorominusmonth are positive and vessel is positive
# # Model 2 without spatial term SSTspatialavg and SSTminusspacetime are negative, chlorospatialavg is positive, chlorominusspacetime is negative and vessel is positive
# # Spatial model makes chlorominusmonth negative now!! Spatial Confounding. NO LONGER TRUE
# 
# # obtain plotting files using HPC and load
# results_simple1 = readRDS('results_simple1.rds')
# 
# plot(x = 1:6, y = simple_fit_covars1_monthly$summary.random$MONTH_INLA[,2])
# lines(x = 1:6, y = simple_fit_covars1_monthly$summary.random$MONTH_INLA[,2])
# #plot(x = 1:6, y = simple_fit_covars1_monthly$summary.fixed['SSTmonthavg',1]*(SST_monthlymeans - SST_overallmean))
# pl1_simple <- ggplot() + 
#   gg( SpatialPixelsDataFrame(results_simple1$pixels_file, data = data.frame(results_simple1$Post_mean))[1] ) + 
#   gg(COAST_transformed) + 
#   ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
#   coord_fixed() +
#   colsc(c(-10, 5))
#   #colsc(results_simple1$Post_mean)
# 
# pl1_simple_SD <- ggplot() + 
#   gg(SpatialPixelsDataFrame(results_simple1$pixels_file, data = data.frame(results_simple1$Post_SD))[1] ) + 
#   gg(COAST_transformed) + 
#   ggtitle("SD of LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
#   coord_fixed() +
#   colsc(results_simple1$Post_SD)
# 
# pl2_simple <- ggplot() + 
#   gg(SpatialPixelsDataFrame(results_simple1$pixels_file, data = data.frame(results_simple1$Post_funct_mean))[1]) + 
#   gg(COAST_transformed) + 
#   ggtitle("LGCP fit to Points", subtitle = "(Response Scale)") + 
#   coord_fixed() +
#   colsc(results_simple1$Post_funct_mean)
# 
# multiplot(pl2_simple, pl1_simple, cols = 2)
# multiplot(pl1_simple, pl1_simple_SD, cols = 2)
# 
# par_plot(barrier_fit)
# 
# multiplot(pl1_barrier, pl2, cols = 2)
# 
# # Different parametrisation of the ST covariates
# simple_fit_covars2_monthly = inla(y ~ -1 + f(simple.field, model = simple.model) + factor(POD_INLA) +
#                                     SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
#                                     SSTminusspacetime + chlorominusspacetime + vessel +  
#                                     f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) ,
#                                   data = inla.stack.data(Big_Stack),
#                                   family = 'poisson',
#                                   E = e,
#                                   control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                   control.mode = list(theta = c(2.9952, 1.4274, 1.7127), restart=TRUE ), 
#                                   control.inla = list(int.strategy = 'eb'),
#                                   verbose = T)
# # -6111.59
# # DIC -5294.85 after refined SI
# # DIC -5101.40
# # chloro is negative!
# summary(simple_fit_covars2_monthly)
# saveRDS(simple_fit_covars2_monthly, 'simple_fit_covars2_monthly.rds')
# 
# plot(x = 1:6, y = simple_fit_covars2_monthly$summary.random$MONTH_INLA[,2])
# lines(x = 1:6, y = simple_fit_covars2_monthly$summary.random$MONTH_INLA[,2])
# plot(x = 1:6, y = simple_fit_covars2_monthly$summary.fixed['SSTmonthavg',1]*(SST_monthlymeans - SST_overallmean))
# 
# # compare DICs
# DIC_table=cbind(c(simple_fit_covars0_monthly$dic$dic,simple_fit_covars1_monthly$dic$dic), c(simple_fit_covars01_monthly$dic$dic,simple_fit_covars2_monthly$dic$dic))
# colnames(DIC_table) = c('monthly_centered','spacetime_centered')
# rownames(DIC_table) = c('no spatial','spatial')
# xtable(DIC_table)
# 
# deltaIC(simple_fit_covars0_monthly,simple_fit_covars01_monthly,simple_fit_covars1_monthly,simple_fit_covars2_monthly)
# deltaIC(simple_fit_covars0_monthly,simple_fit_covars01_monthly,simple_fit_covars1_monthly,simple_fit_covars2_monthly, criterion = 'WAIC')
# 
# 
# # print results and save
# xtable(simple_fit_covars0_monthly$summary.fixed)
# xtable(simple_fit_covars01_monthly$summary.fixed)
# xtable(simple_fit_covars1_monthly$summary.fixed)
# xtable(simple_fit_covars2_monthly$summary.fixed)
# xtable(cbind(no_spatial = simple_fit_covars0_monthly$dic$dic, monthly_centered = simple_fit_covars1_monthly$dic$dic, spacetime_centered = simple_fit_covars2_monthly$dic$dic)
# )



# What about the spatio-temporal mode?
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7)) # 0.7 probability correlation > 0.7
#, group = barrier.field.group,
#control.group = list(model="ar1", hyper = list(theta = pcrho))
#inla.setOption("pardiso.license", "~/ownCloud/Whale Project/Final Project Scripts/licenses/pardiso.lic")
#inla.setOption("pardiso.license", "~/licenses/pardiso.lic")

# Complete spatial randomness model
simple_fit_CSR = inla(y ~ -1 + factor(POD_INLA),
                      data = inla.stack.data(Big_Stack),
                      family = 'poisson',
                      E = e,
                      control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                      control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                      #control.mode = list(theta = c(1.7076, 3.9398, 1.1800), restart=TRUE ), 
                      control.inla = list(int.strategy = 'eb'),
                      verbose = T)
summary(simple_fit_CSR)
#
# 3614.16
simple_fit_CSR$summary.fixed


# Spatial-only model
simple_fit_spatial0_monthly = inla(y ~ -1 + factor(POD_INLA) + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
                                     + f(simple.field, model = simple.model),
                                   data = inla.stack.data(Big_Stack),
                                   family = 'poisson',
                                   E = e,
                                   control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                   control.mode = list(theta = c(1.7076, 3.9398, 1.1800), restart=TRUE ), 
                                   control.inla = list(int.strategy = 'eb'),
                                   verbose = T)
# DIC -1632.92
summary(simple_fit_spatial0_monthly)
saveRDS(simple_fit_spatial0_monthly, 'simple_fit_spatial0_monthly.rds')

simple_fit_spatial0_monthly$summary.fixed
plot(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[7:12,2], pch = 'l',col='red')
lines(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[7:12,2], pch = 'l',col='red')
points(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[13:18,2], pch = 'l',col='blue')
lines(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[13:18,2], pch = 'l',col='blue')
points(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[1:6,2], pch = 'l')
lines(x = 1:6, y = simple_fit_spatial0_monthly$summary.random$MONTH_INLA[1:6,2], pch = 'l')

rm(simple_fit_spatial0_monthly)

# spacetime 
simple_fit_spacetime0_monthly = inla(y ~ -1 + factor(POD_INLA) + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
                                       + f(simple.field, model = simple.model, group = simple.field.group,
                                           control.group = list(model="ar1", hyper = list(theta = pcrho))),
                                     data = inla.stack.data(Big_Stack),
                                     family = 'poisson',
                                     E = e,
                                     control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                     control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                     control.mode = list(theta = c(3.1999, 4.0576, 1.2389, 5.1926), restart=TRUE ), 
                                     control.inla = list(int.strategy = 'eb'),
                                     verbose = T)
# -1730.16 - improvement
summary(simple_fit_spacetime0_monthly)
saveRDS(simple_fit_spacetime0_monthly, file = 'simple_fit_spacetime0_monthly.rds')
rm(simple_fit_spacetime0_monthly)

# Covariate-only models
# Test it without a field at all
# First use depth as a linear term
simple_fit_covars03_monthly = inla(y ~ -1 + factor(POD_INLA) +
                                     SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                     SSTminusspacetime + chlorominusspacetime + #depth +  
                                     f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate=POD_INLA),
                                   data = inla.stack.data(Big_Stack),
                                   family = 'poisson',
                                   E = e,
                                   control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                   control.mode = list(theta = c(1.36), restart=TRUE ), 
                                   control.inla = list(int.strategy = 'eb'),
                                   verbose = T)
summary(simple_fit_covars03_monthly)
saveRDS(simple_fit_covars03_monthly,file = 'simple_fit_covars03_monthly_newe.rds')
# DIC with depth 2244.67
# DIC without depth 2707.15
xtable(simple_fit_covars03_monthly$summary.fixed[,1:5])

# Test it without a field at all - now allowing a spline for depth
depth_spde = inla.spde2.pcmatern(depth_mesh, constr = T, prior.range = c(5, 0.95), prior.sigma = c(1,0.1))
Big_Stack_depthspde = inla.stack(PP_stack_J_depth, PP_stack_K_depth, PP_stack_L_depth)
Big_Stack_depthspde_podSPDE = inla.stack(PP_stack_J_depth_podSPDE, PP_stack_K_depth_podSPDE, PP_stack_L_depth_podSPDE)


simple_fit_covars03_monthly2 = inla(y ~ -1 + factor(POD_INLA) +
                                      SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                      SSTminusspacetime + chlorominusspacetime + 
                                      f(depth_smooth, model = depth_spde) +  
                                      f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA), #+
                                      #f(simple.field, model = simple.model),
                                    data = inla.stack.data(Big_Stack_depthspde),
                                    family = 'poisson',
                                    E = e,
                                    control.predictor = list(A=inla.stack.A(Big_Stack_depthspde),compute=FALSE),
                                    control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                    control.mode = list(theta = c(0.1841, 1.4274, 1.3694 ), restart=TRUE ), 
                                    control.inla = list(int.strategy = 'eb'),
                                    verbose = T)
summary(simple_fit_covars03_monthly2)
# DIC is 
saveRDS(simple_fit_covars03_monthly2,file = 'simple_fit_covars03_monthly.rds')

plot(simple_fit_covars03_monthly2$summary.random$depth_smooth$mean, ylab='Effect size', xlab='depth index')
# MASSIVE OVERFITTING - STICK WITH LINEAR

# Finally without depth but without spacetime covariates. Note our spatial models are unable to fit with depth in.
simple_fit_covars05_monthly = inla(y ~ -1 + factor(POD_INLA) +
                                     SSTmonthavg + SSTminusmonth + chloromonthavg + chlorominusmonth + #depth +
                                     #SSTminusspacetime + chlorominusspacetime  +  
                                     f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),
                                   data = inla.stack.data(Big_Stack),
                                   family = 'poisson',
                                   E = e,
                                   control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                   control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                   control.mode = list(theta = c(1.36), restart=TRUE ), 
                                   control.inla = list(int.strategy = 'eb'),
                                   verbose = T)
summary(simple_fit_covars05_monthly)
saveRDS(simple_fit_covars05_monthly, 'simple_fit_covars05_monthly.rds')
# 2349.28 with depth - spacetime covariates helped
# 2842.76 without depth 

## Spatial models with covariates
# Allow a changing monthly effect per pod. Reports are that J pod stays local more than other pods
simple_fit_covars3_monthly_podmonth = inla(y ~ -1 + f(simple.field, model = simple.model) + factor(POD_INLA) +
                                             SSTmonthavg + chloromonthavg + SSTspatialavg + chlorospatialavg +
                                             SSTminusspacetime + chlorominusspacetime + depth + #SSTminusmonth + chlorominusmonth +# + #depth +  
                                             f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                           data = inla.stack.data(Big_Stack),
                                           family = 'poisson',
                                           E = e,
                                           control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                           control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                           control.mode = list(theta = c(3.7170, 0.97, 1.26), restart=TRUE ), 
                                           control.inla = list(int.strategy = 'eb'),
                                           #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                           verbose = T)
summary(simple_fit_covars3_monthly_podmonth)
saveRDS(simple_fit_covars3_monthly_podmonth,file = 'simple_fit_covars3_monthly_newe_podmonth.rds')
# DIC is HUGE with depth - this is due to south west corner of the map with very high depth, but no observations or search effort
# DIC -1842.17 without depth and without spacetime covars. 
# With spacetime covars it is -1851.80



## What about a new spatial field per pod
# Do we need spacetime field?
simple_fit_covars3_monthly_podSPDE = inla(y ~ -1 + f(simple.field, model = simple.model) + # shared by all pods
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
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.7445, 0.5402, 2.0720, 1.1770, 1.1343), restart=TRUE ), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars3_monthly_podSPDE)
saveRDS(simple_fit_covars3_monthly_podSPDE,file = 'simple_fit_covars3_monthly_newe_podSPDE_dup_removed.rds')
# DIC With only a unique L field (following Hauser) (best)
# -1939.73

rw2_df = data.frame(Month = rep(c("May","June","July","August","September","October"),3),
                    x = rep(1:6,3),
                    y = simple_fit_covars3_monthly_podSPDE$summary.random$MONTH_INLA$mean,
                    ymax = simple_fit_covars3_monthly_podSPDE$summary.random$MONTH_INLA$`0.975quant`,
                    ymin = simple_fit_covars3_monthly_podSPDE$summary.random$MONTH_INLA$`0.025quant`,
                    pod = factor(rep(c('J','K','L'),each=6)))
ggplot(rw2_df, aes(x = x, y=y, ymax=ymax, ymin=ymin)) + geom_point() + geom_line() + 
  geom_errorbar() + facet_wrap(~pod) + 
  ggtitle('A plot of the posterior estimated (sum-to-zero constrained) RW2 effects per pod',
          subtitle = '95% posterior credible intervals are shown') +
  xlab('Month') + ylab('Effect Size') + scale_x_continuous(labels = c("May","June","July","August","September","October"),
                                                           breaks = 1:6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

proj_meshpoly <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/proj_meshpoly.rds")
pixels_meshpoly = readRDS('pixels_meshpoly.rds')
pixels_meshpoly = SpatialPixelsDataFrame(pixels_meshpoly, data.frame(L=as.numeric(proj_meshpoly %*% simple_fit_covars3_monthly_podSPDE$summary.random$simple.field.L$mean)))
ggplot() + gg(pixels_meshpoly['L']) + gg(COAST_transformed) + 
  colsc(pixels_meshpoly$L) + ggtitle('Estimated L pod contrast field', subtitle = 'Linear Predictor Scale')

# Repeat but without spacetime covariates
simple_fit_covars2_monthly_podSPDE = inla(y ~ -1 + f(simple.field, model = simple.model) + # shared by all pods
                                            #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                            f(simple.field.L, model = simple.model) + # L contrast
                                            factor(POD_INLA) +
                                            SSTmonthavg + SSTminusmonth + chloromonthavg + chlorominusmonth +
                                            #SSTminusspacetime + chlorominusspacetime + #depth +  
                                            f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                          data = inla.stack.data(Big_Stack_podSPDE),
                                          family = 'poisson',
                                          E = e,
                                          control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE),compute=FALSE),
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.6656, 1.0375, 5.4853, -0.0936, 1.8485), restart=TRUE ), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars2_monthly_podSPDE)
saveRDS(simple_fit_covars2_monthly_podSPDE,file = 'simple_fit_covars2_monthly_newe_podSPDE.rds')
# -1930.69 - not quite as good
#

# try again with depth indicator
# depth crashes model
# from exploratory analysis 100m depth might be a good cutoff.

#Big_Stack_podSPDE$effects$data$depth = Big_Stack_podSPDE$effects$data$topo
#Big_Stack_podSPDE$effects$data$depth[Big_Stack_podSPDE$effects$data$depth>=51] = 51
#Big_Stack_podSPDE$effects$data$depth = log(abs(Big_Stack_podSPDE$effects$data$depth - 52))
#hist(Big_Stack_podSPDE$effects$data$depth[Big_Stack_podSPDE$effects$data$depth >0] )

# can we include depth now?
simple_fit_covars_depthmesh_monthly_podSPDE = inla(y ~ -1 + f(simple.field, model = simple.model) + # shared by all pods
                                            #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                            f(simple.field.L, model = simple.model) + # L contrast
                                            factor(POD_INLA) + f(depth_smooth, model = depth_spde) + 
                                            SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                            SSTminusspacetime + chlorominusspacetime +  
                                            f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                          data = inla.stack.data(Big_Stack_depthspde_podSPDE),
                                          family = 'poisson',
                                          E = e,
                                          control.predictor = list(A=inla.stack.A(Big_Stack_depthspde_podSPDE),compute=FALSE),
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.3510, 1.2235, 5.8004, -0.3818, 1.5931), restart=TRUE ), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars_depthmesh_monthly_podSPDE)
# Appears to be unstable - Massive WAIC - DIC may be 
# Range of spatial fields are now HUGE. Thus little spatial residual remains after including depth spline
# This explains the poorer model fit.
plot(simple_fit_covars_depthmesh_monthly_podSPDE$summary.random$depth_smooth$mean, ylab='Effect size', xlab='depth index')
rug(covariates_pp_new$depth)

## Spacetime models with covariates

# Repeat best model with spacetime field. 
simple_fit_spacetime_covars1_monthly = inla(y ~ -1 + factor(POD_INLA) + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
                                       + f(simple.field, model = simple.model, group = simple.field.group,
                                           control.group = list(model="ar1", hyper = list(theta = pcrho))) +
                                         f(simple.field.L, model = simple.model) +
                                         SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                         SSTminusspacetime + chlorominusspacetime,# +
                                       #vessel,
                                     data = inla.stack.data(Big_Stack_podSPDE),
                                     family = 'poisson',
                                     E = e,
                                     control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE),compute=FALSE),
                                     control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                     control.mode = list(theta = c(3.1999, 4.0576, 1.2389, 5.1926,5,0), restart=TRUE ), 
                                     control.inla = list(int.strategy = 'eb'),
                                     verbose = T)

# DIC not available - model does not converge
summary(simple_fit_spacetime_covars1_monthly)
saveRDS(simple_fit_spacetime_covars1_monthly, file = 'simple_fit_spacetime_covars1_monthly.rds')
rm(simple_fit_spacetime_covars1_monthly)

# simple_fit_spacetime_covars2_monthly = inla(y ~ -1 + factor(POD_INLA) + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
#                                               + f(simple.field, model = simple.model, group = simple.field.group,
#                                                   control.group = list(model="ar1", hyper = list(theta = pcrho))) +
#                                               SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
#                                               SSTminusspacetime + chlorominusspacetime, #+ vessel,
#                                             data = inla.stack.data(Big_Stack),
#                                             family = 'poisson',
#                                             E = e,
#                                             control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                             control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                             control.mode = list(theta = c( 5.7726935, 3.2314243, 1.3880040, 0.9741558), restart=TRUE ), 
#                                             control.inla = list(int.strategy = 'eb'),
#                                             verbose = T)
# 
# # DIC -17776.55 after revised SI and combining BG WW - best model
# # chloro spatial average still negative
# summary(simple_fit_spacetime_covars2_monthly)
# saveRDS(simple_fit_spacetime_covars2_monthly, file = 'simple_fit_spacetime_covars2_monthly.rds')
# rm(simple_fit_spacetime_covars2_monthly)


# simple_fit_spacetime_covars3_monthly = inla(y ~ -1 + factor(POD_INLA) + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
#                                               + f(simple.field, model = simple.model, group = simple.field.group,
#                                                   control.group = list(model="ar1", hyper = list(theta = pcrho))) +
#                                               SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
#                                               SSTminusspacetime + chlorominusspacetime, #+ depth,
#                                             data = inla.stack.data(Big_Stack),
#                                             family = 'poisson',
#                                             E = e,
#                                             control.predictor = list(A=inla.stack.A(Big_Stack)),
#                                             control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
#                                             control.mode = list(theta = c( 9.9618, 4.4545, 2.3765, 3.4666), restart=TRUE ), 
#                                             control.inla = list(int.strategy = 'eb'),
#                                             verbose = T)
# 
# summary(simple_fit_spacetime_covars3_monthly)
# saveRDS(simple_fit_spacetime_covars3_monthly, file = 'simple_fit_spacetime_covars3_monthly.rds')
# DIC is -21321.32 - BEST MODEL!!
# Numerical instabilities! Will not fit. Model too complex.
# DIC is 110279.43 after revising e

# Does fitting a seperate spacetime GRF for each pod improve model fit?
simple_fit_spacetime_covars3_monthly_sep = inla(y ~ -1 + Intercept + f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T) +
                                              + f(simple.field, model = simple.model, group = simple.field.group,
                                                  replicate = simple.field.repl,
                                                  control.group = list(model="ar1", hyper = list(theta = pcrho))) +
                                              SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                              SSTminusspacetime + chlorominusspacetime, #+ vessel,
                                            data = inla.stack.data(Big_Stack),
                                            family = 'poisson',
                                            E = e,
                                            control.predictor = list(A=inla.stack.A(Big_Stack),compute=FALSE),
                                            control.compute = list(config=F, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                            control.mode = list(theta = c(3.5718, 3.7912, 2.4645, 0.9), restart=TRUE,return.marginals=F ), 
                                            control.inla = list(int.strategy = 'eb'),
                                            verbose = T)

# Model is not identifiable. This is understandable given that we have
sum(PP_stack_J$data$data$y,na.rm=T)
sum(PP_stack_K$data$data$y,na.rm=T)
sum(PP_stack_L$data$data$y,na.rm=T)
# limited observations per pod, per month
# I have verified this 

# Finally, prediction using vessel density is questionable biologically.
# Try swapping vessel for depth
# Note we mis-coded vessel - all values >51 are on land and not considered in likelihood
# Big_Stack$effects$data$depth = Big_Stack$effects$data$topo
# Big_Stack$effects$data$depth[Big_Stack$effects$data$depth>=51] = 51
# Big_Stack$effects$data$depth = log(abs(Big_Stack$effects$data$depth - 52))
# hist(Big_Stack$effects$data$depth[Big_Stack$effects$data$depth >0] )

# Spatial only field

# Remove the confounding
# remove the confounding by integral constraint
simple_fit_covars3_monthly_RSR = inla(y ~ -1 + f(simple.field.constr, model = simple.model.constr) + 
                                        factor(POD_INLA) + SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                        SSTminusspacetime + chlorominusspacetime + #depth +  
                                        f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate=POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                      data = inla.stack.data(Big_Stack_constr),
                                      family = 'poisson',
                                      E = e,
                                      control.predictor = list(A=inla.stack.A(Big_Stack_constr),compute=FALSE),
                                      control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                      control.mode = list(theta = c(3.887927, 1.120609, 2.006264), restart=TRUE,return.marginals=F ), 
                                      control.inla = list(int.strategy = 'eb'),
                                      #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                      verbose = T)
summary(simple_fit_covars3_monthly_RSR)
saveRDS(simple_fit_covars3_monthly_RSR,file = 'simple_fit_covars3_monthly_RSR.rds')
# DIC 

# remove the confounding by regular constraint at obs locations
simple_fit_covars3_monthly_RSR_obs = inla(y ~ -1 + f(simple.field.constr, model = simple.model.constr.obs) + 
                                            factor(POD_INLA) + SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                            SSTminusspacetime + chlorominusspacetime + #depth +  
                                            f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                          data = inla.stack.data(Big_Stack_constr),
                                          family = 'poisson',
                                          E = e,
                                          control.predictor = list(A=inla.stack.A(Big_Stack_constr),compute=FALSE),
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso',return.marginals=F),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.887927, 1.120609, 2.006264), restart=TRUE), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars3_monthly_RSR_obs)
saveRDS(simple_fit_covars3_monthly_RSR_obs,file = 'simple_fit_covars3_monthly_RSR_obs.rds')
# DIC 


## Repeat the best models using barrier mesh and barrier model
simple_fit_covars3_monthly_podSPDE_barrier = inla(y ~ -1 + f(barrier.field, model = simple.model.barriermesh) + # shared by all pods
                                            #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                            f(barrier.field.L, model = simple.model.barriermesh) + # L contrast
                                            factor(POD_INLA) +
                                            SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                            SSTminusspacetime + chlorominusspacetime + #depth +  
                                            f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                          data = inla.stack.data(Big_Stack_podSPDE_barrier),
                                          family = 'poisson',
                                          E = e,
                                          control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE_barrier)),
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.9261, 1.1722, 5.2743, -0.2405, 1.2619), restart=TRUE ), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars3_monthly_podSPDE_barrier)
saveRDS(simple_fit_covars3_monthly_podSPDE_barrier,'simple_fit_covars3_monthly_podSPDE_barrier.rds')
# DIC 
# something went wrong. 
# The problem is that I did not smooth the effort layer for the barrier mesh - leading to a very strange likelihood to optimise
# Try Gaussian blur to the WW effort layer for the barrier mesh.
# Fixed
rw2_df2 = data.frame(Month = rep(c("May","June","July","August","September","October"),3),
                    x = rep(1:6,3),
                    y = simple_fit_covars3_monthly_podSPDE_barrier$summary.random$MONTH_INLA$mean,
                    ymax = simple_fit_covars3_monthly_podSPDE_barrier$summary.random$MONTH_INLA$`0.975quant`,
                    ymin = simple_fit_covars3_monthly_podSPDE_barrier$summary.random$MONTH_INLA$`0.025quant`,
                    pod = factor(rep(c('J','K','L'),each=6)))
ggplot(rw2_df2, aes(x = x, y=y, ymax=ymax, ymin=ymin)) + geom_point() + geom_line() + 
  geom_errorbar() + facet_wrap(~pod) + 
  ggtitle('A plot of the posterior estimated (sum-to-zero constrained) RW2 effects per pod',
          subtitle = '95% posterior credible intervals are shown') +
  xlab('Month') + ylab('Effect Size') + scale_x_continuous(labels = c("May","June","July","August","September","October"),
                                                           breaks = 1:6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
proj_barrier_meshpoly <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/proj_barrier_meshpoly.rds")
pixels_meshpoly = readRDS('pixels_meshpoly.rds')
pixels_meshpoly = SpatialPixelsDataFrame(pixels_meshpoly, data.frame(L=as.numeric(proj_barrier_meshpoly %*% simple_fit_covars3_monthly_podSPDE_barrier$summary.random$barrier.field$mean)))
ggplot() + gg(pixels_meshpoly['L']) + gg(COAST_transformed) + 
  colsc(pixels_meshpoly$L) + ggtitle('Estimated L pod contrast field', subtitle = 'Linear Predictor Scale')
# Plots look good. Perhaps the effort e is not scaled correctly. The intercepts are tiny in the model.

# Now use the barrier model
simple_fit_covars3_monthly_podSPDE_barrier_model = inla(y ~ -1 + f(barrier.field, model = barrier.model) + # shared by all pods
                                                    #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                                    f(barrier.field.L, model = barrier.model) + # L contrast
                                                    factor(POD_INLA) +
                                                    SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                                    SSTminusspacetime + chlorominusspacetime + #depth +  
                                                    f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA),#, hyper = list(theta = list(prior="pc.prec", param=c(1,0.001)))) ,
                                                  data = inla.stack.data(Big_Stack_podSPDE_barrier),
                                                  family = 'poisson',
                                                  E = e,
                                                  control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE_barrier)),
                                                  control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                                  control.mode = list(theta = c(2.9351, 2.1030, 5.1639, 0.4024, 1.1176), restart=TRUE ), 
                                                  control.inla = list(int.strategy = 'eb'),
                                                  #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                                  verbose = T)
summary(simple_fit_covars3_monthly_podSPDE_barrier_model)

# Finally - does adding IID noise to capture small-scale variability improve the model fit?
# Constrain with pc priors
simple_fit_covars3_monthly_podSPDE_IID = inla(y ~ -1 + f(simple.field, model = simple.model) + # shared by all pods
                                            #f(simple.field.K, copy = 'simple.field', fixed = F) + # K contrast
                                            f(simple.field.L, model = simple.model) + # L contrast
                                            factor(POD_INLA) +
                                            SSTmonthavg + SSTspatialavg + chloromonthavg + chlorospatialavg +
                                            SSTminusspacetime + chlorominusspacetime + #depth +  
                                            f(MONTH_INLA, model = 'rw2', constr = T, scale.model = T, replicate = POD_INLA) +
                                            f(simple.field.J, model = 'iid', replicate = simple.field.J.group,
                                              hyper = list(theta = list(prior="pc.prec", param=c(3,0.01)))) + # J IID terms with pc prior prob of 0.01 it has sd > 1
                                            f(simple.field.K, replicate = simple.field.K.group, copy = 'simple.field.J', fixed=T) + #model = 'iid', hyper = list(theta = list(prior="pc.prec", param=c(1,0.01)))) + # K IID terms with pc prior prob of 0.01 it has sd > 1
                                            f(simple.field.L2, replicate = simple.field.L2.group, copy = 'simple.field.J', fixed=T), #model = 'iid', hyper = list(theta = list(prior="pc.prec", param=c(1,0.01)))), # L IID terms with pc prior prob of 0.01 it has sd > 1,
                                          data = inla.stack.data(Big_Stack_podSPDE),
                                          family = 'poisson',
                                          E = e,
                                          control.predictor = list(A=inla.stack.A(Big_Stack_podSPDE)),
                                          control.compute = list(config=T, dic = T, cpo = F, waic = T, smtp = 'pardiso'),#openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = F, waic = F, smtp = 'pardiso'), 
                                          control.mode = list(theta = c(3.7445, 0.5402, 2.0720, 1.1770, 1.1343,4), restart=TRUE ), 
                                          control.inla = list(int.strategy = 'eb'),
                                          #control.fixed = list(mean = list(Intercept = 0), prec = list(Intercept = 0.1)),
                                          verbose = T)
summary(simple_fit_covars3_monthly_podSPDE_IID)
saveRDS(simple_fit_covars3_monthly_podSPDE_IID,file = 'simple_fit_covars3_monthly_newe_podSPDE_dup_removed_IID.rds')
# DIC -1975.15 With only a unique L field (following Hauser) AND IID terms with pc prior upper bound of 1 (slightly better than before)
# upper bound changed to 2


#################################################################

#################################################################

## Plot the fitted models 
# plotting and stuff
# predict based on 1000 posterior samples - change by adding n.samples=///
# To increase the resolution of the plot - add nx and ny arguments to pixels function 

# Plotting done on HPC due to HUGE memory requirements
#results_barrier1 = readRDS('results_barrier1.rds')
#results_simple2 <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/results_spacetime_covars2_newborders.rds")
#results_simple3 <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/results_spacetime_covars3.rds")
#results_simple2 <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/results_space_covars3_newe_podSPDE_monthly2.rds")
results_simple2 <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/results_space_covars3_MCMC_podSPDE_monthly_dup_removed.rds")

source('plotting_scripts.R')
#par_plot(barrier_fit_covars1)
#par_plot(simple_fit_covars3_monthly)

colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}

#colsc2 <- function(...) {
#  scale_fill_gradient(low = 'blue', high = 'red')
#}
colsc2 <- function(limits,...) {
scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                     limits = limits)
}

# Plot the linear predictor

results_simple2_linearpredictor = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                                          data = data.frame(results_simple2$Post_mean))
results_simple2_linearpredictorSD = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                                            data = data.frame(results_simple2$Post_SD))
results_simple2_explinearpredictor = SpatialPixelsDataFrame(results_simple2$combined_results$pixels_file,
                                                             data = data.frame(results_simple2$combined_results$Post_funct_median[,1]))
ggplot() + gg(results_simple2_explinearpredictor) + gg(COAST_transformed) + 
  colsc(results_simple2$combined_results$Post_funct_median)
results_simple2_D5 = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                                         data = data.frame(results_simple2$Post_Bigger_D5))
results_simple2_D6 = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                            data = data.frame(results_simple2$Post_Bigger_D6))
results_simple2_D7 = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                            data = data.frame(results_simple2$Post_Bigger_D7))
results_simple2_D8 = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                            data = data.frame(results_simple2$Post_Bigger_D8))
results_simple2_D9 = SpatialPixelsDataFrame(results_simple2$pixels_file,
                                            data = data.frame(results_simple2$Post_Bigger_D9))

#COAST_plotting_transformed = spTransform(COAST_plotting, results_simple2$pixels_file@proj4string)
COAST_transformed <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/COAST_transformed.rds")

# Plotting function allowing user to specify the definition of 'high density' and a threshold
high_density_map_plot = function(results_file, pod = 'All', decile = 9, threshold = 0, poly,
                                 mean_plot = NA, SD_plot = NA, HDR_plot = NA,
                                 LDR_plot = NA, lower_decile = 1, threshold_lower = 0,
                                 LDR_plot_monthly = NA, HDR_plot_monthly = NA,
                                 GMRF_Mean_plot = NA, GMRF_SD_plot = NA,
                                 return_plots = F, mean_plot_lower = NA,
                                 Prob_Pod_plot = NA){
  # results file should be compiled from HPC_mapping_script.R
  # decile should be an integer between 5 and 9 stating the definition of high density
  # threshold should be a probability, specifying the lower cutoff of the plot.
  # poly is a polygon with the same CRS as results_file
  # mean_plot = 1 will plot posterior means
  # sd_plot = 1 will plot posterior sd
  # HDR_plot = 1 (default) will plot posterior probabilities of exceeding decile
  # GMRF_Mean_plot = 1 will plot posterior means of GMRF
  # GMRF_SD_plot = 1 will plot posterior sd of GMRF
  # return_plots is a logical determining if we want to return the plots in a list
  # mean_plot_range is a scalar specifying the lower limit of the density plot.
  plot_list1 = list()
  plot_list2 = list()
  plot_list3 = list()
  plot_list4 = list()
  plot_list5 = list()
  plot_list6 = list()
  plot_list7 = list()
  plot_list8 = list()
  plot_list9 = list()
  
  if(!(pod %in% c('All','J','K','L')))
  {
    stop('pod must be either All, J, K or L')
  }
  if(pod == 'All')
  {
    results_file = results_file$combined_results
    pod_plot_title = 'total SRKW'
  }
  if(pod == 'J')
  {
    results_file = results_file$J_results
    pod_plot_title = 'J pod'
  }
  if(pod == 'K')
  {
    results_file = results_file$K_results
    pod_plot_title = 'K pod'
  }
  if(pod == 'L')
  {
    results_file = results_file$L_results
    pod_plot_title = 'L pod'
  }
  if(!is.na(HDR_plot))
  {
    if(!(decile %in% c(5,6,7,8,9)))
    {
      stop('decile needs to be an integer between 5 and 9')
    }
    upper_val = 100 - (decile * 10) 
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[[paste('Post_Bigger_D',decile,sep='')]]))
    for(i in 1:dim(results_temp@data)[2])
    {
      if(threshold == 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the top ", upper_val, "%",sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc(results_temp@data)
      }
      if(threshold > 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the top ", upper_val, "%",sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc2(limits = c(threshold,1),results_temp@data[,6])
      }
      plot_list1[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list1, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list1)
    }
  }
  if(!is.na(LDR_plot))
  {
    if(!(lower_decile %in% c(1,2,3,4,5)))
    {
      stop('lower_decile needs to be an integer between 1 and 5')
    }
    lower_val = lower_decile * 10 
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[[paste('Post_Smaller_D',lower_decile,sep='')]]))
    for(i in 1:dim(results_temp@data)[2])
    {
      if(threshold_lower == 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the bottom ", lower_val, "%",sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc(results_temp@data)
      }
      if(threshold_lower > 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the bottom ", lower_val, "%",sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc2(limits = c(threshold_lower,1),results_temp@data[,6])
      }
      plot_list6[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list6, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list6)
    }
  }
  if(!is.na(HDR_plot_monthly))
  {
    if(!(decile %in% c(5,6,7,8,9)))
    {
      stop('decile needs to be an integer between 5 and 9')
    }
    upper_val = 100 - (decile * 10) 
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[[paste('Post_Bigger_D',decile,'_Monthly',sep='')]]))
    names(results_temp@data) = c('V1','V2','V3','V4','V5','V6')
    for(i in 1:dim(results_temp@data)[2])
    {
      if(threshold == 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the top ", upper_val, "%", " for that month" ,sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc(results_temp@data)
      }
      if(threshold > 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the top ", upper_val, "%", " for that month",sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc2(limits = c(threshold,1),results_temp@data[,6])
      }
      plot_list7[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list7, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list7)
    }
  }
  if(!is.na(LDR_plot_monthly))
  {
    if(!(lower_decile %in% c(1,2,3,4,5)))
    {
      stop('lower_decile needs to be an integer between 1 and 5')
    }
    lower_val = lower_decile * 10 
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[[paste('Post_Smaller_D',lower_decile,'_Monthly',sep='')]]))
    names(results_temp@data) = c('V1','V2','V3','V4','V5','V6')
    for(i in 1:dim(results_temp@data)[2])
    {
      if(threshold_lower == 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the bottom ", lower_val, "%",' for that month',sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc(results_temp@data)
      }
      if(threshold_lower > 0)
      {
        plot_temp = ggplot() + 
          gg(results_temp[i]) + 
          gg(poly) + labs(fill = "Probability") + 
          ggtitle(paste("Posterior Probability the ", pod_plot_title, " density \nis in the bottom ", lower_val, "%",' for that month',sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
          coord_fixed() +
          colsc2(limits = c(threshold_lower,1),results_temp@data[,6])
      }
      plot_list8[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list8, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list8)
    }
  }
  if(!is.na(mean_plot))
  {
    plot_list = list()
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[['Post_mean']]))
    for(i in 1:dim(results_temp@data)[2])
    {
      plot_temp = ggplot() + 
        gg(results_temp[i]) + 
        gg(poly) + labs(fill = "Mean") + 
        ggtitle(paste("Posterior Mean of the ", pod_plot_title, " density", sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
        coord_fixed() +
        colsc(c(max(c(min(results_temp@data,na.rm=T), mean_plot_lower), na.rm=T), max(results_temp@data,na.rm=T)))
      plot_list2[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list2, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list2)
    }
  }
  if(!is.na(SD_plot))
  {
    plot_list = list()
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[['Post_SD']]))
    for(i in 1:dim(results_temp@data)[2])
    {
      plot_temp = ggplot() + 
        gg(results_temp[i]) + 
        gg(poly) + labs(fill = "SD") + 
        ggtitle(paste("Posterior SD of the ", pod_plot_title, " density", sep=''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
        coord_fixed() +
        colsc(results_temp@data)
      plot_list3[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list3, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list3)
    }
  }
  if(!is.na(GMRF_Mean_plot))
  {
    plot_list = list()
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[['Post_mean_GMRF']]))
    for(i in 1:dim(results_temp@data)[2])
    {
      plot_temp = ggplot() + 
        gg(results_temp[i]) + 
        gg(poly) + labs(fill = "Mean") + 
        ggtitle("Posterior Mean of the GMRF", subtitle = paste("Month of ",months_vec[i],sep='')) + 
        coord_fixed() +
        colsc(results_temp@data)
      plot_list4[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list4, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list4)
    }
  }
  if(!is.na(GMRF_SD_plot))
  {
    plot_list = list()
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(results_file[['Post_SD_GMRF']]))
    for(i in 1:dim(results_temp@data)[2])
    {
      plot_temp = ggplot() + 
        gg(results_temp[i]) + 
        gg(poly) + labs(fill = "SD") + 
        ggtitle("Posterior SD of the GMRF", subtitle = paste("Month of ",months_vec[i],sep='')) + 
        coord_fixed() +
        colsc(results_temp@data)
      plot_list5[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list5, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list5)
    }
  }
  if(!is.na(Prob_Pod_plot))
    {
    plot_list = list()
    months_vec = c("May","June","July","August","September","October")
    results_temp = SpatialPixelsDataFrame(results_file$pixels_file,
                                          data = data.frame(t(results_file[[grep('Post_Prob_',names(results_file), value=T)]])))
    Pod_letter = substr(grep('Post_Prob_',names(results_file), value=T),11,11)
    for(i in 1:dim(results_temp@data)[2])
    {
      plot_temp = ggplot() + 
        gg(results_temp[i]) + 
        gg(poly) + labs(fill = "Probability") + 
        ggtitle(paste("Posterior probability a sighting in each pixel is pod ",Pod_letter,sep = ''), subtitle = paste("Month of ",months_vec[i],sep='')) + 
        coord_fixed() +
        colsc(c(0,1))
      plot_list9[[i]] = plot_temp
    }
    print(multiplot(plotlist = plot_list9, cols = 2, layout = matrix(c(1:6), nrow=3, ncol =2, byrow = T)))
    if(return_plots == T)
    {
      return(plot_list9)
    }
    }
}
# First plot HDR region
high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 8, threshold = 0.95, poly = COAST_transformed, HDR_plot = T)
# First plot LDR region
high_density_map_plot(results_file = results_simple2, pod = 'J', lower_decile = 5, threshold_lower = 0.8, poly = COAST_transformed, LDR_plot = T)
# First plot HDR region, with the exceedance value unique for each month
high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 9, threshold = 0.9, poly = COAST_transformed, HDR_plot_monthly = T)
# First plot LDR region, with the exceedance value unique for each month
high_density_map_plot(results_file = results_simple2, pod = 'J', lower_decile = 5, threshold_lower = 0.9, poly = COAST_transformed, LDR_plot_monthly = T)
# Next plot the Posterior Means
high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, mean_plot = T, HDR_plot = NA)
# Next plot the Posterior SD
high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, SD_plot = T, HDR_plot = NA)
# Next plot the Posterior Mean of the GMRF
high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, GMRF_Mean_plot = T, HDR_plot = NA)
# Next plot the Posterior SD of the GMRF
high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, GMRF_SD_plot = T, HDR_plot = NA)
# Plot probability of the region being a specific pod
high_density_map_plot(results_file = results_simple2, pod = 'L', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA)# Plot the random walk effects

#### Plots for paper ####
# First plot HDR region
plot_1 <- high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 7, threshold = 0, poly = COAST_transformed, HDR_plot = T, return_plots = T)[[1]] +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL)
plot_2 <- high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 7, threshold = 0.95, poly = COAST_transformed, HDR_plot = T, return_plots = T)[[1]] +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL)
# First plot HDR region, with the exceedance value unique for each month
plot_3 <- high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 7, threshold = 0, poly = COAST_transformed, HDR_plot_monthly = T, return_plots = T)[[1]] +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL)
plot_4 <- high_density_map_plot(results_file = results_simple2, pod = 'J', decile = 7, threshold = 0.95, poly = COAST_transformed, HDR_plot_monthly = T, return_plots = T)[[1]] +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL)
# Plot probability of the region being a specific pod
plot_5 <- high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA, return_plots = T)[[1]] +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL) 
plot_6 <- high_density_map_plot(results_file = results_simple2, pod = 'K', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA, return_plots = T)[[1]]+
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL) 
plot_7 <- high_density_map_plot(results_file = results_simple2, pod = 'L', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA, return_plots = T)[[1]]+
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.ticks.length = unit(0, "pt")) + 
  labs(x = NULL, y = NULL) 
# Plot them all
library(cowplot)
ggdraw() +
  draw_plot(plot_1, 0, 2/3, 4/11, 4/11) +
  draw_plot(plot_2, 0.5, 2/3, 4/11, 4/11) +
  draw_plot(plot_3, 0, 1/3, 4/11, 4/11) +
  draw_plot(plot_4, 0.5, 1/3, 4/11, 4/11) +
  draw_plot(plot_5, 0, 0, 4/13, 4/13) +
  draw_plot(plot_6, 1/3, 0, 4/13, 4/13) +
  draw_plot(plot_7, 2/3, 0, 4/13, 4/13) 
#+
#  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

plot_1_1 <- multiplot(plotlist = list(plot_1,plot_2,plot_3,plot_4), cols = 2, layout = matrix(c(1,2,3,4), nrow=2, ncol =2, byrow = T))
multiplot(plotlist = list(plot_1,plot_2,plot_3,plot_4,plot_5,plot_6,plot_7), cols = 3, layout = matrix(c(1,2,NA,3,4,NA,5,6,7), nrow=3, ncol =3, byrow = T))

# loop over the pods and produce the plots and save them to a directory
plot_generator = function()
{
  setwd('~/ownCloud/Whale Project/Plots for final project/Posterior Density plots/Spatial Model no Depth no vessel WITH MC Removed L Duplicates Neworder/')
  
  for(pod_val in c('All','J','K','L'))
  {
    for(upper_decile_val in c(5,6,7,8,9))
    {
      for(threshold in c(0, 0.9))
      {
        # First plot HDR region
        ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, decile = upper_decile_val, threshold = threshold, poly = COAST_transformed, HDR_plot = T), 
               file=paste0("Upper_",10-upper_decile_val , "0_p_greater_", threshold*100, "_Pod", pod_val,".png"), width = 28, height = 20, units = "cm")
        
        # plot HDR region, with the exceedance value unique for each month
        ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, decile = upper_decile_val, threshold = threshold, poly = COAST_transformed, HDR_plot_monthly = T), 
               file=paste0("Upper_",10-upper_decile_val , "0_p_greater_", threshold*100, "_Pod", pod_val,"_monthly.png"), width = 28, height = 20, units = "cm")
      }
    }
    for(lower_decile_val in c(5))
    {
      for(threshold in c(0, 0.9))
      {
        # First plot LDR region
        ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, decile = upper_decile_val, threshold = threshold, poly = COAST_transformed, LDR_plot = T), 
               file=paste0("Lower_",lower_decile_val , "0_p_greater_", threshold*100, "_Pod", pod_val,".png"), width = 28, height = 20, units = "cm")
        
        # plot LDR region, with the exceedance value unique for each month
        ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, decile = upper_decile_val, threshold = threshold, poly = COAST_transformed, LDR_plot_monthly = T), 
               file=paste0("Lower_",lower_decile_val , "0_p_greater_", threshold*100, "_Pod", pod_val,"_monthly.png"), width = 28, height = 20, units = "cm")
      }
    }
    # Next plot the Posterior Means
    ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, poly = COAST_transformed, mean_plot = T, HDR_plot = NA),
           file=paste0("Posterior_Mean_Pod", pod_val,".png"), width = 28, height = 20, units = "cm")
    # Next plot the Posterior SD
    ggsave(high_density_map_plot(results_file = results_simple2, pod = pod_val, poly = COAST_transformed, SD_plot = T, HDR_plot = NA),
           file=paste0("Posterior_SD_Pod", pod_val,".png"), width = 28, height = 20, units = "cm")
  }
  
  # Plot probability of the region being a specific pod
  ggsave(high_density_map_plot(results_file = results_simple2, pod = 'J', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA),
         file=paste0("Posterior_Prob_PodJ.png"), width = 28, height = 20, units = "cm")
  ggsave(high_density_map_plot(results_file = results_simple2, pod = 'K', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA),
         file=paste0("Posterior_Prob_PodK.png"), width = 28, height = 20, units = "cm")
  ggsave(high_density_map_plot(results_file = results_simple2, pod = 'L', poly = COAST_transformed, Prob_Pod_plot = T, HDR_plot = NA),
         file=paste0("Posterior_Prob_PodL.png"), width = 28, height = 20, units = "cm")
  
  # Next plot the Posterior Mean of the GMRF
  ggsave(high_density_map_plot(results_file = results_simple2, pod = 'All', poly = COAST_transformed, GMRF_Mean_plot = T, HDR_plot = NA),
         file=paste0("Posterior_Mean_GMRF.png"), width = 28, height = 20, units = "cm")
  # Next plot the Posterior SD of the GMRF
  ggsave(high_density_map_plot(results_file = results_simple2, pod = 'All', poly = COAST_transformed, GMRF_SD_plot = T, HDR_plot = NA),
         file=paste0("Posterior_SD_GMRF.png"), width = 28, height = 20, units = "cm")
  
}
plot_generator()

rw2_df = data.frame(Month = rep(c("May","June","July","August","September","October"),3),
                    x = rep(1:6,3),
                    y = c(results_simple2$parameters$MONTH_INLA_J$Mean,results_simple2$parameters$MONTH_INLA_K$Mean,results_simple2$parameters$MONTH_INLA_L$Mean),
                    ymax = c(results_simple2$parameters$MONTH_INLA_J$UCL_95,results_simple2$parameters$MONTH_INLA_K$UCL_95,results_simple2$parameters$MONTH_INLA_L$UCL_95),
                    ymin = c(results_simple2$parameters$MONTH_INLA_J$LCL_95,results_simple2$parameters$MONTH_INLA_K$LCL_95,results_simple2$parameters$MONTH_INLA_L$LCL_95),
                    pod = factor(rep(c('J','K','L'),each=6)))
ggplot(rw2_df, aes(x = x, y=y, ymax=ymax, ymin=ymin)) + geom_point() + geom_line() + 
  geom_errorbar() + facet_wrap(~pod) + 
  ggtitle('A plot of the posterior estimated (sum-to-zero constrained) RW2 effects per pod',
          subtitle = '95% posterior credible intervals are shown') +
  xlab('Month') + ylab('Effect Size') + scale_x_continuous(labels = c("May","June","July","August","September","October"),
                                                           breaks = 1:6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at the distribution of the DIC values over the 1000 MC samples
DIC_values <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/DIC_values_dup_removed.rds")
DIC_values = data.frame(DIC = DIC_values - mean(DIC_values))
ggplot(DIC_values, aes(x = DIC)) + geom_histogram() + xlab(expression(paste(Delta, 'DIC'))) + ggtitle('A histogram of the DIC values from the 1000 Monte Carlo model fits', subtitle='Values are relative to the mean of the 1000 values')

# Plot the coefficient estimates with and without MC search effort error
coeffs_df = results_simple2$parameters$pars
coeffs_df$MCMC = 'Yes'
rownames(coeffs_df)[1:3] = c('Pod J', 'Pod K', 'Pod L')
coeffs_df$Par = rownames(coeffs_df)
simple_fit_covars3_monthly_newe_podSPDE2 <- readRDS("~/ownCloud/Whale Project/Final Project Scripts/simple_fit_covars3_monthly_newe_podSPDE_dup_removed.rds")
coeffs_df2 = simple_fit_covars3_monthly_newe_podSPDE2$summary.fixed[,c(1,2,3,5,4)]
coeffs_df2$MCMC = 'No'
rownames(coeffs_df2)[1:3] = c('Pod J', 'Pod K', 'Pod L')
coeffs_df2$Par = rownames(coeffs_df2)
names(coeffs_df2) = names(coeffs_df)
coeffs_df = rbind(coeffs_df,coeffs_df2)
coeffs_df$MCMC = factor(coeffs_df$MCMC)
coeffs_df$Par = factor(coeffs_df$Par, ordered = T, levels=as.character(rownames(coeffs_df2)))

ggplot(coeffs_df, aes(x = Par, y = Mean, ymax = UCL_95, ymin = LCL_95, colour = MCMC)) + 
  geom_point(position=position_dodge(width=0.5)) + 
  geom_errorbar(position=position_dodge(width=0.5),width=0.2) +
  ggtitle('A side-by-side plot of the parameters from the model with and without search effort error',
          subtitle = '95% posterior credible intervals are shown') +
  xlab('Parameter') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0)

# Plot the variance components of the LP, covariates and the GMRF
results_simple2_LPvar = SpatialPixelsDataFrame(results_simple2$combined_results$pixels_file,
                                                data = data.frame(results_simple2$combined_results$pixels_linearpred_spatial_variance))
results_simple2_Covvar = SpatialPixelsDataFrame(results_simple2$combined_results$pixels_file,
                                                 data = data.frame(results_simple2$combined_results$pixels_covariates_spatial_variance))
results_simple2_GMRFvar = SpatialPixelsDataFrame(results_simple2$combined_results$pixels_file,
                                                  data = data.frame(results_simple2$combined_results$pixels_field_spatial_variance))
results_simple2_GMRF_Cov_Covariance = SpatialPixelsDataFrame(results_simple2$combined_results$pixels_file,
                                                 data = data.frame(results_simple2$combined_results$pixels_covariates_field_spatial_covariance))


pl1_simple2_LPvar <- ggplot() + 
  gg(results_simple2_LPvar[1]) + 
  gg(COAST_transformed) + 
  ggtitle("Posterior Variance of Linear Predictor", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(c(results_simple2_Covvar@data[,1],
          results_simple2_LPvar@data[,1],
          results_simple2_GMRFvar@data[,1],
          results_simple2_GMRF_Cov_Covariance@data[,1])) +
  theme(legend.position = "none")

pl1_simple2_Covvar <- ggplot() + 
  gg(results_simple2_Covvar[1]) + 
  gg(COAST_transformed) + 
  ggtitle("Posterior Variance of Covariates", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(c(results_simple2_Covvar@data[,1],
        results_simple2_LPvar@data[,1],
        results_simple2_GMRFvar@data[,1],
        results_simple2_GMRF_Cov_Covariance@data[,1])) +
  theme(legend.position = "none")

pl1_simple2_GMRFvar <- ggplot() + 
  gg(results_simple2_GMRFvar[1]) + 
  gg(COAST_transformed) + 
  ggtitle("Posterior Variance of GMRF", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(c(results_simple2_Covvar@data[,1],
          results_simple2_LPvar@data[,1],
          results_simple2_GMRFvar@data[,1],
          results_simple2_GMRF_Cov_Covariance@data[,1])) +
  theme(legend.position = "none")

pl1_simple2_GMRF_Cov_covariance <- ggplot() + 
  gg(results_simple2_GMRF_Cov_Covariance[1]) + 
  gg(COAST_transformed) + 
  ggtitle("Posterior covariance of GMRF and covariates", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(c(results_simple2_Covvar@data[,1],
          results_simple2_LPvar@data[,1],
          results_simple2_GMRFvar@data[,1]),
           results_simple2_GMRF_Cov_Covariance@data[,1]) +
  theme(legend.position = "none")

multiplot(pl1_simple2_LPvar, pl1_simple2_Covvar, pl1_simple2_GMRFvar, pl1_simple2_GMRF_Cov_covariance, cols = 2)

pl1_simple_D6 <- ggplot() + 
  gg(results_simple2_D6['X1']) + 
  gg(COAST_transformed) + 
  ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_simple2_D6@data[,1])
pl1_simple_D7 <- ggplot() + 
  gg(results_simple2_D7['X1']) + 
  gg(COAST_transformed) + 
  ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_simple2_D7@data[,1])
pl1_simple_D8 <- ggplot() + 
  gg(results_simple2_D8['X1']) + 
  gg(COAST_transformed) + 
  ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_simple2_D8@data[,1])
pl1_simple_D9 <- ggplot() + 
  gg(results_simple2_D9['X1']) + 
  gg(COAST_transformed) + 
  ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_simple2_D9@data[,1])
multiplot(pl1_simple_D5,pl1_simple_D6,pl1_simple_D7,pl1_simple_D8,pl1_simple_D9, cols = 3)

# Posterior Predictive checks for number of sightings
observed_sightings_monthly = data.frame(sightings = as.numeric(t(readRDS('observed_sightings_monthly_dup_removed.rds'))))
observed_sightings_monthly$LCL = as.numeric(results_simple2$combined_results$predicted_sightings_LCL)
observed_sightings_monthly$UCL = as.numeric(results_simple2$combined_results$predicted_sightings_UCL)
observed_sightings_monthly$Median = as.numeric(results_simple2$combined_results$predicted_sightings_Median)
observed_sightings_monthly$Pod = rep(c('J','K','L'), times = 6)
observed_sightings_monthly$Month = rep(1:6, each = 3)
observed_sightings_monthly

ggplot(data = observed_sightings_monthly, aes(x = Month, y = sightings)) +
  geom_point() + geom_errorbar(aes(ymax = UCL, ymin = LCL)) +
  ggtitle('A plot of the observed vs posterior predicted SRKW sightings per Month', 
          subtitle = 'Points are observed values, error bars are 95% posterior credible intervals') +
  geom_hline(yintercept = c(31*8, 30*8))  +
  facet_wrap(~Pod)

###########################################

# Plot the GMRF only


# Make an animating plot using gganimate
library(gganimate)

# reshape the dataframe
results_barrier1_animate = data.frame(Post_mean = as.numeric(results_barrier1$Post_mean),
                                      Post_SD = as.numeric(results_barrier1$Post_SD),
                                      x = rep(coordinates(results_barrier1$pixels_file)[,1],times = dim(results_barrier1$Post_mean)[2]) ,
                                      y = rep(coordinates(results_barrier1$pixels_file)[,2],times = dim(results_barrier1$Post_mean)[2]) ,
                                      month = rep( rep((5:10), each = dim(results_barrier1$Post_mean)[1]), times = 8),
                                      year = rep(2009:2016, each = dim(results_barrier1$Post_mean)[1]*6) )
results_barrier1_animate$YearMonth = paste(results_barrier1_animate$year, results_barrier1_animate$month, sep='-')

anim_plot_mean = ggplot(results_barrier1_animate) + gg(COAST_plotting_transformed) + geom_tile(results_barrier1_animate,
                                                                               mapping = aes_string(x = 'x',y = 'y', fill = 'Post_mean')) +
  coord_fixed() +
  colsc(results_barrier1_animate$Post_mean) +
  labs(title = 'Posterior Mean at Year-Month: {closest_state}') +
  transition_states(YearMonth, transition_length = 2, state_length = 1) +
  ease_aes('linear')
gif_object = animate(anim_plot_mean, renderer = gifski_renderer(),
                     width = 1920, height = 1080)
mpeg_object = animate(anim_plot_mean, renderer = ffmpeg_renderer())
save_animation(mpeg_object,'attempted_video.mp4')
save_animation(gif_object,'attempted_gif.gif')


# what does the plot of the linear predictor averaged over the months and years look like?
results_barrier1_linearpredictor_avg = SpatialPixelsDataFrame(results_barrier1$pixels_file,
                                                              data = data.frame(rowMeans(results_barrier1$Post_mean)))
results_barrier1_linearpredictorSD_avg = SpatialPixelsDataFrame(results_barrier1$pixels_file,
                                                                data = data.frame(rowMeans(results_barrier1$Post_SD)))
pl1_barrier1_avg <- ggplot() + 
  gg(results_barrier1_linearpredictor_avg[1]) + 
  gg(COAST_plotting_transformed) + 
  ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_barrier1_linearpredictor_avg@data[,1])

pl1_barrier1SD_avg <- ggplot() + 
  gg(results_barrier1_linearpredictorSD_avg[1]) + 
  gg(COAST_plotting_transformed) + 
  ggtitle("SD of LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
  coord_fixed() +
  colsc(results_barrier1_linearpredictorSD_avg@data[,1])

multiplot(pl1_barrier1_avg, pl1_barrier1SD_avg, cols = 2)

### Spatial confounding analysis
# investigate the first few eigenvectors of the Q matrix for confounding with the covariates
# NOTE THIS IS SLOW - SHOULD USE SPARSE MATRIX SOLVER IF YOU HAVE TIME 
# eigen_Q = eigen(Q)
# eigen_Q = readRDS('eigen_Q.rds')
low_freq = which(eigen_Q$values %in% sort(eigen_Q$values, decreasing = FALSE)[1:5])
low_freq_evectors = eigen_Q$vectors[,low_freq]
low_freq_evalues = eigen_Q$values[low_freq] # decreasing order

high_freq = which(eigen_Q$values %in% sort(eigen_Q$values, decreasing = TRUE)[1:5])
high_freq_evectors = eigen_Q$vectors[,high_freq]
high_freq_evalues = eigen_Q$values[high_freq] # decreasing order

# compute the correlation between the covariates and the low frequency e'vectors
correlation_lowfreq = cbind(
  cor(cbind(low_freq_evectors, chloro_dmesh_barrier_month$weighted_mean@data))[1:5,6:11],
  cor(cbind(low_freq_evectors, log(abs(topo_dmesh_barrier$weighted_mean$V1-51.3))))[1:5,6],
  cor(cbind(low_freq_evectors, log(vessel_dmesh_barrier$weighted_mean$band1 + 1)))[1:5,6],
  
  cor(cbind(low_freq_evectors, SST_dmesh_barrier_month$weighted_mean@data))[1:5,6:11] )
correlation_lowfreq = rbind( correlation_lowfreq, apply(correlation_lowfreq, 2, FUN = function(x){mean(abs(x))}) )
colnames(correlation_lowfreq) = c('May chloro','Jun chloro','July chloro','Aug chloro','Sep chloro','Oct chloro',
                                  'depth',
                                  'vessel',
                                  'May SST','Jun SST','July SST','Aug SST','Sep SST','Oct SST')
rownames(correlation_lowfreq) = c('5th lowest freq','4th lowest freq','3rd lowest freq','2nd lowest freq','1st lowest freq', 'average abs correlation')
round(correlation_lowfreq,2)
write.csv(round(correlation_lowfreq, 2), 'correlations_lowfreq.csv')
#xtable(correlation_lowfreq)
# strong correlations exist between the lowest frequency e'vectors and ALL the covariates

# compute the correlation between the covariates and the high frequency e'vectors
cor(cbind(high_freq_evectors, chloro_dmesh_barrier_month$weighted_mean@data))[1:5,6:11]
cor(cbind(high_freq_evectors, topo_dmesh_barrier$weighted_mean@data))[1:5,6]
cor(cbind(high_freq_evectors, vessel_dmesh_barrier$weighted_mean@data$band1))[1:5,6]
cor(cbind(high_freq_evectors, SST_dmesh_barrier_month$weighted_mean@data))[1:5,6:11]

# plot the low frequency e'vectors
evect_plot_5 = plot.field(low_freq_evectors[,1],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(low_freq_evectors) + ggtitle('5th smallest evector')
evect_plot_4 = plot.field(low_freq_evectors[,2],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(low_freq_evectors) + ggtitle('4th smallest evector')
evect_plot_3 = plot.field(low_freq_evectors[,3],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(low_freq_evectors) + ggtitle('3rd smallest evector')
evect_plot_2 = plot.field(low_freq_evectors[,4],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(low_freq_evectors) + ggtitle('2nd smallest evector')
evect_plot_1 = plot.field(low_freq_evectors[,5],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(low_freq_evectors) + ggtitle('1st smallest evector')
multiplot(evect_plot_1, evect_plot_2, evect_plot_3, evect_plot_4, evect_plot_5, cols = 2)

# plot the high frequency e'vectors
evect_plot_5_2 = plot.field(high_freq_evectors[,1],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(high_freq_evectors) + ggtitle('5th smallest evector')
evect_plot_4_2 = plot.field(high_freq_evectors[,2],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(high_freq_evectors) + ggtitle('4th smallest evector')
evect_plot_3_2 = plot.field(high_freq_evectors[,3],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(high_freq_evectors) + ggtitle('3rd smallest evector')
evect_plot_2_2 = plot.field(high_freq_evectors[,4],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(high_freq_evectors) + ggtitle('2nd smallest evector')
evect_plot_1_2 = plot.field(high_freq_evectors[,5],COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(high_freq_evectors) + ggtitle('1st smallest evector')
multiplot(evect_plot_1_2, evect_plot_2_2, evect_plot_3_2, evect_plot_4_2, evect_plot_5_2, cols = 2)

# plot the covariates (in their mapped and filtered mesh version)
SST_plot1 = plot.field(SST_dmesh_barrier_month$weighted_mean$May,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('May SST')
SST_plot2 = plot.field(SST_dmesh_barrier_month$weighted_mean$June,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('Jun SST')
SST_plot3 = plot.field(SST_dmesh_barrier_month$weighted_mean$July,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('July SST')
SST_plot4 = plot.field(SST_dmesh_barrier_month$weighted_mean$Aug,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('Aug SST')
SST_plot5 = plot.field(SST_dmesh_barrier_month$weighted_mean$Sep,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('Sep SST')
SST_plot6 = plot.field(SST_dmesh_barrier_month$weighted_mean$Oct,COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)) + ggtitle('Oct SST')
multiplot(SST_plot1,SST_plot2,SST_plot3,SST_plot4,SST_plot5,SST_plot6, cols = 2)

chloro_plot1 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$May),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('May chloro')
chloro_plot2 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$June),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('Jun chloro')
chloro_plot3 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$July),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('July chloro')
chloro_plot4 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$Aug),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('Aug chloro')
chloro_plot5 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$Sep),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('Sep chloro')
chloro_plot6 = plot.field(log(chloro_dmesh_barrier_month$weighted_mean$Oct),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.matrix(log(chloro_dmesh_barrier_month$weighted_mean@data))) + ggtitle('Oct chloro')
multiplot(chloro_plot1,chloro_plot2,chloro_plot3,chloro_plot4,chloro_plot5,chloro_plot6, cols = 2)

vessel_plot = plot.field(log(vessel_dmesh_barrier$weighted_mean$band1 + 1),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.numeric(log(vessel_dmesh_barrier$weighted_mean$band1 + 1))) + ggtitle('Vessel Density')
depth_plot = plot.field(log(abs(topo_dmesh_barrier$weighted_mean$V1-51.3)),COAST_transformed, mesh_barrier_transformed, pixels_plotting) + colsc(as.numeric(log(abs(topo_dmesh_barrier$weighted_mean$V1-51.3)))) + ggtitle('Depth')
multiplot(vessel_plot, depth_plot, cols = 1)



