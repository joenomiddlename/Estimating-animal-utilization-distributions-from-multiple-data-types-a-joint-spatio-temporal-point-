# Encounter simulation results script

# load libraries
library(ggplot2)
library(abind)
library(plyr)
library(mgcv)
library(mgcViz)

# load results
Encounter_Sim_Results <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New.rds")
Encounter_Sim_Results3 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_overlapcorrect.rds")
Encounter_Sim_Results4 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_overlapcorrect2.rds")
Encounter_Sim_Results5 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_overlapcorrect2_fastanimal.rds")
Encounter_Sim_Results6 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_fastanimal.rds")
Encounter_Sim_Results7 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_overlapcorrect_fastanimal.rds")
Encounter_Sim_Results8 <- 
  readRDS("~/ownCloud/Whale Project/Encounters Simulations/Encounter_Sim_Results_New_thinnedpoints2_30.rds")


# Glennie et al 2015 - if observers look ahead, then estimates of intensity biased
# assuming animals move. " bias is considerably smaller in line transect sampling than strip transect sampling"
# Thus in our results, we witness the bias from animal movement.
# We are not assuming anything about the animal's movement process
# Thus, we do not divide by the probability measure over all (measurable)
# paths that the animal could have taken as Glennie et al 2020 JASA
# Remember- we are estimating effort so there is large error with only 1 observer!

dimnames(Encounter_Sim_Results)

# collapse the results into a dataframe
Encounter_Sim_Results <- adply(Encounter_Sim_Results, c(1:7))
names(Encounter_Sim_Results) <- c(names(Encounter_Sim_Results)[1:7],'Value')

Encounter_Sim_Results3 <- adply(Encounter_Sim_Results3, c(1:7))
names(Encounter_Sim_Results3) <- c(names(Encounter_Sim_Results3)[1:7],'Value')

Encounter_Sim_Results4 <- adply(Encounter_Sim_Results4, c(1:7))
names(Encounter_Sim_Results4) <- c(names(Encounter_Sim_Results4)[1:7],'Value')

Encounter_Sim_Results5 <- adply(Encounter_Sim_Results5, c(1:7))
names(Encounter_Sim_Results5) <- c(names(Encounter_Sim_Results5)[1:7],'Value')

Encounter_Sim_Results6 <- adply(Encounter_Sim_Results6, c(1:7))
names(Encounter_Sim_Results6) <- c(names(Encounter_Sim_Results6)[1:7],'Value')

Encounter_Sim_Results7 <- adply(Encounter_Sim_Results7, c(1:7))
names(Encounter_Sim_Results7) <- c(names(Encounter_Sim_Results7)[1:7],'Value')

Encounter_Sim_Results8 <- adply(Encounter_Sim_Results8, c(1:7))
names(Encounter_Sim_Results8) <- c(names(Encounter_Sim_Results8)[1:7],'Value')

head(Encounter_Sim_Results)
head(Encounter_Sim_Results3)

# set the baseline levels to ideal values
Encounter_Sim_Results$Detect_surface <- factor(Encounter_Sim_Results$Detect_surface,levels=c(T, F))
Encounter_Sim_Results$Detect_Range <- factor(Encounter_Sim_Results$Detect_Range,levels=c('10','2','50'))
Encounter_Sim_Results$N_Trips <- factor(Encounter_Sim_Results$N_Trips,levels=c('150','300'))
Encounter_Sim_Results$Sigma_Observer <- factor(Encounter_Sim_Results$Sigma_Observer,
                                               levels=levels(Encounter_Sim_Results$Sigma_Observer)[c(2,1)],
                                               labels = c('High Observer Bias','Low Observer Bias')[c(2,1)])
Encounter_Sim_Results$Quantity <- factor(Encounter_Sim_Results$Quantity,
                                         levels=levels(Encounter_Sim_Results$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results$Observers <- factor(Encounter_Sim_Results$Observers,
                                          levels=levels(Encounter_Sim_Results$Observers)[c(1,2,4,3)],
                                          labels=c('1 Mobile','20 Static','1 Mobile + 20 Static','20 Mobile'))
# create labelling functions
N_Trips_labels <- c(
  '300' = "300 Trips",
  '150' = "150 Trips")
Detect_surface_labels <- c(
  'FALSE' = 'Perfect Detectability assumed',
  'TRUE' = 'Detectability surface modeled'
)
Detect_Range_labels <- c(
  '10' = 'Correct detection range',
  '2' = 'Underestimated',
  '50' = 'Overestimated'
)

### AGAIN for overlap-corrected data
# set the baseline levels to ideal values
Encounter_Sim_Results3$Detect_surface <- factor(Encounter_Sim_Results3$Detect_surface,levels=c(T))
Encounter_Sim_Results3$Detect_Range <- factor(Encounter_Sim_Results3$Detect_Range,levels=c('10'))
Encounter_Sim_Results3$N_Trips <- factor(Encounter_Sim_Results3$N_Trips,levels=c('150','300'))
Encounter_Sim_Results3$Sigma_Observer <- factor(Encounter_Sim_Results3$Sigma_Observer,
                                                levels=levels(Encounter_Sim_Results3$Sigma_Observer),
                                                labels = c('High Observer Bias')[c(1)])
Encounter_Sim_Results3$Quantity <- factor(Encounter_Sim_Results3$Quantity,
                                          levels=levels(Encounter_Sim_Results3$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results3$Observers <- factor(Encounter_Sim_Results3$Observers,
                                           levels=levels(Encounter_Sim_Results3$Observers)[c(1)],
                                           labels=c('20 Mobile'))

Encounter_Sim_Results4$Detect_surface <- factor(Encounter_Sim_Results4$Detect_surface,levels=c(T))
Encounter_Sim_Results4$Detect_Range <- factor(Encounter_Sim_Results4$Detect_Range,levels=c('10'))
Encounter_Sim_Results4$N_Trips <- factor(Encounter_Sim_Results4$N_Trips,levels=c('150','300'))
Encounter_Sim_Results4$Sigma_Observer <- factor(Encounter_Sim_Results4$Sigma_Observer,
                                                levels=levels(Encounter_Sim_Results4$Sigma_Observer),
                                                labels = c('Low Observer Bias')[c(1)])
Encounter_Sim_Results4$Quantity <- factor(Encounter_Sim_Results4$Quantity,
                                          levels=levels(Encounter_Sim_Results4$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results4$Observers <- factor(Encounter_Sim_Results4$Observers,
                                           levels=levels(Encounter_Sim_Results4$Observers)[c(1)],
                                           labels=c('20 Mobile'))

Encounter_Sim_Results5$Detect_surface <- factor(Encounter_Sim_Results5$Detect_surface,levels=c(T))
Encounter_Sim_Results5$Detect_Range <- factor(Encounter_Sim_Results5$Detect_Range,levels=c('10'))
Encounter_Sim_Results5$N_Trips <- factor(Encounter_Sim_Results5$N_Trips,levels=c('150','300'))
Encounter_Sim_Results5$Sigma_Observer <- factor(Encounter_Sim_Results5$Sigma_Observer,
                                                levels=levels(Encounter_Sim_Results5$Sigma_Observer),
                                                labels = c('Low Observer Bias')[c(1)])
Encounter_Sim_Results5$Quantity <- factor(Encounter_Sim_Results5$Quantity,
                                          levels=levels(Encounter_Sim_Results5$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results5$Observers <- factor(Encounter_Sim_Results5$Observers,
                                           levels=levels(Encounter_Sim_Results5$Observers)[c(1)],
                                           labels=c('20 Mobile'))
#levels(Encounter_Sim_Results5$Quantity) <- paste0(levels(Encounter_Sim_Results5$Quantity), '_fastanimal')

# set the baseline levels to ideal values
Encounter_Sim_Results6$Detect_surface <- factor(Encounter_Sim_Results6$Detect_surface,levels=c(T, F))
Encounter_Sim_Results6$Detect_Range <- factor(Encounter_Sim_Results6$Detect_Range,levels=c('10','2','50'))
Encounter_Sim_Results6$N_Trips <- factor(Encounter_Sim_Results6$N_Trips,levels=c('150','300'))
Encounter_Sim_Results6$Sigma_Observer <- factor(Encounter_Sim_Results6$Sigma_Observer,
                                               levels=levels(Encounter_Sim_Results6$Sigma_Observer)[c(2,1)],
                                               labels = c('High Observer Bias','Low Observer Bias')[c(2,1)])
Encounter_Sim_Results6$Quantity <- factor(Encounter_Sim_Results6$Quantity,
                                         levels=levels(Encounter_Sim_Results6$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results6$Observers <- factor(Encounter_Sim_Results6$Observers,
                                          levels=levels(Encounter_Sim_Results6$Observers)[c(1,2,4,3)],
                                          labels=c('1 Mobile','20 Static','1 Mobile + 20 Static','20 Mobile'))

Encounter_Sim_Results7$Detect_surface <- factor(Encounter_Sim_Results7$Detect_surface,levels=c(T))
Encounter_Sim_Results7$Detect_Range <- factor(Encounter_Sim_Results7$Detect_Range,levels=c('10'))
Encounter_Sim_Results7$N_Trips <- factor(Encounter_Sim_Results7$N_Trips,levels=c('150','300'))
Encounter_Sim_Results7$Sigma_Observer <- factor(Encounter_Sim_Results7$Sigma_Observer,
                                                levels=levels(Encounter_Sim_Results7$Sigma_Observer),
                                                labels = c('High Observer Bias')[c(1)])
Encounter_Sim_Results7$Quantity <- factor(Encounter_Sim_Results7$Quantity,
                                          levels=levels(Encounter_Sim_Results7$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results7$Observers <- factor(Encounter_Sim_Results7$Observers,
                                           levels=levels(Encounter_Sim_Results7$Observers)[c(1)],
                                           labels=c('20 Mobile'))

Encounter_Sim_Results8$Detect_surface <- factor(Encounter_Sim_Results8$Detect_surface,levels=c(T))
Encounter_Sim_Results8$Detect_Range <- factor(Encounter_Sim_Results8$Detect_Range,levels=c('10'))
Encounter_Sim_Results8$N_Trips <- factor(Encounter_Sim_Results8$N_Trips,levels=c('150','300'))
Encounter_Sim_Results8$Sigma_Observer <- factor(Encounter_Sim_Results8$Sigma_Observer,
                                                levels=levels(Encounter_Sim_Results8$Sigma_Observer),
                                                labels = c('Low Observer Bias')[c(1)])
Encounter_Sim_Results8$Quantity <- factor(Encounter_Sim_Results8$Quantity,
                                          levels=levels(Encounter_Sim_Results8$Quantity)[c(1,2,3,4,5,7,6,8)])
Encounter_Sim_Results8$Observers <- factor(Encounter_Sim_Results8$Observers,
                                           levels=levels(Encounter_Sim_Results8$Observers)[c(1)],
                                           labels=c('20 Mobile'))
# create labelling functions
N_Trips_labels <- c(
  '300' = "300 Trips",
  '150' = "150 Trips")
Detect_surface_labels <- c(
  'FALSE' = 'Perfect Detectability assumed',
  'TRUE' = 'Detectability surface modeled'
)
Detect_Range_labels <- c(
  '10' = 'Correct detection range',
  '2' = 'Underestimated',
  '50' = 'Overestimated'
)


# What percentage overlap in observer effort?
mean(Encounter_Sim_Results3$Value[801:900]) # 9% high bias
mean(Encounter_Sim_Results4$Value[801:900]) # 5% low bias

#####

# Slow animal results

####


# merge the overlap-corrected data - (run code in section below to shape data)
Encounter_Sim_Results_merged <- rbind(Encounter_Sim_Results3,Encounter_Sim_Results4)
#Encounter_Sim_Results_merged <- Encounter_Sim_Results3
Encounter_Sim_Results_merged$Quantity <- as.character(Encounter_Sim_Results_merged$Quantity)
Encounter_Sim_Results_merged$Quantity[Encounter_Sim_Results_merged$Quantity=='PV'] <- 'PV3'
Encounter_Sim_Results_merged$Quantity[Encounter_Sim_Results_merged$Quantity=='UD_center_bias_y'] <- 'UD_center_bias3_y'
Encounter_Sim_Results_merged3 <- Encounter_Sim_Results8
#Encounter_Sim_Results_merged <- Encounter_Sim_Results3
Encounter_Sim_Results_merged3$Quantity <- as.character(Encounter_Sim_Results_merged3$Quantity)
Encounter_Sim_Results_merged3$Quantity[Encounter_Sim_Results_merged3$Quantity=='PV'] <- 'PV4'
Encounter_Sim_Results_merged3$Quantity[Encounter_Sim_Results_merged3$Quantity=='UD_center_bias_y'] <- 'UD_center_bias4_y'
Encounter_Sim_Results_merged2 <- Encounter_Sim_Results
Encounter_Sim_Results_merged2$Quantity <- as.character(Encounter_Sim_Results_merged2$Quantity)

Encounter_Sim_Results_merged <- rbind(Encounter_Sim_Results_merged,
                                      Encounter_Sim_Results_merged2, 
                                      Encounter_Sim_Results_merged3)
Encounter_Sim_Results_merged$Quantity <- factor(Encounter_Sim_Results_merged$Quantity)
Encounter_Sim_Results_merged$Sigma_Observer <- factor(Encounter_Sim_Results_merged$Sigma_Observer,
                                                      levels = levels(Encounter_Sim_Results_merged$Sigma_Observer)[c(2,1)])
Encounter_Sim_Results_merged$Measure <- 'MSPE' 
Encounter_Sim_Results_merged$Measure[grepl(x=substr(Encounter_Sim_Results_merged$Quantity,1,2),pattern = 'UD')] <- 'Bias' 
Encounter_Sim_Results_merged$Measure <- as.factor(Encounter_Sim_Results_merged$Measure)

##### NOTE - PV is effort-corrected, PV2 is not. 

# function for computing mean, and 95% CI values
# min.mean.sd.max <- function(x) {
#   r <- c(mean(x) - 2*sd(x)/sqrt(length(x)), mean(x), mean(x) + 2*sd(x)/sqrt(length(x)))
#   names(r) <- c("ymin", "y", "ymax")
#   r
# }
# function for computing robust medan, and intervals based on mad, scaled to ensure 
# asymptotic consistency with normally-distributed standard deviation

min.mean.sd.max <- function(x) {
  r <- c(median(x) - 2*mad(x)/sqrt(length(x)), median(x), median(x) + 2*mad(x)/sqrt(length(x)))
  names(r) <- c("ymin", "y", "ymax")
  r
}

#Observers_labels <- c('1 Mobile','20 Static','20 Mobile','1 Mobile + 20 Static')
# PLot the results
ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y')
                             & Encounter_Sim_Results$Detect_Range=='10' 
                             & Encounter_Sim_Results$Detect_surface=='TRUE',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) + 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y')
                             & Encounter_Sim_Results$Sigma_Observer=='High Observer Bias' 
                             & Encounter_Sim_Results$N_Trips=='150',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(Detect_Range),vars(Detect_surface),
             labeller = labeller(Detect_Range=Detect_Range_labels,
                                 Detect_surface=Detect_surface_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# Prediction Variance
ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PV','PV2','PV3')
                             & Encounter_Sim_Results$Detect_Range=='10' 
                             & Encounter_Sim_Results$Detect_surface=='TRUE',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PV','PV2')
                             & Encounter_Sim_Results$Sigma_Observer=='High Observer Bias' 
                             & Encounter_Sim_Results$N_Trips=='150',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(Detect_Range),vars(Detect_surface),
             labeller = labeller(Detect_Range=Detect_Range_labels,
                                 Detect_surface=Detect_surface_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results_merged[Encounter_Sim_Results_merged$Quantity %in% c('PV','PV2','PV3','PV4','UD_center_bias_y','UD_center_bias2_y','UD_center_bias3_y','UD_center_bias4_y')
                             & Encounter_Sim_Results_merged$Detect_Range=='10' 
                             & Encounter_Sim_Results_merged$Detect_surface=='TRUE'
                             & Encounter_Sim_Results_merged$Observers %in% c('20 Mobile'),],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = rep(c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model","Thinned and OC"),2)) +
  scale_linetype_discrete(name = "Quantity", labels = rep(c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model","Thinned and OC"),2)) +
  facet_grid(vars(Measure),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results_merged[Encounter_Sim_Results_merged$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y','UD_center_bias3_y','UD_center_bias4_y')
                             & Encounter_Sim_Results_merged$Detect_Range=='10' 
                             & Encounter_Sim_Results_merged$Detect_surface=='TRUE'
                             & Encounter_Sim_Results_merged$Observers %in% c('20 Mobile')
                             & Encounter_Sim_Results_merged$Sigma_Observer == 'Low Observer Bias',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected", "Uncorrected","Overlap-Corrected","Thinned + O-C")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected", "Uncorrected","Overlap-Corrected","Thinned + O-C")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales='free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# compare the fast animal results
ggplot(Encounter_Sim_Results_merged[Encounter_Sim_Results_merged$Quantity %in% c('PV','PV2','PV3','PV_fastanimal','PV2_fastanimal')
                             & Encounter_Sim_Results_merged$Detect_Range=='10' 
                             & Encounter_Sim_Results_merged$Detect_surface=='TRUE' 
                             & Encounter_Sim_Results_merged$Observers == '20 Mobile',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model","Overlap-Corrected Fast Animal", "Uncorrected Model","Uncorrected Fast Animal","Overlap-Corrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# Prediction Bias
ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PB','PB2')
                             & Encounter_Sim_Results$Detect_Range=='10' 
                             & Encounter_Sim_Results$Detect_surface=='FALSE',],
       aes(x=Observers, y=Value, colour=Quantity)) +
  geom_boxplot() + ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(N_Trips),vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PB','PB2')
                             & Encounter_Sim_Results$Sigma_Observer=='High Observer Bias' 
                             & Encounter_Sim_Results$N_Trips=='150',],
       aes(x=Observers, y=Value, colour=Quantity)) +
  geom_boxplot() + ylab('Average Prediction Bias') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(Detect_Range),vars(Detect_surface),
             labeller = labeller(Detect_Range=Detect_Range_labels,
                                 Detect_surface=Detect_surface_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# Does the magnitude of the observer bias impact the bias-corrected methods?
ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PV')
                             & Encounter_Sim_Results$Detect_Range=='10' 
                             & Encounter_Sim_Results$Detect_surface=='TRUE',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  stat_summary(fun.data = min.mean.sd.max, geom = "point",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model")) +
  facet_grid(vars(N_Trips),vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

#####

# Fast animal results

####


min.mean.sd.max <- function(x) {
  r <- c(median(x) - 2*mad(x)/sqrt(length(x)), median(x), median(x) + 2*mad(x)/sqrt(length(x)))
  names(r) <- c("ymin", "y", "ymax")
  r
}

#Observers_labels <- c('1 Mobile','20 Static','20 Mobile','1 Mobile + 20 Static')
# PLot the results
ggplot(Encounter_Sim_Results6[Encounter_Sim_Results6$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y')
                             & Encounter_Sim_Results6$Detect_Range=='10' 
                             & Encounter_Sim_Results6$Detect_surface=='TRUE',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) + 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results6[Encounter_Sim_Results6$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y')
                             & Encounter_Sim_Results6$Sigma_Observer=='High Observer Bias' 
                             & Encounter_Sim_Results6$N_Trips=='150',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(Detect_Range),vars(Detect_surface),
             labeller = labeller(Detect_Range=Detect_Range_labels,
                                 Detect_surface=Detect_surface_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# Prediction Variance
ggplot(Encounter_Sim_Results6[Encounter_Sim_Results6$Quantity %in% c('PV','PV2','PV3')
                             & Encounter_Sim_Results6$Detect_Range=='10' 
                             & Encounter_Sim_Results6$Detect_surface=='TRUE',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results6[Encounter_Sim_Results6$Quantity %in% c('PV','PV2')
                             & Encounter_Sim_Results6$Sigma_Observer=='High Observer Bias' 
                             & Encounter_Sim_Results6$N_Trips=='150',],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(Detect_Range),vars(Detect_surface),
             labeller = labeller(Detect_Range=Detect_Range_labels,
                                 Detect_surface=Detect_surface_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# merge the overlap-corrected data - (run code in section below to shape data)
Encounter_Sim_Results6_merged <- rbind(Encounter_Sim_Results5, Encounter_Sim_Results7)
#Encounter_Sim_Results6_merged <- Encounter_Sim_Results63
Encounter_Sim_Results6_merged$Quantity <- as.character(Encounter_Sim_Results6_merged$Quantity)
Encounter_Sim_Results6_merged$Quantity[Encounter_Sim_Results6_merged$Quantity=='PV'] <- 'PV3'
Encounter_Sim_Results6_merged$Quantity[Encounter_Sim_Results6_merged$Quantity=='UD_center_bias_y'] <- 'UD_center_bias3_y'
Encounter_Sim_Results6_merged2 <- Encounter_Sim_Results6
Encounter_Sim_Results6_merged2$Quantity <- as.character(Encounter_Sim_Results6_merged2$Quantity)

Encounter_Sim_Results6_merged <- rbind(Encounter_Sim_Results6_merged,
                                       Encounter_Sim_Results6_merged2)
Encounter_Sim_Results6_merged$Quantity <- factor(Encounter_Sim_Results6_merged$Quantity)
Encounter_Sim_Results6_merged$Sigma_Observer <- factor(Encounter_Sim_Results6_merged$Sigma_Observer,
                                                       levels = levels(Encounter_Sim_Results6_merged$Sigma_Observer)[c(2,1)])
ggplot(Encounter_Sim_Results6_merged[Encounter_Sim_Results6_merged$Quantity %in% c('UD_center_bias_y','UD_center_bias2_y','UD_center_bias3_y')
                                     & Encounter_Sim_Results6_merged$Detect_Range=='10' 
                                     & Encounter_Sim_Results6_merged$Detect_surface=='TRUE'
                                     & Encounter_Sim_Results6_merged$Observers %in% c('20 Mobile'),],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab(expression(paste('Bias of estimates of UDs ', mu[Y]))) +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales='free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

ggplot(Encounter_Sim_Results6_merged[Encounter_Sim_Results6_merged$Quantity %in% c('PV','PV2','PV3')
                                     & Encounter_Sim_Results6_merged$Detect_Range=='10' 
                                     & Encounter_Sim_Results6_merged$Detect_surface=='TRUE'
                                     & Encounter_Sim_Results6_merged$Observers %in% c('20 Mobile'),],
       aes(x=Observers, y=Value, colour=Quantity, linetype=Quantity)) +
  #geom_boxplot() + 
  stat_summary(fun.data = min.mean.sd.max, geom = "errorbar",size=1) +
  ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model")) +
  scale_linetype_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model","Overlap-Corrected Model")) +
  facet_grid(,vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels),
             scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')

# Prediction Bias
ggplot(Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% c('PB','PB2')
                             & Encounter_Sim_Results$Detect_Range=='10' 
                             & Encounter_Sim_Results$Detect_surface=='FALSE',],
       aes(x=Observers, y=Value, colour=Quantity)) +
  geom_boxplot() + ylab('Mean Squared Prediction Error') +
  xlab('Observer Type') + 
  geom_hline(yintercept = 0) +
  scale_color_discrete(name = "Quantity", labels = c("Bias-Corrected Model", "Uncorrected Model")) +
  facet_grid(vars(N_Trips),vars(Sigma_Observer),
             labeller = labeller(N_Trips=N_Trips_labels)) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))+ 
  theme(legend.position = 'top')


### fit linear models 
f = Value ~ Observers + Sigma_Observer + N_Trips + 
  Detect_surface + Detect_Range
#as.formula(paste("Value ~", paste(names(Encounter_Sim_Results)[!names(Encounter_Sim_Results) %in% c("Value",'N_rep')], collapse = " * Quantity + ")))

plot_list <- vector(mode = "list", length = length(levels(Encounter_Sim_Results$Quantity))/2)
mod_list_biascorrect <- vector(mode = "list", length = length(levels(Encounter_Sim_Results$Quantity))/2)
mod_list_biasuncorrect <- vector(mode = "list", length = length(levels(Encounter_Sim_Results$Quantity))/2)

fancy_scientific <- function(decimals=0) {
  
  function(x) format(x,digits = decimals,scientific = T)
  
}

for( i in 1:(length(levels(Encounter_Sim_Results$Quantity))/2))
{
  model_biascorrect <- gam(formula=f, 
                           data=Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% levels(Encounter_Sim_Results$Quantity)[((i-1)*2)+1],],
                           family = 'gaussian')
  #anova.gam(model_biascorrect)
  mod_list_biascorrect[[i]] <- model_biascorrect
  
  model_biasuncorrect <- gam(formula=f, 
                             data=Encounter_Sim_Results[Encounter_Sim_Results$Quantity %in% levels(Encounter_Sim_Results$Quantity)[((i-1)*2)+2],],
                             family = 'gaussian')
  #anova.gam(model_biasuncorrect)
  mod_list_biasuncorrect[[i]] <- model_biasuncorrect
  
  plot_list[[i]] <- list(plots=c(plot( getViz(model_biascorrect), allTerms = T )[[1]], plot( getViz(model_biasuncorrect), allTerms = T )[[1]])[c(1,6,2,7,3,8,4,9,5,10)],
                         empty=T)
  
  # add the intercept to each plot and scale axes
  for(j in 1:length(plot_list[[i]]$plots))
  {
    if(round(j/2) == (j/2)) #uncorrected
    {
      plot_list[[i]]$plots[[j]]$ggObj$data[c('y','ty')] <- 
        plot_list[[i]]$plots[[j]]$ggObj$data[c('y','ty')] +
        model_biasuncorrect$coefficients[1]
      
      plot_list[[i]]$plots[[j]]$data$fit[c('y','ty')] <- 
        plot_list[[i]]$plots[[j]]$data$fit[c('y','ty')] +
        model_biasuncorrect$coefficients[1]
      
      plot_list[[i]]$plots[[j]]$ggObj$labels$y <-
        paste('Average response across the groups of',plot_list[[i]]$plots[[j]]$ggObj$labels$x)
      
      limits <- range(range(plot_list[[i]]$plots[[j]]$data$fit[c('y','ty')],
                            plot_list[[i]]$plots[[j-1]]$data$fit[c('y','ty')]) + 
                        c(-2.1,2.1)*max(plot_list[[i]]$plots[[j]]$data$fit$se, 
                                        plot_list[[i]]$plots[[j-1]]$data$fit$se))
      
      
      plot_list[[i]]$plots[[j]]$ggObj <- plot_list[[i]]$plots[[j]]$ggObj +
        scale_y_continuous(limits = limits,
                           labels = fancy_scientific(2))
      plot_list[[i]]$plots[[j-1]]$ggObj <- plot_list[[i]]$plots[[j-1]]$ggObj +
        scale_y_continuous(limits = limits,
                           labels = fancy_scientific(2))
    }
    if(round(j/2) != (j/2)) # corrected
    {
      plot_list[[i]]$plots[[j]]$ggObj$data[c('y','ty')] <- 
        plot_list[[i]]$plots[[j]]$ggObj$data[c('y','ty')] +
        model_biascorrect$coefficients[1]
      
      plot_list[[i]]$plots[[j]]$data$fit[c('y','ty')] <- 
        plot_list[[i]]$plots[[j]]$data$fit[c('y','ty')] +
        model_biascorrect$coefficients[1]
      
      plot_list[[i]]$plots[[j]]$ggObj$labels$y <-
        paste('Average response across the groups of',plot_list[[i]]$plots[[j]]$ggObj$labels$x)
      
    }
  }
  
  class(plot_list[[i]]) <- c("plotGam","gg")
}
#names(plot_list) <- c('Prediction Variance','Prediction Bias','UD Center Bias X','UD Center Bias Y')

print(plot_list[[1]], pages=5, ask = F, top='APV for corrected vs uncorrected models')
print(plot_list[[2]], pages=5, ask = F, top='APB for corrected vs uncorrected models')
print(plot_list[[3]], pages=5, ask = F, top='Bias of UD_x for corrected vs uncorrected models')
print(plot_list[[4]], pages=5, ask = F, top='Bias of UD_y for corrected vs uncorrected models')

summary(mod_list_biascorrect[[1]])
summary(mod_list_biasuncorrect[[1]])

summary(mod_list_biascorrect[[2]])
summary(mod_list_biasuncorrect[[2]])

summary(mod_list_biascorrect[[3]])
summary(mod_list_biasuncorrect[[3]])

summary(mod_list_biascorrect[[4]])
summary(mod_list_biasuncorrect[[4]])


# create labelling functions
#N_Trips_labels <- c(
#  '300' = "300 Trips",
#  '150' = "150 Trips")
#Detect_surface_labels <- c(
#  'TRUE' = 'Detectability surface modeled'
#)
#Detect_Range_labels <- c(
#  '10' = 'Correct detection range'
#)
