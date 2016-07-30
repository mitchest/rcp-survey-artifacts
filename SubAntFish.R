#########################################################################################
### Code to run modelsused in Foster, Hill & Lyons 2016 paper 		                  ###
### Written by N. Hill Nov 2015.                                                      ###
### Works with RCPmod version 2.142                                                   ###
#########################################################################################


####################
## Data prep
####################

# load libraries and source file and datafile
library(RCPmod)
library(raster)
library(rasterVis)
library(tidyr)


# contains additional functions for data transformations, plotting etc.
source("RCP_Helper_Functions.R")


#read in ID, biological and environmental data file. 
#Note: this file only contains species and environmental variables used in RCP models. 
#Original data and metadata can be found here:
#Hill, N. and Lamb, T. (2015) HIMI Demersal Fish Update and Environmental Covariates. 
#Australian Antarctic Data Centre 
#doi:10.4225/15/5671FDEC717B4

fish<-read.csv("SubAntFish_bioenv.csv")
species <-names(fish)[9:23]

#generate datafile with orthogonal polynomial terms
rcp_env_vars<-c("Long_MP", "log_depth", "caisom_floor_temperature")

rcp_poly<-poly_data(poly_vars=rcp_env_vars, degree=c(2,2,2), 
                     id_vars="HaulIndex",sample_vars="Season", 
                     species_vars=species, data=fish)

rcp_data<-rcp_poly$rcp_data


#Load rasters and create dataframe of prediction space
pred_masked<-brick("pred_masked")


#convert rasters to dataframe and log transform depth
pred_space_rcp<-as.data.frame(rasterToPoints(
  subset(pred_masked, c("Long_MP", "bathymetry", "caisom_floor_temperature"))))
pred_space_rcp<-na.omit(pred_space_rcp)
pred_space_rcp$log_depth<-log(pred_space_rcp$bathymetry* -1)


# Transform using stored polys, predict and plot results
rcp_poly_pred<-poly_pred_space(pred_space_rcp, rcp_poly$poly_output, 
                                sampling_vals="Autumn/Winter", 
                                sampling_name="Season", sampling_factor_levels = c("Autumn/Winter","Spring","summer"))

#create RCP formula
form<- as.formula(paste("cbind(",paste(species, collapse=", "),")~",paste(names(rcp_data)[18:23], collapse="+")))


########################
## Run RCPs
########################

# With Season/Year as sampling effect----

#Note: No seed was set so if running code from here will get slightly different results. 
#This bit  will take a while on a single core machine
# to get the same results and save some time, load the "RCPsamp_fin.RDS" and skip to line 102
nstarts<-1000
max.nRCP<-6

nRCPs_samp <- list()
for( ii in 1:max.nRCP)	
  nRCPs_samp[[ii]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=ii, 
                                       inits="random2", nstart=nstarts, dist="NegBin", mc.cores=1)

#get BICs
RCPsamp_BICs <- sapply( nRCPs_samp, function(x) sapply( x, function(y) y$BIC))
#Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
RCPsamp_minPosteriorSites <- cbind( 181, sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
RCPsamp_ObviouslyBad <- RCPsamp_minPosteriorSites < 2
RCPsamp_BICs[RCPsamp_ObviouslyBad] <- NA

#plot minimum BIC for each nRCP
RCPsamp_minBICs <- apply( RCPsamp_BICs, 2, min, na.rm=TRUE)
plot( 1:max.nRCP, RCPsamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep( 1:max.nRCP, each=nrow( RCPsamp_BICs)), RCPsamp_BICs, pch=20)

#choose 3 RCPs and run best model from above (to get additonal tidbits in model output for later steps)
RCPsamp_goodun <- which.min( RCPsamp_BICs[,3])

control <- list( optimise=FALSE, quiet=FALSE)
RCPsamp_fin<-regimix(form.RCP=form, form.spp=~Season, 
                     nRCP=3, data=rcp_data, dist="NegBin", inits = unlist( nRCPs_samp[[3]][[RCPsamp_goodun]]$coef), control=control)

rm(RCPsamp_BICs,RCPsamp_minPosteriorSites, RCPsamp_ObviouslyBad, RCPsamp_minBICs, RCPsamp_goodun, control)

#plot model diagnostics
# residual plots
plot.regimix(RCPsamp_fin, type="RQR", fitted.scale="log") #looks OK

#Cooks Distance Plots
#takes a while
stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,3,4,5,6,7,8,9,10,20,30,40,50), mc.cores=1, times=500)

# examine dispersion parameter for negative Binomial
hist(RCPsamp_fin$coefs$disp, xlab="Dispersion Parameter", 
     main="Negative Binomial Model", col="grey", cex.main=0.8, cex=0.8, cex.lab=0.8 )

#generate bootstrap estimates of parameters
#again may a while and get slighlty different results. To avoid load "RCPsamp_boots.RDS"
rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=1000, mc.cores=1)

#### Average, SD and CI of species abundances in each  RCP
#for some reason will(?) produce warnings.
RCP_abund_samp<-Sp_abund_all(rcpsamp_boots)

#Get autumn values and make pretty
aut_samp<-as.data.frame(matrix(data=paste0(sprintf("%.2f", round(RCP_abund_samp$autumn$mean,2)), " (", sprintf("%.2f",round(RCP_abund_samp$autumn$sd,2)), ")"),
               ncol=3, nrow=15, byrow=TRUE))
names(aut_samp)<-paste0("RCP", 1:3)
rownames(aut_samp)<- gsub("."," ", species, fixed=TRUE)

## plot of sampling factor effects
sampling_dotplot2(RCPsamp_fin,rcpsamp_boots,legend_fact=c("Spring", "Summer"), col=c("black", "red"), lty=c(1,2))

#Spatial Predictions
RCPsamp_SpPreds<-predict.regimix(object=RCPsamp_fin, object2=rcpsamp_boots, newdata=rcp_poly_pred)
predict_maps2_SDF2(RCPsamp_SpPreds, pred_space=pred_space_rcp, pred_crop=pred_masked, nRCP=3)

###Run models without sampling effect------

#Note: No seed was set so if running code from here will get slightly different results. 
#This bit  will take a while on a single core machine
# to get the same results and save some time, load the "RCPNosamp_fin.RDS" and skip to line 160
nstarts<-1000
max.nRCP<-6

nRCPs_NoSamp <- list()
for( ii in 1:max.nRCP)	
  nRCPs_NoSamp[[ii]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=ii, 
                                         inits="random2", nstart=nstarts, dist="NegBin", mc.cores=1)

#get BICs
RCPNoSamp_BICs <- sapply( nRCPs_NoSamp, function(x) sapply( x, function(y) y$BIC))
#Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
RCPNoSamp_minPosteriorSites <- cbind( 181, sapply( nRCPs_NoSamp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
RCPNoSamp_ObviouslyBad <- RCPNoSamp_minPosteriorSites < 2
RCPNoSamp_BICs[RCPNoSamp_ObviouslyBad] <- NA

#plot minimum BIC for each nRCP
RCPNoSamp_minBICs <- apply( RCPNoSamp_BICs, 2, min, na.rm=TRUE)
plot( 1:max.nRCP, RCPNoSamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep( 1:max.nRCP, each=nrow( RCPNoSamp_BICs)), RCPNoSamp_BICs, pch=20)

#choose 3 RCPs ---
RCPNoSamp_goodun <- which.min( RCPNoSamp_BICs[,3])

control <- list( optimise=FALSE, quiet=FALSE)
RCPNoSamp_fin<-regimix(form.RCP=form, 
                       nRCP=3, data=rcp_data, dist="NegBin", inits = unlist( nRCPs_NoSamp[[3]][[RCPNoSamp_goodun]]$coef), control=control)

rm(RCPNoSamp_BICs,RCPNoSamp_minPosteriorSites, RCPNoSamp_ObviouslyBad, RCPNosamp_minBICs, RCPNoSamp_goodun)

#plot model diagnostics
#residual plot
plot.regimix(RCPNoSamp_fin, type="RQR", fitted.scale="log") 

#Cooks Distance Plots
#will take a while to run
stability.regimix(RCPNoSamp_fin, oosSizeRange=c(1,2,3,4,5,6,7,8,9,10,20,30,40,50), mc.cores=1, times=500)

#generate bootstrap estimates of parameters
#again may a while and get slighlty different results. To avoid load "RCPNoSamp_boots.RDS"
rcpNoSamp_boots<-regiboot(RCPNoSamp_fin, type="BayesBoot", nboot=1000, mc.cores=1)

#Spatial Predictions
RCPNoSamp_SpPreds<-predict.regimix(object=RCPNoSamp_fin, object2=rcpNoSamp_boots, newdata=rcp_poly_pred)
predict_maps2_SDF2(RCPNoSamp_SpPreds, pred_space=pred_space_rcp, pred_crop=pred_masked, nRCP=3)
