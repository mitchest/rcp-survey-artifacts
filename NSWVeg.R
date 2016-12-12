## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## from Foster, Hill & Lyons, 2017, JRSS (Series C, Applied Stats)
## No licence, though please attribute us and publish any derivitives open source =)

library(RCPmod) # v2.128
library(dplyr)


# load data ---------------------------------------------------------------

# # An .RData file is a good option if using more often - uncompressed .csv is a little unwieldly in this case
load("NSWVeg_covariates_species.RData")
# # .csv if that's what you prefer
# covariates.species <- read.csv("NSWVeg_covariates_species.csv", header=T, stringsAsFactors=F)



# parameratise models -----------------------------------------------------
# If you are seriously wanting to reproduce this analysis, you will need access to a computing cluster
# You should plan on using ~60-70 days of CPU time, or ~8 hours on ~200 cores
# If that's something you might want to do, then head over to this repo:
# https://github.com/mitchest/rcp-cluster-parallelisation

# Alternatively you can use the subsetting controls below to run with a subset of data
# In that case you will want to set the arr.job=FALSE, and subset.data=TRUE with subset.size=n

arr.job <- TRUE # is this an array job or a single job
if (arr.job) {
  # get job number from array job env. var.
  job <- Sys.getenv("PBS_ARRAYID")
  job <- as.numeric(job)
  # should the models be fit with species formula?
  # should be one of "species", "nospecies", "both"
  # should be passed from the PBS batch job
  if ("both" %in% commandArgs()) {species.model <- "both"}
  if ("species" %in% commandArgs()) {species.model <- "species"}
  if ("nospecies" %in% commandArgs()) {species.model <- "nospecies"}
  # define nRCP (number of communities/RCPs)
  starts <- rep(2:30,100)
  nRCP <- starts[job]
} else {
  # choose your nRCP/model type if not array job
  nRCP <- 3
  species.model <- "nospecies"
}


# use all communities?
subset.data <- FALSE # specify T/F to subset data
subset.size <- 1000 # specifcy random subset size


# conditional processing to subset data to subset - ensure sitename column is correctly specified
if (subset.data) {
  # specify the  subset to use in the analysis
  set.seed(subset.size)
  sample.sites <- sample(covariates.species$SiteNo, subset.size)
  covariates.species <- covariates.species[covariates.species$SiteNo %in% sample.sites,]
  print(paste0("Successfully subsetted [",subset.size,"] random sites"))
} else {
  print("No subsetting performed")
}

# record site order
site.names <- covariates.species$SiteNo

# define where abundance data starts in covariates.species, test with names(covariates.species)[n.abund]
n.abund <- 18

# choose species to model with
# create a list of species with occurance > n to use for modelling
species.n <- 50 # change this to desired minimum occurance
species.count <- data.frame(count=sort(colSums(covariates.species[,n.abund:ncol(covariates.species)]), decreasing=T))
species.count$species <- row.names(species.count)
model.species.vector <- species.count$species[species.count$count>species.n]
model.species.string <- paste0("cbind(", paste(model.species.vector, collapse=","),")")

# choose which covariates to use
model.covariates.string <- "tempmtcp.1+tempmtcp.2+precipseas.1+precipseas.2+precipann.1+precipann.2+tempmtwp.1+tempmtwp.2+rough500.1+rough500.2+bd200.1+bd200.2+ph200.1+ph200.2"
model.covariates.vector <- c("tempmtcp","precipseas","precipann","tempmtwp","rough500","bd200","ph200")

covar.data <- covariates.species[,model.covariates.vector]
# calculate quadratic polynomial cols
covar.data <- data.frame(poly(covar.data$tempmtcp, 2),
                        poly(covar.data$precipseas, 2),
                        poly(covar.data$precipann, 2),
                        poly(covar.data$tempmtwp, 2),
                        poly(covar.data$rough500, 2),
                        poly(covar.data$bd200, 2),
                        poly(covar.data$ph200, 2))
names(covar.data) <- c("tempmtcp.1","tempmtcp.2",
                      "precipseas.1","precipseas.2",
                      "precipann.1","precipann.2",
                      "tempmtwp.1","tempmtwp.2",
                      "rough500.1","rough500.2",
                      "bd200.1","bd200.2",
                      "ph200.1","ph200.2")


# generate model data
model.data <- data.frame(covariates.species[,model.species.vector], covar.data)

# define model form
RCP.form <- paste0(model.species.string,"~","1","+",model.covariates.string)

# add column containing factor for form.spp
# e.g. model.data$observer <- as.factor(covariates.species$Observers)
model.data$score.method <- as.factor(covariates.species$Species.score.method)
model.data$date.int <- scale(as.integer(covariates.species$Date), center=T, scale=T)

# define species form
# e.g. "~Observer"
species.form <- "~score.method+date.int"

# clear unused variables
rm(covariates.species, species.count)
gc()



# fit mixture models ------------------------------------------------------

# The code below is designed to fit one model, based on the logical pareterisation above
# Since lots of models are fit, just the key information is saved
# The final models, after model seleciton, can then be fit with the saved parameters (without optimisation)
# See below this block for a non-cluster alternative

my.cont <- list(maxit=3000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)

if (species.model %in% c("both", "species")) {
  tic <- proc.time()
  fit.regi <- regimix(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=nRCP, 
                     dist="Bernoulli", control=my.cont, inits="noPreClust", titbits=TRUE)
  toc <- proc.time()
  
  # write model fit stats
  modelStats=list(sites=site.names, covariates=model.covariates.vector, species=model.species.vector, 
                  SppMin=species.n, SppN=fit.regi$S, nRCP=fit.regi$nRCP, runtime=round((toc-tic)[3]/60),
                  AIC=fit.regi$AIC, BIC=fit.regi$BIC, postProbs=fit.regi$postProbs, logl=fit.regi$logl, 
                  coefs=fit.regi$coefs, species.form=species.form, penalties=unlist(my.cont))
  save(modelStats, file=paste0("spec/RegimixStats.n",fit.regi$n,
                               ".rcp",fit.regi$nRCP,".s",fit.regi$S,round(fit.regi$logl),".RData"))
}

if (species.model %in% c("both", "nospecies")) {
  tic <- proc.time()
  fit.regi <- regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=nRCP, 
                     dist="Bernoulli", control=my.cont, inits="noPreClust", titbits=TRUE)
  toc <- proc.time()
  
  # write model fit stats
  modelStats=list(sites=site.names, covariates=model.covariates.vector, species=model.species.vector, 
                  SppMin=species.n, SppN=fit.regi$S, nRCP=fit.regi$nRCP, runtime=round((toc-tic)[3]/60),
                  AIC=fit.regi$AIC, BIC=fit.regi$BIC, postProbs=fit.regi$postProbs, logl=fit.regi$logl, 
                  coefs=fit.regi$coefs, species.form=species.form, penalties=unlist(my.cont))
  save(modelStats, file=paste0("nospec/RegimixStats.n",fit.regi$n,
                               ".rcp",fit.regi$nRCP,".s",fit.regi$S,round(fit.regi$logl),".RData"))
}


# Alternatively, if you think you have subsetted to a small enough number, you can try regimix.multifit()
# That will allow you to do multiple starts (parralellised, if your system allows it) on the same model specs,
# and then choose the best fit. See ?regimix.multifit() for a worked example, but here's an example for this data:

# my.cont <- list(maxit=3000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)
# regi.fits <- regimix.multifit(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=nRCP, 
#                    dist="Bernoulli", control=my.cont, inits="noPreClust", nstart=30, titbits=FALSE)
# # Sometimes the model 'mis-fits' and one or more of the RCP groups has no sites associated
# # with it. These need to be removed (based on the colSums of the posterior probabilities)
# postProbSums <- t(sapply(regi.fits, function(x) colSums(x$postProbs)))
# # Identify those models with zero posterior prob classes
# allGoodUns <- apply(postProbSums, 1, function(x) all(x!=0))
# # Subset the fits
# regi.fits.clean <- regi.fits[allGoodUns]
# # Choose the model with the lowest BIC
# goodUn <- which.min(sapply(regi.fits.clean, BIC))
# # Using the 'best' model, save the key information needed for re-fitting
# modelStats <- list(nRCP=regi.fits.clean[[goodUn]]$nRCP, coef=regi.fits.clean[[goodUn]]$coefs)

# Now you can re-join the workflow in the "fit the best model" section below



# load results ------------------------------------------------------------

# Once you have lots of model fit results, these functions will allow you to quickly load and plot them

regimix.results <- function(path, pattern, controls="", remove.misfit=T, plot.only=T) {
  controls <- controls
  
  files <- list.files(path=paste0(path,controls),full.names=T,pattern=pattern)
  nRCP <- as.numeric(rep(NA,length(files)))
  BIC <- as.numeric(rep(NA,length(files)))
  AIC <- as.numeric(rep(NA,length(files)))
  runtime <- as.numeric(rep(NA,length(files)))
  minPostProb <- as.numeric(rep(NA,length(files)))
  maxPostProb <- as.numeric(rep(NA,length(files)))
  logl <- as.numeric(rep(NA,length(files)))
  
  for (i in 1:length(files)) {
    load(files[i])
    nRCP[i] <- modelStats$nRCP
    BIC[i] <- modelStats$BIC
    AIC[i] <- modelStats$AIC
    runtime[i] <- modelStats$runtime
    logl[i] <- modelStats$logl
    minPostProb[i] <- min(colSums(modelStats$postProbs))
    maxPostProb[i] <- max(colSums(modelStats$postProbs))
    rm(modelStats)
  }
  
  nRCP.plot <- data.frame(nRCP, BIC, AIC, runtime, logl, minPP=round(minPostProb, 3), maxPP=round(maxPostProb, 3))
  nRCP.plot <- nRCP.plot[!is.na(nRCP),]
  
  if (remove.misfit) {
    # Sometimes the model 'mis-fits' and one or more of the RCP groups has no sites associated
    # with it. These need to be removed (based on the colSums of the posterior probabilities)
    nRCP.plot <- nRCP.plot[nRCP.plot$minPP>2,]
  }
  if (plot.only) {
    return(nRCP.plot[,c("nRCP","BIC")])
  } else {
    return(nRCP.plot)
  }
}

plot.regimix.results <- function(regimix.results) {
  print(plot(regimix.results$BIC~regimix.results$nRCP, main=paste0("BIC")))
  print(abline(v=regimix.results$nRCP[regimix.results$BIC==min(regimix.results$BIC)], lty=2))
}

# Once you have plotted and chosen the number of RCPs that minimises BIC, you can directly choose the file
# since the fitting code above saved the log-Lik value directly into the file name, so you can sort by best fit!



# fit the best model ------------------------------------------------------

# Note here that you will still have to have the data loaded from the model fitting stage
# Can do this for model both with and without 

load("RegimixStats.best.fit.model.with.species.dependency.RData") # load the saved model data (or use existing modelStats object)
my.cont <- list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
params <- unlist(modelStats$coefs)
fit.regi.sp <- regimix(form.RCP=RCP.form, form.spp=species.form, data=model.data, nRCP=modelStats$nRCP,
                      dist="Bernoulli", control=my.cont, inits=params, titbits=TRUE)
# save(fit.regi.sp, file="fit.regi.sp.RData") # save if you want to, but it fits very fast without optimisation

load("RegimixStats.best.fit.model.without.species.dependency.RData") # load the saved model data
my.cont <- list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
params <- unlist(modelStats$coefs)
fit.regi.nosp <- regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=modelStats$nRCP,
                      dist="Bernoulli", control=my.cont, inits=params, titbits=TRUE)
# save(fit.regi.nosp, file="fit.regi.nosp.RData") # save if you want to, but it fits very fast without optimisation


# diagnostic plots
# note that if you did end up fitting big models, you MIGHT need to consider another approach here ;)

# residual plots
plot.regimix(fit.regi.sp, type="RQR", fitted.scale="log")
plot.regimix(fit.regi.nosp, type="RQR", fitted.scale="log")
# cooks distance plots
stability.regimix(fit.regi.sp, oosSizeRange=c(1,10,50,100,200,300,400,500), mc.cores=1, times=700)
stability.regimix(fit.regi.nosp, oosSizeRange=c(1,10,50,100,200,300,400,500), mc.cores=1, times=700)



# prepare covariates for spatial prediciton -------------------------------

# load covariates data
covars <- read.dbf("NSW_covars_2km.dbf", as.is=T)

# make covar data and location data
locs <- covars[,c("lat","long")]

# get poly() coefficients from model fit and predict to newdata
model.covars <- covariates.species[,model.covariates.vector]

covars <- data.frame(predict(poly(model.covars$tempmtcp, 2), newdata=covars$tempmtcp),
                    predict(poly(model.covars$precipseas, 2), newdata=covars$precipseas),
                    predict(poly(model.covars$precipann, 2), newdata=covars$precipann),
                    predict(poly(model.covars$tempmtwp, 2), newdata=covars$tempmtwp),
                    predict(poly(model.covars$rough500, 2), newdata=covars$rough500),
                    predict(poly(model.covars$bd200, 2), newdata=covars$bd200),
                    predict(poly(model.covars$ph200, 2), newdata=covars$ph200))
names(covars) <- c("tempmtcp.1","tempmtcp.2",
                  "precipseas.1","precipseas.2",
                  "precipann.1","precipann.2",
                  "tempmtwp.1","tempmtwp.2",
                  "rough500.1","rough500.2",
                  "bd200.1","bd200.2",
                  "ph200.1","ph200.2")

# need to choose the species model levels
covars$score.method <- model.data$score.method[1]
scale.att <- attributes(scale(as.integer(covariates.species$Date), center=T, scale=T))
covars$date.int <- as.numeric(scale(as.integer(as.Date("2015-01-01")),
                                   center=scale.att$`scaled:center`, scale=scale.att$`scaled:scale`))



# boostrap CIs and make predicitons ---------------------------------------
# note that if you did end up fitting big models, you WILL need to consider another approach here ;)

# build bootstrap sampple distribution
fit.regiboot.sp <- regiboot(fit.regi.sp, nboot=1000)
fit.regiboot.nosp <- regiboot(fit.regi.nosp, nboot=1000)

# predict out to geographic covariate space
predicted.sp <- predict.regimix(object=fit.regiboot.sp, object2=fit.regiboot.sp, newdata=covars)
predicted.nosp <- predict.regimix(object=fit.regiboot.nosp, object2=fit.regiboot.nosp, newdata=covars)

# add back geographic location information
predicted.sp.btpreds <- data.frame(locs, predicted.sp)
predicted.sp.low <- data.frame(locs, predicted.sp)
predicted.sp.up <- data.frame(locs, predicted.sp)

predicted.nosp.btpreds <- data.frame(locs, predicted.nosp)
predicted.nosp.low <- data.frame(locs, predicted.nosp)
predicted.nosp.up <- data.frame(locs, predicted.nosp)

# You can now plot these out, for example with rasterFromXYZ() from the {raster} package
