#########################################################################################
### Helper Functions for RCPmod used in Foster, Hill & Lyons 2016 JRSS paper          ###
### Written by N. Hill July 2015.                                                     ###
### Works with RCPmod version 2.88                                                    ###
#########################################################################################



#### poly_data----
# Generate orthogonal polynomials for RCPmod input. Useful to avoid convergence problems
# Save polynomial bases to transform prediction space

poly_data<-function(poly_vars,      #vector of predictor variable names to create orthogonal polynomials. Matched to column names in 'data'
                    degree,         #vector of same length as poly_vars specifying the polynomial degree for each predictor variable in poly_vars
                    id_vars,        #vector of ID variable names (Not transformed). Matched to column names in 'data'
                    sample_vars=NULL, #vector of sampling level variable names (e.g. Gear type). Matched to column names in 'data'
                    species_vars,    #vector of response variables names. Matched to column names in 'data'
                    offset=NULL,     #name of offset variable.  Matched to column names in 'data'
                    data, ...)
  {
  store_polys<-list()
  for(i in 1:length(poly_vars)){
    store_polys[[i]]<-poly(data[,poly_vars[i]], degree=degree[i])
    dimnames(store_polys[[i]])[[2]]<-paste0(poly_vars[i],seq(1:degree[i]))
  }
  names(store_polys)<-poly_vars
  
  rcp_data<-na.omit(cbind(subset(data, select= c(id_vars, sample_vars, offset, species_vars)), 
                          do.call(cbind, store_polys)))
  
  return(list(rcp_data=rcp_data, poly_output=store_polys))
}


### poly_pred_space----
# Transforming prediction space using same basis as orthogonal polynomials used to build models
# Note: only accomodates one sampling term at the moment
# Note: offset isn't actually used in predict function that predicts RCP membership, but will keep as it might be useful to predict expected abundance of species at site.

poly_pred_space<-function(pred_space,             #dataframe containing variables
                          poly_output,            #extracted list of stored polynomial attribtutes from 'poly_data' function
                          offset_val=NULL,        #an offset value. Possibly mean of offset used in model building. Will be logged within function
                          offset_name=NULL,       #name of offset used in RCP model
                          sampling_vals=NULL,     #level of sampling factor for prediction 
                          sampling_name=NULL,     #name of sampling factor used in RCP model
                          sampling_factor_levels=NULL )  #levels of sampling factor used in RCP model
  {
  
  # transform predictors using saved orthogonal polynomial attributes
  pred_polys<-list()
  vars<-names(poly_output)
  for( i in 1: length(vars)){
    pred_polys[[i]]<- predict( poly_output[[i]], pred_space[, names(pred_space) %in% vars[i]]) 
    dimnames(pred_polys[[i]])[[2]]<-dimnames(poly_output[[i]])[[2]]
  }
  pred_polys_df<-as.data.frame(do.call(cbind, pred_polys))
  
  #create offset term
  if(!is.null(offset_val)){
    pred_polys_df$offset<-log(offset_val)
    names(pred_polys_df)[ncol(pred_polys_df)]<-paste0("log(", offset_name, ")")
  }
  
  # create sampling variable. 
  # only accommodates one sampling factor
  if(!is.null(sampling_vals)){
    reps<- length(sampling_vals)
    pred_polys_df$sampling<-factor(sampling_vals,levels=sampling_factor_levels)
    names(pred_polys_df)[ncol(pred_polys_df)]<-sampling_name
  }
  
  return(pred_polys_df)
}



#### Calculate average SD and CI of species abundances in each RCP in each level of sampling factor.
#Note: not generalised beyond HIMI case
# boot_obj= regiboot object

Sp_abund_all<-function(boot_obj) {
  require(tidyr)
  
  #set up coefficient extraction
  taus<-grepl("tau",dimnames(boot_obj)[[2]])
  alphas<-grepl("alpha",dimnames(boot_obj)[[2]])
  gammas<-grepl("gamma",dimnames(boot_obj)[[2]])
  
  autumn_res<-spring_res<-summer_res<-list()
  
  #run loop for each row in boot object
  for(i in 1:dim(boot_obj)[1]){
    
    #extract and reformat coeficients 
    #tau
    temp_tau<-data.frame(value=boot_obj[i,taus])
    temp_tau$species<-sapply(strsplit(rownames(temp_tau),"_"), "[", 1)
    temp_tau$coef<-sapply(strsplit(rownames(temp_tau),"_"), "[", 3)
    
    tau_aut<-spread(temp_tau[,1:3], species, value)
    tau_aut<-rbind( tau_aut[,-1], -colSums( tau_aut[,-1]))
    
    #alpha- OK as is
    temp_alphas<-boot_obj[i,alphas]
    
    #gamma
    temp_gamma<-data.frame(value=boot_obj[i,gammas])
    temp_gamma$species<-sapply(strsplit(rownames(temp_gamma),"_"), "[", 1)
    temp_gamma$coef<-sapply(strsplit(rownames(temp_gamma),"_"), "[", 3)
    
    #not generalised beyond HIMI case
    gamma_spring<-temp_gamma[1:15,1]
    gamma_summer<-temp_gamma[16:30,1]
    
    #calulate autumn values
    lps <- sweep( tau_aut, 2, temp_alphas, "+") 
    #we don't have an offset
    autumn_res[[i]] <- as.matrix(round(exp( lps),2))
    
    #calculate spring values
    lps_spring<-sweep( lps, 2, gamma_spring, "+") 
    spring_res[[i]] <- as.matrix(round(exp( lps_spring),2))
    
    #calculate summer values
    lps_summer<-sweep( lps, 2, gamma_summer, "+") 
    summer_res[[i]] <- as.matrix(round(exp( lps_summer),2))
  }
  
  
  #compile summaries of results
  autumn_summary<-list(mean=apply(simplify2array(autumn_res), c(1,2), mean),
                       sd=apply(simplify2array(autumn_res), c(1,2), sd),
                       lower=apply(simplify2array(autumn_res), c(1,2), function(x) quantile(x, probs=0.025)),
                       upper=apply(simplify2array(autumn_res), c(1,2), function(x) quantile(x, probs=0.975)))
  
  spring_summary<-list(mean=apply(simplify2array(spring_res), c(1,2), mean),
                       sd=apply(simplify2array(spring_res), c(1,2), sd),
                       lower=apply(simplify2array(spring_res), c(1,2), function(x) quantile(x, probs=0.025)),
                       upper=apply(simplify2array(spring_res), c(1,2), function(x) quantile(x, probs=0.975)))
  
  summer_summary<-list(mean=apply(simplify2array(summer_res), c(1,2), mean),
                       sd=apply(simplify2array(summer_res), c(1,2), sd),
                       lower=apply(simplify2array(summer_res), c(1,2), function(x) quantile(x, probs=0.025)),
                       upper=apply(simplify2array(summer_res), c(1,2), function(x) quantile(x, probs=0.975)))
  return(list(autumn=autumn_summary,spring=spring_summary,summer=summer_summary))
}





#### sampling_dotplot---
# Produce dotplot of sampling-level coefficients (and 95% CI) for each species

sampling_dotplot2<-function(best_mod,              # output of regimix function for final model
                            boot_obj,              # output of regiboot function for final model
                            legend_fact,           # levels of categorical sampling variable to plot. Coefficients relative to first level of factor. 
                                                   # So usually 2:n levels(sampling_variable)
                            CI=c(0.025, 0.975),    # confidence interval to plot
                            col="black",           # colour/s of dots and CI lines. Specified in same way as col is usually specified
                            lty=1)                 # lty= line type of CI lines. Specified in same way as lty is usually specified
{
  
  require(lattice)
  # gammas<- paste0("gamma", 1: (length(Species)*length(sampling_names)))
  gammas<-grepl("gamma",dimnames(boot_obj)[[2]])
  temp_dat<-boot_obj[,gammas]
  
  temp<-data.frame(avs=as.numeric(unname(colMeans(temp_dat))),
                   t(apply(temp_dat, 2, quantile, probs=CI)),
                   #sampling_var=rep(sampling_names, each=length(length(best_mod$names$spp))),
                   sampling_var=sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 3),
                   #Species=factor(rep(best_mod$names$spp,length(sampling_names))))
                   Species=factor(sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 1)))
  
  names(temp)[2:3]<-c("lower", "upper")
  temp$Species<-gsub("."," ", temp$Species, fixed=TRUE) #get rid of '.' in species names
  temp$Species<-as.factor(temp$Species) #convert back to factor
  temp$Species <- factor(temp$Species, levels=rev(levels(temp$Species))) 
  
  trellis.par.set(superpose.symbol=list(pch=16,col=col, cex=1.2),
                  superpose.line=list(col="transparent"))
  dotplot(Species ~ avs, groups=sampling_var, data=temp, cols=col, lty=lty, low=temp$lower, high=temp$upper, subscript=TRUE, 
          auto.key=list(space="top", columns=2, cex=1.4, text=legend_fact),
          ylab=list(cex=1.4), xlab=list("Coefficient",cex=1.4),
          scales = list(tck = c(1, 0), x=list(cex=1.2), y=list(cex=1.2)),
          prepanel = function(x, y, ...) { list(xlim=range(temp$lower, temp$upper)) }, 
          panel=panel.superpose,
          panel.groups=function(x, y, subscripts, group.number, cols, low, high, ...)
          {
            if(group.number==1) jiggle <- 0.1 else jiggle <- -0.1
            panel.abline(v=0, lty=2)
            panel.abline(h=1:length(best_mod$names$spp), col.line="light grey", lty=1)
            panel.dotplot(x, y+jiggle, group.number, ...)
            panel.arrows(low[subscripts], y+jiggle, high[subscripts], y+jiggle, code=3, angle=90, 
                         length=0.05, col=cols[group.number], lty=lty[group.number])
            #panel.segments(temp$lower, y+jiggle, + temp$upper, y+jiggle, lty = lty, col =col, lwd=2, cex=1.2)
          })            
  
}



#####predict_maps2_SDF2----
# Plot RCP predictions
# S. Foster version, modified for vertical instead of horizontal RCP layout

predict_maps2_SDF2 <- function(predictions,     #output from predict.regimix
                               pred_space,      #dataframe containing coordinates for the prediction space
                               pred_crop,       #raster of extent of prediction space (used in rasterize function)
                               nRCP,            #the number of RCPs
                               my.ylim=NULL,    #ylimit to plot
                               my.xlim=NULL,    #x limit to plot
                               my.asp=1)        #aspect for plotting
  {
  require(rasterVis)

  preds<-predictions
  av_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootPreds)), pred_crop)
  av_pred <- dropLayer( av_pred, 1)
  low_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootCIs[1:nrow(pred_space),1:nRCP,1])), pred_crop)
  low_pred <- dropLayer( low_pred, 1)
  upp_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootCIs[1:nrow(pred_space),1:nRCP,2])), pred_crop)
  upp_pred <-  dropLayer( upp_pred, 1)
  
  colour <- c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d")#,"#000000")
  breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  
  nPlots <- nRCP*3
  layMat <- matrix( NA, ncol=3+1, nrow=nRCP+2)
  layMat[1,1] <- 0
  layMat[1,-1] <- 1:3
  layMat[-1,1] <- 3+1:(nRCP+1)
  layMat[nrow( layMat),1] <- 0
  tmp <- matrix( 2 + (nrow( layMat)-1) + 1:(3*nRCP), ncol=3, byrow=TRUE)
  layMat[-(c(1,nrow( layMat))),-1] <- tmp
  layMat[nrow( layMat), -1] <- max( layMat, na.rm=TRUE) + 1
  #good grief!
  par( oma=c(1,1,1,1))
  layout( layMat, heights=c(0.15, rep(1, nRCP), 0.5), widths=c(0.5, rep( 1, 3)))
  #first plot (empty)
  
  #second plot (lower CI label)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Lower CI",cex=1.5)
  #third plot (pt pred)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Point Prediction",cex=1.5)
  #fourth plot (upper CI label)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Upper CI",cex=1.5)
  #Next row labels
  for( ii in 1:nRCP){
    par( mar=c( 0,0,0,0))
    plot.new()
    text( 0.5,0.5, paste("RCP", ii, sep=" "), srt=90, cex=1.5)
  }
  #  plot.new()
  
  par( mar=c( 0.7,0.7,0,0)+0.1)
  for( ii in 1:nRCP){
    my.ylab <- c( paste( "RCP ", ii, sep=""), "", "")
    my.yaxt <- c( "s", "n", "n")
    my.xlab <- rep( "", 3)
    if( ii == 1)
      my.xaxt <- c( "n", "n", "n")
    else if( ii==nRCP)
      my.xaxt <- rep( "s", 3)
    else
      my.xaxt <- rep( "n", 3)
    tmpDat <- subset( low_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[1], yaxt=my.yaxt[1], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
    tmpDat <- subset( av_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[2], yaxt=my.yaxt[2], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
    tmpDat <- subset( upp_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[3], yaxt=my.yaxt[3], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
  }
 
  plot.new()
  require( fields)#for image.plot
  image.plot( x=breaks, y=breaks, z=matrix( seq( from=0, to=1, length=length( breaks)^2), ncol=length( breaks)), horizontal=TRUE, legend.only=TRUE, legend.shrink=0.8, col=colour, smallplot=c(0.05,0.95,0.55,0.7), axis.args=list( at=seq( from=-0.055, to=1.055, length=6),#c(-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05),# c(0,0.25,0.5,0.75,1),
                                                                                                                                                                                                                                    labels=seq( from=0, to=1, by=0.2),#c(0,0.25,0.5,0.75,1),
                                                                                                                                                                                                                                    cex.axis=1.1), legend.args=list(text='Probability of RCP membership', side=1, line=1.9, cex=0.8))
  
  return( NULL)
}
