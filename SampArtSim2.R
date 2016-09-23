########################################################################
####	Illustration of effects of sampling artefacts				####
####	Scott Foster	July 2016	Hobart							####
####	Appendical material											####
########################################################################

#setting up where the server should look for packages -- you may not need this line
.libPaths( c(.libPaths(),"~/lib/R/library"))

library( RCPmod)

rm( list=ls())	#cleaning before populating

#how many cores to run the script on (for bootstrapping and multifits)
my.mc.cores <- 14	#will be ignored on windows machines (effectively set to running in serial)

for( seedseed in c(12,112,212,312,412,512,612,712,812,912,1012)){	# an aribtrary set of numbers.  I like 12 today!
	
	cat( "Seed: ", seedseed, "\n")
	#number of sampling sites
	n <- 150
	#number of species
	S <- 20
	#number of RCPs
	K <- 3
	#single covariate (quadratic)
	set.seed( seedseed)
	X <- sort( runif( n=n, min=-1, max=1))
	X <- poly( X, degree=2)
	X <- cbind( 1, sweep( X, MARGIN=2, STATS=apply( X, 2, sd), FUN="/"))	#just to make interpreting coefficients easier.
	colnames( X) <- c("Intercept","lin","quad")
	#sampling artefact
	W <- model.matrix( ~1+factor( sample( c("gear1","gear2"), n, replace=TRUE)))[,-1,drop=FALSE]
	colnames( W) <- "w"
	#Species Intercepts
	alpha <- rnorm( S, 0, 0.5)
	#species x RCP effects 
	tau <- matrix( rnorm( (K-1)*S), nrow=K-1, ncol=S)
	#Covariate effect on RCP group
	beta <- matrix( c( -1,-2,0, 0,0,-2), nrow=K-1, ncol=ncol(X), byrow=TRUE)
	#Sampling artefact effect
	set.seed( seedseed+1)
	gamma <- matrix( rnorm( 1*S, mean=0, sd=5), nrow=S, ncol=1)	#there will only be 1 sampling effect (level 2 of a factor)
	
	nsim <- 5
	gear.effect.multiplier <- seq( from=0, to=1, length=nsim)
	my.form.RCP <- paste( paste( paste(
	'cbind(', paste( paste( 'spp', 1:S, sep=''), collapse=','), sep=''),
		')',sep=''),
		'~lin+quad',sep='')
	my.form.spp <- ~ w
	preds <- fm.final <- list()
	preds.nosppMod <- fm.final.nosppMod <- list()
	set.seed( seedseed+2)
	for( gg in 1:nsim){
		cat( "Simulation No. ",gg,"\n")
		#temper (or exaggerate) sampling effect
		gamma1 <- gear.effect.multiplier[gg] * gamma
		#generate the data for this amount of sample confounding
		dat <- simRCPdata( nRCP=K, S=S, n=n, p.x=ncol( X), p.w=ncol( W), alpha=alpha, tau=tau, beta=beta, gamma=gamma1, X=X, W=W, dist="Bernoulli")
		#fit the RCP models (code taken from example( regimix.multifit)) #only 50 restarts as model should be pretty well identifiable
		tmp <- system.time( fm.final[[gg]] <- regimix( form.RCP=my.form.RCP, form.spp=my.form.spp, data=dat, nRCP=K, dist="Bernoulli", control=list(quiet=TRUE), init=c(alpha, tau, beta, gamma)))
		print( tmp)
		tmp <- system.time( fm.final.nosppMod[[gg]] <- regimix( form.RCP=my.form.RCP, form.spp=NULL, data=dat, nRCP=K, dist="Bernoulli", control=list(quiet=TRUE), init=c(alpha, tau, beta)))
		print( tmp)
		fm.boot <- regiboot( fm.final[[gg]], nboot=1000, mc.cores=my.mc.cores, quiet=TRUE)	
		fm.boot.nosppMod <- regiboot( fm.final.nosppMod[[gg]], nboot=1000, mc.cores=my.mc.cores, quiet=TRUE)
		preds[[gg]] <- predict.regimix( fm.final[[gg]], fm.boot, mc.cores=my.mc.cores)
		preds.nosppMod[[gg]] <- predict.regimix( fm.final.nosppMod[[gg]], fm.boot.nosppMod, mc.cores=my.mc.cores)
		#order the predictions to match the simulated data (and all simulated data sets)
		my.RCPs <- attr( dat, "RCPs")
	}
	
	#The fitting is done. Now for plotting.
	
	#first the models with the artefact terms
	bitmap( paste( paste( "sim1_withTerms_",seedseed,sep=""),".png",sep=""), type="png16m", res=350, height=12, width=8)
	
	par( mfrow=c(nsim, K), mar=c(2.1,2.1,2.1,1.1), oma=c(2,2,2.5,0))
	for( gg in 1:nsim){	#loop over the simulations
		for( kk in 1:K){	#loop over the RCP groups
			my.main <- ""
			my.xlab <- ""
			my.ylab <- ""
			my.pos <- "left"
			#the point predictions
			p.preds <- preds[[gg]]$ptPreds[,kk]	#for the model with sampling
			#the upper CIs
			uCI <- preds[[gg]]$bootCIs[,kk,"upper"]
			#lower CIs
			lCI <- preds[[gg]]$bootCIs[,kk,"lower"]
			#plotting it all
			plot( c(X[,2],rev(X[,2])), c(lCI,rev(uCI)), type='n', ylab=my.ylab, xlab=my.xlab, main=my.main, ylim=c(0,1))
			polygon( x=c(X[,2],rev(X[,2])), y=c(lCI,rev(uCI)), col=grey(0.75), border=NA)	#the confidence region
			lines( X[,2], p.preds, lwd=2)	#the point predictions
			lines( X[,2], attr( dat, "pis")[,kk], col=2, lwd=2, lty=3)
			if( kk==1) 
				mtext( bquote( rho==.(round( gear.effect.multiplier[gg],2))), side=2, line=2, outer=FALSE, adj=0.5, cex=1.25)
			if( gg==1)
				mtext( paste("RCP", kk,sep=" "), side=3, line=0, outer=FALSE, adj=0.5)
		}
	}
	mtext( "Sampling Artefacts Included in Model", side=3, outer=TRUE, adj=0.5, cex=1.9)
	dev.off()
	
	#now the models *without* the artefact terms
	bitmap( paste( paste( "sim1_noTerms_",seedseed,sep=""),".png",sep=""), type="png16m", res=350, height=12, width=8)
	
	par( mfrow=c(nsim, K), mar=c(2.1,2.1,2.1,1.1), oma=c(2,2,2.5,0))
	for( gg in 1:nsim){	#loop over the simulations
		for( kk in 1:K){	#loop over the RCP groups
			my.main <- ""
			my.xlab <- ""
			my.ylab <- ""
			my.pos <- "left"
			#the point predictions
			p.preds <- preds.nosppMod[[gg]]$ptPreds[,kk]	#for the model with sampling
			#the upper CIs
			uCI <- preds.nosppMod[[gg]]$bootCIs[,kk,"upper"]
			#lower CIs
			lCI <- preds.nosppMod[[gg]]$bootCIs[,kk,"lower"]
			#plotting it all
			plot( c(X[,2],rev(X[,2])), c(lCI,rev(uCI)), type='n', ylab=my.ylab, xlab=my.xlab, main=my.main, ylim=c(0,1))
			polygon( x=c(X[,2],rev(X[,2])), y=c(lCI,rev(uCI)), col=grey(0.75), border=NA)	#the confidence region
			lines( X[,2], p.preds, lwd=2)	#the point predictions
			lines( X[,2], attr( dat, "pis")[,kk], col=2, lwd=2, lty=3)
			if( kk==1) 
				mtext( bquote( rho==.(round( gear.effect.multiplier[gg],2))), side=2, line=2, outer=FALSE, adj=0.5, cex=1.25)
			if( gg==1)
				mtext( paste("RCP", kk,sep=" "), side=3, line=0, outer=FALSE, adj=0.5)
		}
	}
	mtext( "Sampling Artefacts NOT Included in Model", side=3, outer=TRUE, adj=0.5, cex=1.9)
	dev.off()
	
}

rm( list=ls())
