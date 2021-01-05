##############################################################################################################
####################################logistic model fit########################################################
##############################################################################################################

### NAs can be included in x

logismodfit <- function(x,years=1:length(x),sd=NULL,catch=rep(0,length(x)),cov1=rep(0,length(x)),cov2=rep(0,length(x)),
		cov3=rep(0,length(x)),nboot=0,nsim=0,name="",plot=FALSE,hessian=FALSE){
# x is population size
# sd is demographic variance
# catch: harvest after census
# census before recruitment

  ### log-likelihood function
	log.gauss <- function(x, my, sig2){

		log.gauss <- -0.5*log(sig2)-0.5*((x-my)^2)/sig2
		return(log.gauss)
	}

	###########################################
	######### likelihood function #############
	###########################################
	lik <- function(par){
		br <- par[1] ### coef for s (i.e. r - 0.5 sigma.e2)
		bk <- par[2] ### coef for beta
		b1 <- par[3] ### coef for cov1
		b2 <- par[4] ### coef for cov2
		b3 <- par[5] ### coef for cov3
		sig2 <- exp(par[6]) # env var
		my <- xc - 0.5*sd*exp(-xc) + br + bk*exp(xc) + b1*cov1 + b2*cov2 + b3*cov3 ### expected population size
		s <- sig2 + sd*exp(-xc) ### expected variance on the log scale
		a <- -sum(log.gauss(x[-1], my[-length(xc)], s[-length(xc)]), na.rm = TRUE)
		return(a)
	}
	############################################

	### Log-likelihood optimization
	estparam <- function(par,hessian=FALSE){
		test <- optim(par,lik, control=list(trace=F,maxit=3000),hessian=hessian)
		return(test)
	}



	#######################################################################################################
	##########     set initial values for optimization                                        #############
	#######################################################################################################

	### Harvesting
	if (any(na.omit(x)<=0)) stop("N less than or equal to 0!!!")
	if(any(diff(years)!=1)){
		temp <- cbind(years,x,catch,cov1,cov2,cov3)
		years <- data.frame(years=min(years):max(years))
		test <- merge(years,temp,all.x=TRUE)
		years<-test$years; x<-test$x; catch<-test$catch; cov1<-test$cov1; cov2<-test$cov2; cov3<-test$cov3
	}
	xc <- log(x-catch) # population size after harvest
	x <- log(x)
	xback <- x

	### Set starting values
	fit <- lm(diff(x)~exp(x[-length(x)]))
	est.r <- coef(fit)[1]
	est.k <- -coef(fit)[1]/coef(fit)[2]

	par <- vector("numeric",6)
	par[1] <- est.r
	par[2] <- -est.r/est.k
	par[3] <- 0
	par[4] <- 0
	par[5] <- 0
	par[6] <- log(0.074118)

	#########select part of covariates in use #############
	NAS <- cbind(x,cov1,cov2,cov3)
	NAS <- is.na(NAS)
	sel <- c((apply(NAS,1,sum)==0),FALSE)
	sel[length(sel)] <- FALSE
	if (sum(sel)<5) stop("too few cases")

	# Centering covariates
	cov1 <- cov1-mean(cov1[sel])
	cov2 <- cov2-mean(cov2[sel])
	cov3 <- cov3-mean(cov3[sel])

	cov1boot <- cov1
	cov2boot <- cov2
	cov3boot <- cov3
	cov1boot[is.na(cov1boot)] <- 0
	cov2boot[is.na(cov2boot)] <- 0
	cov3boot[is.na(cov3boot)] <- 0

	ncov <- sum(c(any(na.omit(cov1)!=0),any(na.omit(cov2)!=0),any(na.omit(cov3)!=0)))

	##################################


	########################################################################################################



	########################################################################################################
	#########################      estimation of parameters        #########################################
	########################################################################################################

	### Uncertainty estimation using hessian
	N <- exp(x); names(N) <- as.character(years)
	estimat <- estparam(par,hessian=hessian)
	hessval <- NA
	if (hessian & ncov==1) hessval <- sqrt(diag(solve(estimat$hessian[-(4:5),-(4:5)])))[3]
	if (hessian & ncov==2) hessval <- sqrt(diag(solve(estimat$hessian[-5,-5])))[3:4]
	if (hessian & ncov==3) hessval <- sqrt(diag(solve(estimat$hessian)))[3:5]

	### Returning outputs
	br <- estimat$par[1]
	bk <- estimat$par[2]
	b1.est <- estimat$par[3]
	b2.est <- estimat$par[4]
	b3.est <- estimat$par[5]
	sig2 <- exp(estimat$par[6])   # this is the environmental variance component that is not due to covariates
	npar <- 3+sum(c(any(na.omit(cov1)!=0),any(na.omit(cov2)!=0),any(na.omit(cov3)!=0)))
	nobs <- sum(sel)
	aic <- 2*estimat$value+ 2*npar
	aicc <- 2*estimat$value+2*npar+2*npar*(npar+1)/(nobs-npar+1)
	r <- br+b1.est*mean(cov1[sel])+b2.est*mean(cov2[sel])+b3.est*mean(cov3[sel])
#r <- br+b1.est*mean(cov1)+b2.est*mean(cov2)+b3.est*mean(cov3)
	k <- -r/bk
	predicted <- c(NA,xc-0.5*sd*exp(-xc)+br+bk*exp(xc)+b1.est*cov1+b2.est*cov2+b3.est*cov3)
	residual <-  c(x,NA)-predicted
#demcomp <- c(sd*exp(-x)[-length(x)+1],NA)
	demcomp <- sd*exp(-xc)
	predicted <- predicted[-length(predicted)]
	residual <- residual[-1]
#demcomp <- demcomp[-1]
	names(demcomp) <- names(residual) <- names(predicted) <- as.character(years)
	#########################################################################################################


	#######################################################################################################
	##########              simulate data based on estimated parameters                       #############
	#######################################################################################################
# adjust for harvest after census

	### Simulate time series repeatedly using the estimated parameters
	### to bootstrap the uncertainty

	simdat <- function(b1=b1.est,b2=b2.est,b3=b3.est){
		fir <- min(which(!is.na(x)))
		xsim <- rep(NA,length(x))
		xsim[fir] <- x[fir] # initial value, census first year before hunting
		tall <- 0
		CH <- catch
		if (any(is.na(CH))) CH[which(is.na(CH))] <- round(approx(x=1:length(CH),y=CH,xout=which(is.na(CH)))$y)

		while (any(is.na(xsim))| any(xsim==0)){
			tall <- tall+1
			if (tall>1000) stop("unable to simulate time series")
			for (i in (fir+1):length(x)){
				AH <- log(exp(xsim[i-1])-CH[i-1]) # pop after hunting
				if (!is.finite(AH) | AH<0) AH <- 0
				xsim[i] <- AH-0.5*sd*exp(-AH)+br+bk*exp(AH)+b1*cov1boot[i-1]+b2*cov2boot[i-1]+
						b3*cov3boot[i-1]+sqrt(sig2+sd*exp(-AH))*rnorm(1)
				if (xsim[i]<0) xsim[i] <- 0
			}
		}
		xsim[is.na(x)] <- NA
		return(xsim)
	}
	######################################################################################################


	########################################################################################################
	###############################    parametric bootstrapping       ######################################
	########################################################################################################

	### Apply bootstrapping procedure on simulated data

	bootres <- NULL
	x <- xback
	if (nboot>0){
		sim.matrix <- matrix(NA,length(x),nboot)
		for (i in 1:nboot){
			sim.matrix[,i] <- simdat()
		}

		bootres <- matrix(NA,nboot,6)
		colnames(bootres) <- c("K","r","sig2","b1","b2","b3")
		catval <- seq(0,10000,100)
		for (i in 1:nboot){
			if (any(i==catval)) cat("boot",i, "of", nboot, "\n")
			x <- xc <- sim.matrix[,i]
			estimat <- estparam(par)
			#br <- estimat$par[1]
			#bk <- estimat$par[2]
			bootres[i,4] <- estimat$par[3]
			bootres[i,5] <- estimat$par[4]
			bootres[i,6] <- estimat$par[5]
			bootres[i,3] <- round(exp(estimat$par[6]),7)
			bootres[i,2] <- estimat$par[1]+estimat$par[3]*mean(cov1[sel])+estimat$par[4]*mean(cov2[sel])+estimat$par[5]*mean(cov3[sel])
			bootres[i,1] <- -estimat$par[1]/estimat$par[2]
		}
	}
	########################################################################################################


	#########################################################################################################
	#####################   simulation test of covariates      ##############################################
	#########################################################################################################

	### Simulate time series repeatedly using the estimated parameters
	### while regression coefficients of covariates are set to 0.
	### Compare this to null expectation. Can the chance of the regression
	### coefficient being non-zero be due to chance or is this a real effect?

	p.cov <- NULL
	if (nsim>0){
		catval <- c(1,seq(0,10000,100))

		#ncov <- sum(c(any(cov1!=0),any(cov2!=0),any(cov3!=0)))
		simres <- matrix(NA,nsim,3)

		if (ncov>0){
			if (any(na.omit(cov1)!=0)){
				sim.matrix <- matrix(NA,length(x),nsim)
				x <- xback
				for (i in 1:nsim){
					sim.matrix[,i] <- simdat(b1=0)
				}

				for (i in 1:nsim){
					if (any(i==catval)) cat("sim",i, "of", nsim, "covariate 1 of ",ncov,  "\n")
					x <- xc <- sim.matrix[,i]
					estimat <- estparam(par)
					simres[i,1] <- estimat$par[3]
				}
			}

			if (any(na.omit(cov2)!=0)){
				sim.matrix <- matrix(NA,length(x),nsim)
				x <- xback
				for (i in 1:nsim){
					sim.matrix[,i] <- simdat(b2=0)
				}


				for (i in 1:nsim){
					if (any(i==catval)) cat("sim",i, "of", nsim, "covariate 2 of ",ncov,  "\n")
					x <- xc <- sim.matrix[,i]
					estimat <- estparam(par)
					simres[i,2] <- estimat$par[4]
				}
			}

			if (any(na.omit(cov3)!=0)){
				sim.matrix <- matrix(NA,length(x),nsim)
				x <- xback
				for (i in 1:nsim){
					sim.matrix[,i] <- simdat(b3=0)
				}

				for (i in 1:nsim){
					if (any(i==catval)) cat("sim",i, "of", nsim, "covariate 3 of ",ncov,  "\n")
					x <- xc <- sim.matrix[,i]
					estimat <- estparam(par)
					simres[i,3] <- estimat$par[5]
				}
			}

		}
		p.cov <- matrix(NA,3,3)
		dimnames(p.cov) <- list(c("b1","b2","b3"),c("p.less", "p.greater", "p.diff0"))

		p.cov[1,1] <- sum(simres[,1]<b1.est)/nsim
		p.cov[1,2] <- sum(simres[,1]>b1.est)/nsim
		p.cov[1,3] <- sum(abs(simres[,1])>abs(b1.est))/nsim

		p.cov[2,1] <- sum(simres[,2]<b2.est)/nsim
		p.cov[2,2] <- sum(simres[,2]>b2.est)/nsim
		p.cov[2,3] <- sum(abs(simres[,2])>abs(b2.est))/nsim

		p.cov[3,1] <- sum(simres[,3]<b3.est)/nsim
		p.cov[3,2] <- sum(simres[,3]>b3.est)/nsim
		p.cov[3,3] <- sum(abs(simres[,3])>abs(b3.est))/nsim
	}
	#########################################################################################################
	##### gives actually used vector of the covariate
	cov1[!sel] <- NA
	cov2[!sel] <- NA
	cov3[!sel] <- NA
	x <- xback

	res <- list("s"=r, "K"=k, "sig2"=sig2, "sigd"=sd,"b1"=b1.est, "b2"=b2.est, "b3"=b3.est, "N"=N,"catch"=catch,
			"cov1"=cov1, "cov2"=cov2, "cov3"=cov3, "predicted"=predicted, "residuals"=residual, "demcomp"=demcomp,
			"boot"=bootres, "p.cov"=p.cov,"aic"=aic,"aicc"=aicc,hessval=hessval,"model"="logistic","name"=name)
	if (plot) logistic.plot(res)
	return(res)
}


##############################################################################################################
##############################################################################################################
##############################################################################################################


portman.Q <- function (x, K) {
	# portman.Q uses the cummulative ACF to test for whiteness  of a time series.
	# This is the Ljung-Box version of the the Portemanteau  test for
	# whiteness (Tong 1990). It may in particular be  usefull to test
	# for whiteness in the residuals from time  series models.
	#
	# A vector is returned consisting of the asymtpotic  chi-square
	# value, the associated d.f. and asymptotic  p.val for the test of
	# whiteness. p.val<0.05 -> non-whiteness!
	# Tong, H. (1990) Non-linear time series : a dynamical  system approach. Clarendon Press, Oxford.
	# Author: Ottar N. Bjornstad onb1@psu.edu

	Q <- 0
	n <- length(x)
	p <- acf(x, plot = FALSE, lag.max = K,na.action=na.pass)$acf[2:(K + 1)]
	for (k in 1:K) Q <- Q + p[k]^2/(n - k)
	Q <- n * (n + 2) * Q
	res <- list(chisq = round(Q, 4), df = K, p.val = round(1 -
							pchisq(Q, K), 4))
	unlist(res)
}



######################################################################################################################################################################
######################################################################################################################################################################
logistic.plot <- function(fit){
	if(fit$model!="logistic") stop("The fit is not based on the logistic model!")
	lok <- paste(fit$name)
	year <- as.numeric(names(fit$N))
	obs.growth <- log(c(fit$N[-1],NA))-log((fit$N-fit$catch))
	Nmax <- max(c(fit$K,fit$N),na.rm=T); if (Nmax>2*max(fit$N,na.rm=T)) Nmax <- 1.2*max(fit$N,na.rm=T)
	Nmin <- min(fit$N, na.rm=T)
	Nrange <- seq(trunc(Nmin*0.9),trunc(Nmax*1.1),by=0.1);  Nrange <- Nrange[Nrange>0.9]
	pred <- NA
	pred <- -0.5*fit$sigd*exp(-log(Nrange))+fit$s+(-fit$s/fit$K)*Nrange


	is.R <- !is.null(version$language)  # Check for R/S-Plus
	if (is.R) {     # R
		portman.Q <- function (x, K) {
			# portman.Q uses the cummulative ACF to test for whiteness  of a time series.
			# This is the Ljung-Box version of the the Portemanteau  test for
			# whiteness (Tong 1990). It may in particular be  usefull to test
			# for whiteness in the residuals from time  series models.
			#
			# A vector is returned consisting of the asymtpotic  chi-square
			# value, the associated d.f. and asymptotic  p.val for the test of
			# whiteness. p.val<0.05 -> non-whiteness!
			# Tong, H. (1990) Non-linear time series : a dynamical  system approach. Clarendon Press, Oxford.
			# Author: Ottar N. Bjornstad onb1@psu.edu

			Q <- 0
			n <- length(x)
			p <- acf(x, plot = FALSE, lag.max = K,na.action=na.pass)$acf[2:(K + 1)]
			for (k in 1:K) Q <- Q + p[k]^2/(n - k)
			Q <- n * (n + 2) * Q
			res <- list(chisq = round(Q, 4), df = K, p.val = round(1 -
									pchisq(Q, K), 4))
			unlist(res)
		}

		par(mfrow=c(4,2),mai=c(.7,.7,.3,.2))

		plot(Nrange[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))],pred[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))],type="n",ylab="predicted/observed r",xlab="N",
				ylim=range(c(obs.growth,pred[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))]),na.rm=T),
				sub=paste("s: ", round(fit$s,3),"  K: ", round(fit$K,1), "  (sig_e)^2: ", round(fit$sig2,3)),cex.sub=0.7)
		lines(Nrange[is.finite(pred)],pred[is.finite(pred)],col="red",lwd=2)
		points(fit$N,obs.growth)
		sel <- which((abs(scale(fit$residuals))>1.5))
		if (length(sel)>0) text(fit$N[sel],obs.growth[sel],labels=year[sel],pos=4,cex=0.7)
		if (!is.null(lok)) title(lok)

		plot(year,log(fit$N),type="n",ylim=range(c(fit$predicted,log(fit$N)),na.rm=T)-c(0.1,-0.1),ylab="observed/predicted ln N")
		if (fit$K<max(fit$N,na.rm=T))abline(h=log(fit$K),lty=3,col="red",lwd=1)
		points(year,log(fit$N))
		lines(year,log(fit$N))
		points(year,fit$predicted,col="red",pch=16)

		plot(year,obs.growth,type="n",ylim=range(c(obs.growth,diff(fit$predicted)),na.rm=T),xlim=range(year),ylab="observed/predicted r")
		points(year,obs.growth)
		lines(year,obs.growth)
		points(year,-0.5*fit$sigd*exp(-log(fit$N-fit$catch))+fit$s+(-fit$s/fit$K)*(fit$N-fit$catch),col="red",pch=16)
		abline(h=0,lty=2)

		plot(year,fit$residuals,ylab="residuals")
		lines(year,fit$residuals)
		abline(h=0,lty=2)

		lmfit <- lm(fit$residuals~year,na.action="na.exclude")
		prob <- round(anova(lmfit)[1,5],3)
		plot(year,fit$residuals,ylab="residuals",sub=paste("p-val: ",round(prob,3)))
		if (prob<0.05) abline(lmfit,lwd=2,col="red") else abline(lmfit,lwd=1,lty=2)

		lmfit <- lm(fit$residuals~fit$N,na.action="na.exclude")
		prob <- round(anova(lmfit)[1,5],3)
		plot(fit$N,fit$residuals,ylab="residuals",xlab="N",sub=paste("p-val: ",round(prob,3)))
		if (prob<0.05) abline(lmfit,lwd=2,col="red") else abline(lmfit,lwd=1,lty=2)

		acf(fit$residuals,na.action=na.pass,lag.max=10,main="")
		port <- portman.Q(fit$residuals,K=10)
		plot(c(1,2),c(1,2),type="n", axes=F,ylab="",xlab="")
		text(1.1,1.8,"Test for whiteness of residuals:",cex=1,pos=4)
		text(1.1,1.6,paste("p-val:", port[3],sep=" "),pos=4)
		if (!any(is.na(port))){
			if (port[3]>0.05) text(1.1,1.4,"Residuals are not significantly different",pos=4)
			if (port[3]>0.05) text(1.1,1.3,"from white noise",pos=4)
			if (port[3]<=0.05)text(1.1,1.4,"WARNING: Residuals are significantly",pos=4,col="red",cex=0.7)
			if (port[3]<=0.05)text(1.1,1.3,"different from white noise",pos=4,col="red",cex=0.7)
		}
	}
}






##############################################################################################################
#####################              thetalogistic model                 #######################################
##############################################################################################################

### If missing values, mark those using "year"
# So:
# Year: 2000, 2001, 2002, 2003
# N: 100, 150, NA, 170

# Add
# N <- c(100, 150, 170)
# year <- c(2000, 2001, 2003)


thetalogistic.est <- function(N,year,sd=NULL,nbias=10000,nboot=0,ncorr=0,
		theta=NULL,thetamin=-2,thetamax=20,name=""){

	############################ INPUT ##################################################
# N is the vector of population sizes
# year is the vector of years
# NB: only timesteps of length 1 is included in the estimation procedure !!!!
# sd is the demographic variance
# nbias is number of simulations for correcting bias in r1
# nboot is the number of bootstraps
# ncorr is the number of simulations for correcting bias in r1 within each bootstrap
# theta=NULL -> estimate theta and other parameters,
# theta=value -> estimate other parameters based on theta=value
# thetamin and thetamax are min and max values for theta when optimizing,
# extreme values could lead to optimization problems
#

	########################### OUTPUT #################################################
# - predicted[year t] is predicted growth rate (r) for [year t]->[year t +1]
# - residual[year t] is residual for [year t]->[year t +1]:
#   observed (ln(N_t+1)-ln(N_t))-predicted adjusted for (1/(2*N))*(sigma_d)^2
# - demcomp[year t] is demographic component for [year t]->[year t +1]:
#   (sigma_d)^2*(1/N)

	g <- function(x,theta){
		if (theta!=0) out <- (exp(x*theta)-1)/theta else out <- x
		return(out)
	}

	est <- function(x,theta,r=NULL){     # only timesteps of length 1 are included in the estimation
		uz <- u2 <- v2 <- vz <- uv <- se <- 0
		u <- sqrt(diff(year))
		z <- diff(x)/u
		v <- u*g(x[-length(x)],theta)
		v[is.na(z)] <- NA
		u[is.na(z)] <- NA # s?rg for at u'ene bare er de som inng?r
		uz <- sum(u*z,na.rm=T)
		u2 <- sum(u^2,na.rm=T)
		v2 <- sum(v^2,na.rm=T)
		vz <- sum(v*z,na.rm=T)
		uv <- sum(u*v,na.rm=T)
		b <- (uz*uv-vz*u2)/(uv*uv-u2*v2) ### bias in r1
		if (is.null(r)) r <- (uz-b*uv)/u2
		c <- uv/u2 ### unknown, perhaps the top point of the curve of delta N | N versus N --- see Fig 5.1 in The Book
		nstep <- length(na.omit(z))
		xsel <- x[-length(x)]; xsel[is.na(z)] <- NA
		# varierer mellom n-3 og n-2 i f?rste ledd i div Steinar-skript:
		se <- sum( (z-r*u-b*v)^2/(nstep-2) -sd*exp(-xsel)/(nstep-1) ,na.rm=T)
		if (theta==0) {  K <- exp(-r/b); if (K>1000000) K <- 1000000}
		if (theta!=0) {
			TF <- ((1-r/b*theta)<0.00000001)
			if(TF) K <- 1000000
			if (!TF) K <- exp(log(1-r/b*theta)/theta)
		}
		gamma <- r*theta/(1-exp(-theta*log(K)))
		predicted <- (r+b*v)                               # predicted growth rate
		predicted.adj <- predicted - 0.5*sd*exp(-xsel)     # diffusion approximation
		residuals <- c(z-predicted.adj,NA)
		# residuals: observed (ln(N_t+1)-ln(N_t))-predicted adjusted for (1/(2*N))*(sigma_d)^2
		predicted.adj <- c(predicted.adj,NA)
		demcomp <- c(sd*exp(-xsel),NA)
		return(list(b=b,r1=r,c=c,se=se,K=K,gamma=gamma,theta=theta,
						predicted=predicted.adj,residuals=residuals,demcomp=demcomp))
	}



	sim.corr <- function(nbias,res){
		sim.mat <- simdata.corr(nbias,res)
		if (is.null(sim.mat)) b <- NA
		if (!is.null(sim.mat)){
			sim.mat[is.na(x),] <- NA
			xmat <- sim.mat
			zmat <- diff(sim.mat)
			vmat <- g(xmat[-nrow(xmat),],res$theta)
			vmat[is.na(diff(x)),1:ncol(vmat)] <- NA
			uz <- colSums(zmat,na.rm=TRUE)
			u2 <- rep(nrow(zmat),ncol(zmat))
			v2 <- colSums(vmat^2,na.rm=TRUE)
			vz <- colSums(vmat*zmat,na.rm=TRUE)
			uv <- colSums(vmat,na.rm=TRUE)
			b <- (uz*uv-vz*u2)/(uv*uv-u2*v2)
		}
		return(b)
	}

	simdata.corr <- function(nbias,res){
		nr.comp.ts <- 0
		fac <- 1.05
		unable <- FALSE
		while (nr.comp.ts<nbias){
			sim.mat <- matrix(0,n,trunc(nbias*fac))
			sim.mat[1,] <- x[1]
			norm.mat <-  matrix(rnorm((n-1)*trunc(nbias*fac)),nrow=(n-1))
			for (i in 2:n){
				sel <- which(sim.mat[i-1,]>0 & is.finite(sim.mat[i-1,]))
				nsel <- length(sel)
				sim.mat[i,sel] <- sim.mat[i-1,sel] +
						(res$r1+res$b*g(sim.mat[i-1,sel],res$theta)) +
						sqrt(res$se+sd*exp(-sim.mat[i-1,sel]))*norm.mat[i-1,sel]
			}
			nr.comp.ts <- sum(sim.mat[n,]>0 & is.finite(sim.mat[n,]))
			fac <- 1/(nr.comp.ts/nbias); fac <- fac + 0.1*fac
			if (fac > 500) {unable <- TRUE;break }
		}
		if (!unable){
			sim.mat <- sim.mat[,(sim.mat[n,]>0 & is.finite(sim.mat[n,]))]
			sim.mat <- sim.mat[,1:nbias]
		}
		if (unable){ sim.mat <- NULL}
		return(sim.mat)
	}

	sim.boot <- function(nboot,ncorr){
		sim.mat <- matrix(NA,n,nboot)
		cateval <- seq(0,10000,by=100)
		cat("simulating time series","\n")
		for (i in 1:nboot){ sim.mat[,i] <- simdata.boot() }
		sim.mat[is.na(x),] <- NA
		cat("fitting parameters","\n")
		# theta is constant in bootstraps !!!
		simres <- apply(sim.mat,2,function (xsim) est(xsim,res$theta))
		if (ncorr>0){
			cat("bias correction of bootstraps (r1) ","\n")
			for (i in 1:length(simres)){
				if (any(cateval==i)) cat("correction ", i, " of ", length(simres),"\n")
				bias <- mean(sim.corr(ncorr,simres[[i]]))-simres[[i]]$b
				simres[[i]]$r1.corr <- simres[[i]]$r1+bias*simres[[i]]$c
			}
		}
		return(simres)
	}

	simdata.boot <- function(){
		d <- diff(year)
		xx <- rep(0,n)
		xx[1] <- x1
		count <- 0
		if (res$theta==0) b <- -res$r1.corr/log(res$K) else{
			b <- -res$r1.corr*res$theta/(exp(res$theta*log(res$K))-1)
		}

		while (xx[length(xx)]==0){
			for (i in 2:n){
				if (xx[i-1]>0){
					xx[i] <- xx[i-1] +
							(res$r1.corr+b*(exp(xx[i-1]*res$theta)-1)/res$theta) *
							d[i-1] +sqrt((res$se+sd*exp(-xx[i-1]))*d[i-1]) *
							rnorm(1)
				}
			}
			count <- count + 1; if (count==1000) stop("unable to simulate ts-series!")
		}
		return(xx)
	}

	func <- function(par,x){
		res <- est(x,par)
		res$se
	}

	est.theta <- function(x,thetamin,thetamax){
		optimize(func,lower=thetamin, upper=thetamax,x=x)$minimum
	}

	f.stasj.distr.thetalog <- function(n,K,theta,sige2,r1){
		r <- r1/(1-K^-theta)
		alpha <- (2*r/sige2)-1
		f <- exp(((log(abs(theta)) +(alpha/theta)*log((alpha+1)/(theta))) -
							(log(K)+lgamma(alpha/theta))) + (alpha-1)*(log(n)-log(K)) -
						((alpha+1)/theta)*(n/K)^theta)
		f <- f/sum(f)
		return(f)
	}

	################ initial data handling and transformation ############


	if (any(diff(year)!=1)){                # if some rows (years) are missing, insert missing years with NA for N
		yr <- data.frame(year=min(year):max(year))
		newdata <- merge(yr,data.frame(year=year,N=N),all.x=T)
		N <- newdata$N ; year <- newdata$year }

	x <- log(N)
	n <- length(x)       # number of cases
	# copy population size in the first year; starting point for simulating time series
	x1 <- x[1]
	sd <- sd


	################ parameter estimation ################################
	res <- vector("list",1)
	# estimate theta within the range [thetamin,thetamax]
	if (is.null(theta)) res$theta <- est.theta(x,thetamin,thetamax) else res$theta <- theta
	# estimates of other parameters based on estimate of theta
	res <- est(x,theta=res$theta)
	# name residuals with years, example: residual[year t] is residual for [year t]->[year t +1]
	names(res$residuals) <- year
	# name residuals with years, example: residual[year t] is residual for [year t]->[year t +1]
	names(res$demcomp) <- year
	names(res$predicted) <- year
	res$N <- N; names(res$N) <- year
	res$sd <- sd
	cat("theta = ", round(res$theta,3),"\n")
	cat("r1 = ", round(res$r1,3),"\n")
	cat("gamma = ", round(res$gamma,3),"\n")
	cat("K = ", round(res$K,3),"\n")
	cat("se = ", round(res$se,3),"\n")

	################ bias correction for r1 ################################
	if (nbias>0) bias <- mean(sim.corr(nbias,res))-res$b else bias <- 0
	res$r1.corr <- res$r1+bias*res$c
	cat("biascorrected r1 = ", round(res$r1.corr,3), " nsim: ",nbias,"\n")

	#########   compute stationary distribution, assumes demvar = 0 ! ######
	res$f.stat.dist <- res$var.stat.dist <- NA
	if ((2*res$r1/res$se-1)>0){
		f.stationary <- f.stasj.distr.thetalog(1:round(res$K*5),K=res$K,theta=res$theta,sige2=res$se,r1=res$r1)
		if (any(!is.na(f.stationary))){
			if (any(f.stationary!=0)) f.stationary <- f.stationary/sum(f.stationary) else f.stationary <- NULL
			if (!is.null(f.stationary)){
				stat.dist <- sample(1:round(res$K*5),size=1e+6,replace=TRUE,prob=f.stationary)
				res$f.stat.dist <- f.stationary[1:max(stat.dist)]
				res$var.stat.dist <- var(stat.dist)
			}
		}
	}

	################    bootstrapping       ################################

	if ((nboot>0) & (!is.na(res$r1.corr))){
		cat("bootstrapping","\n")
		bootres <- sim.boot(nboot,ncorr)
		res$boot.K <- unlist(lapply(bootres,function(x) x$K))
		res$boot.r1 <- unlist(lapply(bootres,function(x) x$r1))
		res$boot.r1.corr <- unlist(lapply(bootres,function(x) x$r1.corr))
		res$boot.gamma <- unlist(lapply(bootres,function(x) x$gamma))
		res$boot.se<- unlist(lapply(bootres,function(x) x$se))
	}
	if ((nboot>0) & (is.na(res$r1.corr))){
		cat("bootstrapping","\n")
		res$boot.K <- NA
		res$boot.r1 <- NA
		res$boot.r1.corr <- NA
		res$boot.gamma <- NA
		res$boot.se<- NA
	}
	res$model <- "thetalogistic"
	res$name <- name
	return(res)
}

### residuals ->  sigmaE2 + sigmaD2/N
### demcomp -> sigmaD2/N
### var.stat.dist -> variance of stationary distribution
### f.stat.dist -> probability distribution of population sizes

################################################################################
################################################################################

portman.Q <- function (x, K) {
	# portman.Q uses the cummulative ACF to test for whiteness  of a time series.
	# This is the Ljung-Box version of the the Portemanteau  test for
	# whiteness (Tong 1990). It may in particular be  usefull to test
	# for whiteness in the residuals from time  series models.
	#
	# A vector is returned consisting of the asymtpotic  chi-square
	# value, the associated d.f. and asymptotic  p.val for the test of
	# whiteness. p.val<0.05 -> non-whiteness!
	# Tong, H. (1990) Non-linear time series : a dynamical  system approach. Clarendon Press, Oxford.
	# Author: Ottar N. Bjornstad onb1@psu.edu

	Q <- 0
	n <- length(x)
	p <- acf(x, plot = FALSE, lag.max = K,na.action=na.pass)$acf[2:(K + 1)]
	for (k in 1:K) Q <- Q + p[k]^2/(n - k)
	Q <- n * (n + 2) * Q
	res <- list(chisq = round(Q, 4), df = K, p.val = round(1 -
							pchisq(Q, K), 4))
	unlist(res)
}
################################################################################
################################################################################

thetaplot <- function(f, title=NULL){
	fit <- f
	year <- as.numeric(names(fit$N))
	obs.growth <- c(diff(log(fit$N)),NA)

	par(mfrow=c(3,2))
	Nmax <- max(c(fit$K,fit$N),na.rm=T); if (Nmax>2*max(fit$N,na.rm=T)) Nmax <- 1.2*max(fit$N,na.rm=T)
	Nmin <- min(fit$N, na.rm=T)
	Nrange <- seq(trunc(Nmin*0.9),trunc(Nmax*1.1),by=0.1);  Nrange <- Nrange[Nrange>0.9]
	pred <- NA
	if (fit$theta!=0) pred <- fit$r1*(1-((Nrange^fit$theta-1)/(fit$K^fit$theta-1)))-0.5*fit$sd*exp(-log(Nrange))
	if (fit$theta==0) pred <- fit$r1-((fit$r1/log(fit$K-1))*log(Nrange))-0.5*fit$sd*exp(-log(Nrange))

	par(mfrow=c(4,2),mai=c(.7,.7,.3,.2))
	plot(Nrange[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))],pred[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))],type="n",ylab="predicted/observed r",xlab="N",
			ylim=range(c(obs.growth,pred[is.finite(pred)&(pred<max(3*abs(obs.growth),na.rm=T))]),na.rm=T),sub=paste("theta: ",round(fit$theta,2), "  r1: ", round(fit$r1,3), "  K: ", round(fit$K,1),
					"  gamma: ", round(fit$gamma,3), "  (sig_e)^2: ", round(fit$se,3)),cex.sub=0.7)
	lines(Nrange[is.finite(pred)],pred[is.finite(pred)],col="red",lwd=2)
	points(fit$N,obs.growth)
	sel <- which((abs(scale(fit$residuals))>1.5))
	if (length(sel)>0) text(fit$N[sel],obs.growth[sel],labels=year[sel],pos=4,cex=0.7)
	if (!is.null(title)) title(title)

	pred <- c(NA,(exp(fit$predicted+log(fit$N))[-length(fit$N)]))
	plot(year,fit$N,type="n",ylim=range(c(pred,fit$N),na.rm=T)-c(1,-1),ylab="observed/predicted N")
	if (fit$K<max(fit$N,na.rm=T))abline(h=fit$K,lty=3,col="red",lwd=1)
	points(year,fit$N)
	lines(year,fit$N)
	points(year,pred,col="red",pch=16)

	plot(year,obs.growth,type="n",ylim=range(c(obs.growth,fit$predicted),na.rm=T),xlim=range(year),ylab="observed/predicted r")
	points(year,obs.growth)
	lines(year,obs.growth)
	points(year,fit$predicted,col="red",pch=16)
	abline(h=0,lty=2)

	plot(year,fit$residuals,ylab="residuals")
	lines(year,fit$residuals)
	abline(h=0,lty=2)

	lmfit <- lm(fit$residuals~year,na.action="na.exclude")
	prob <- round(anova(lmfit)[1,5],3)
	plot(year,fit$residuals,ylab="residuals",sub=paste("p-val: ",round(prob,3)))
	if (prob<0.05) abline(lmfit,lwd=2,col="red") else abline(lmfit,lwd=1,lty=2)

	lmfit <- lm(fit$residuals~fit$N,na.action="na.exclude")
	prob <- round(anova(lmfit)[1,5],3)
	plot(fit$N,fit$residuals,ylab="residuals",xlab="N",sub=paste("p-val: ",round(prob,3)))
	if (prob<0.05) abline(lmfit,lwd=2,col="red") else abline(lmfit,lwd=1,lty=2)

	acf(fit$residuals,na.action=na.pass,lag.max=5)
	port <- portman.Q(fit$residuals,K=10)
	plot(c(1,2),c(1,2),type="n", axes=F,ylab="",xlab="")
	text(1.1,1.8,"Test for whiteness of residuals:",cex=1,pos=4)
	text(1.1,1.6,paste("p-val:", port[3],sep=" "),pos=4)
	if (port[3]>0.05) text(1.1,1.4,"Residuals are not significantly different from white noise",pos=4)
	if (port[3]<=0.05)text(1.1,1.4,"WARNING: Residuals are significantly different from white noise",pos=4)

}
################################################################################
################################################################################
################################################################################
################################################################################

#
# #
# #--------------------        simuler data    ----------------------------------#
# le <- 30
# r <- 0.5
# b <- -0.002
# sigD2 <- 0.3
# sigE2 <- 0.02
# x <- 1
#
# for (i in 2:le){
# 	x[i] <- x[i-1] + b*exp(x[i-1]) +
# 			r + sqrt(sigE2)*rnorm(1) + sqrt(sigD2/exp(x[i-1]))*rnorm(1) -
# 			0.5*sigE2 - sigD2*exp(-x[i-1])
# }
# ts.plot(x)
# #------------------------------------------------------------------------------#
#
# # estimer parametre
# fit <- thetalogistic.est(exp(x),1:length(x),sd=0.3,nboot=10, nbias=10)
# # se p? residualer m.m
# thetaplot(fit)
# # se p? stasjon?rfordelinga
# ts.plot(fit$f.stat.dist)
#
#
# # logistic model
# fit <- logismodfit(exp(x),1:length(x),sd=0.3)
# logistic.plot(fit)






#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#### logistic ###
m <- function(n, r, K) r*n-(r/K)*n^2
v <- function(n, sigE2, sigD2) sigE2*n^2 + sigD2*n


GreenFunc <- function(x,loww=0, X0, r, K, sigE2, sigD2){

	sv <- function(x) m(x, r, K)/v(x, sigE2, sigD2)
	mx <- function(x) 1/(v(x, sigE2=sigE2, sigD2=sigD2)*s(x))
	s <-  Vectorize( function(y) exp(-2*integrate(sv,lower=loww,upper=y,subdivisions=100)$value))
	S <-  function(x) integrate(s, lower = loww, upper = x)$value
	ifelse(x<=X0,2*mx(x)*S(x),2*mx(x)*S(X0))
}

ET <- function(lower=0,upper=Inf){
	integrate(GreenFunc,lower,upper,loww=lower)$value
}

QuasiStasj <- function(nr=200,from=0.1,to=(K*2),lower=0,upper=Inf, X0, r, K, sigE2, sigD2){
	vec <- seq(from,to,length=nr)
	res <- NA
	for (i in 1:length(vec)){ res[i] <- GreenFunc(vec[i],lower, X0, r, K, sigE2, sigD2) }
	res <- res/(sum(res))/diff(vec[1:2])
	names(res) <- vec
	res
}

simTS <- function(N0, nrYear, nrSeries, delta, r, K, sigE2, sigD2){
	nrSteps <- round(nrYear/delta)
	mat <- matrix(NA, nrSteps, nrSeries)
	mat[1,] <- N0
	U <- matrix(rnorm(ncol(mat)*nrow(mat)), nrow(mat), ncol(mat))
	for (i in 2:nrow(mat)){
		mat[i,] <- mat[i-1,] + m(mat[i-1,], r, K)*delta + sqrt(delta*v(mat[i-1,],sigE2, sigD2))*U[i-1,]
		mat[i,mat[i,]<0] <- NA
	}
	mat
}





stasj.logismod <- function(x){
	r <- x$s + 0.5*x$sig2
	K <- x$K
	sigE2 <- x$sig2
	sigD2 <- x$sigd
	X0 <- K # startverdi
	if (sigD2<0.01){
		nedre <- 1 # nedre grense (=1 hvis sigd2=0, = 0 hvis sigd2>0)
	}else{
		nedre <- 0
	}
	simN <- simTS(N0=K, nrYear=100, nrSeries=1000, delta=0.01, r=r, K=K, sigE2=sigE2, sigD2=sigD2)
	par(mfrow=c(2,1))
	hist(simN[-(1:10),],prob=T, xlim=c(0,max(na.omit(simN))*1.1), nclass=100, main="Kontinuerlig")
# quasi-stasj fordeling diff approx
	qs <- QuasiStasj(to=K*2, lower=nedre, X0=X0, r=r, K=K, sigE2=sigE2, sigD2=sigD2)
	lines(as.numeric(names(qs)), qs, col="red")
	simND <- simTS(N0=K, nrYear=1000, nrSeries=1000, delta=1 , r=r, K=K, sigE2=sigE2, sigD2=sigD2) # delta = 1!
	hist(simND[-(1:10),],prob=T, xlim=c(0,max(na.omit(simN))*1.1), nclass=30, main="Diskret",
			ylim=c(0, max(qs)))
	lines(as.numeric(names(qs)), qs, col="red")
	list("stat.dist"=qs, "var.stat.dist"=var(as.numeric(simN[-(1:10),])))
}


# simuler diffusjon

## simuler diskret (diskret != diff approx stasjon?rfordeling )
#simN <- simTS(N0=K, nrYear=1000, nrSeries=1000, delta=1) # delta = 1!
#hist(simN[-1,],prob=T, xlim=c(0,300), nclass=100, ylim=c(0, 0.025))
## quasi-stasj fordeling diff approx
#qs <- QuasiStasj(to=300, lower=1)
#lines(as.numeric(names(qs)), qs, col="red")

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################



# setwd("C:\\Berntes\\Mod\\DensDepDem\\Blue Tit")
# blue.tit.n<-read.delim("n.blue.tit.VL.txt",header=T)
#
# head (blue.tit.n)
# names (blue.tit.n)
# blue.tit.n$N<-blue.tit.n$Density
#
# blue.tit.n<-blue.tit.n[(blue.tit.n$Year>1963 & blue.tit.n$Year<2012),]
#
#
# sink ("thetalog.blue.tit.total.txt")
# fit.total <- thetalogistic.est(blue.tit.n$N,blue.tit.n$Year,sd=0,nboot=1000, nbias=1000)
# fit.total
# fit.total$var.stat.dist
# postscript("thetalog.blue.tit.total.eps")
# thetaplot(fit.total)
# dev.off()
# sink ()
#
# sink ("logistisk.blue.tit.total.txt")
# fit1.total<-logismodfit(blue.tit.n$N,blue.tit.n$Year,sd=0,nboot=1000)
# fit1.total
# stasj.logismod(fit1.total)
# postscript("logistic.blue.tit.total.eps")
# logistic.plot(fit1.total)
# dev.off()
# sink ()

