
remap <- function(ds,k) {

	a <- -1
	b <- -1
	c <- -1

	if (k == 4) {
		# K = 4
		a <-  1.916047
		b <-  -3.753523
		c <- 2.691373
	} else if (k == 2) {
		# K = 2
		a <- 0.9542213
		b <- -1.878146
		c <- 1.850853
	}
	tp <- - b / (2*a)
	xx <- seq(min(ds$gc),max(ds$gc),len=500)
	yy <- c(c,a,b)%*% rbind(1,xx^2,xx)

	mapped <- data.frame(id=ds$id,xk=numeric(nrow(ds)),yk=numeric(nrow(ds)),x=numeric(nrow(ds)),y=numeric(nrow(ds)),dist=numeric(nrow(ds)),man=numeric(nrow(ds)))
	# print(ds)
	xtt <- apply(ds[,2:3],1,point2poly,a,b,c) # like a loop for each row
	# print(xtt)
	xtm <- matrix(xtt,ncol=2,byrow=T)
	# print(xtm)
	mapped$xk <- ds[,2]  # GC value
	mapped$yk <- ds[,3]	# OFDEG value
	mapped$x <- xtm[,1]	# GC values transformed to new PC 
	mapped$y <- xtm[,2] # OFDEG values transformed to new PC 
	
	signs <- rep(1,nrow(mapped))
	cond <- (mapped[,2] < tp & mapped[,2] < mapped[,4]) | (mapped[,2] > tp & mapped[,2] > mapped[,4]) | (mapped[,3] < mapped[,5])
	signs[cond] <- -1
	xts <- signs * sqrt( (mapped[,2]-mapped[,4])^2 + (mapped[,3]-mapped[,5])^2 ) # calculate the distance eucledian
	mapped$dist <- xts
	s <- sort(mapped$x,index.return=T) # sort by x
	mapped.s <- mapped[s$ix,]  # save the acending order 2d array in to mapped.s
	# print(mapped)
	sxx <- diff(mapped.s$x)^2
	sxy <- diff(mapped.s$y)^2
	dist <- sqrt(sxx + sxy) # this distance is mahalanobis dist i think
	ddx <- cumsum(c(0,dist))
	# print(ddx)
	mapped.s$man <- ddx
	
	mapped.s$xk <- as.numeric(mapped.s$xk)
	mapped.s$yk <- as.numeric(mapped.s$yk)
	return(mapped.s)
}
	
point2poly <- function(xkyk,a,b,c) {
	library(polynom)
	xk <- xkyk[1]   #contain the gc content values
	yk <- xkyk[2]	#contain the ODFEG values
	u <- 2*a^2
	v <- 3*a*b
	w <- b^2 + 2*a*c + 1 - 2*a*yk
	z <- b*c - xk - yk*b
	p <- polynomial(c(z,w,v,u))
	sp <- solve(p) #solve the equation
	sp.im <- unlist(lapply(sp,Im))
	sp.r <- as.double(sp[sp.im==0])
	
	d <- numeric(length(sp.r))
	xt <- numeric(length(sp.r))
	yt <- numeric(length(sp.r))
	for (i in 1:length(d)) {
		xt[i] <- sp.r[i]
		yt[i] <- a*xt[i]^2 + b * xt[i] + c
		d[i] <- sqrt( (xk - xt[i])^2 + (yk - yt[i])^2 ) # eucleadian distances
	}
	min.dist <- which.min(d) # return the index of the first layer that has the min 
	xm <- xt[min.dist]
	ym <- yt[min.dist]
	return(c(xm,ym)) # return the values of new cordinates
}

clusterCleanEM_L1 <- function(toclust, ids, pnoise,k,Gs=1:10) {
	library(mclust)
	library(covRobust)
	if (pnoise > 0) {
		xcve <- cov.nnve(toclust,pnoise=pnoise,k=k) #noise removal using NNVE method
		noise <- xcve$classification == 0
		toclustN1 <- toclust[!noise,] # remove noisy data 
		idsN1 <- ids[!noise] 
	} else {
		toclustN1 <- toclust
		idsN1 <- ids
	}
	m.o <- Mclust(toclustN1,modelNames="VVV",G=Gs) # g - > specify number of clusters 
	# modelNames - > indicating the models to be fitted in the EM phase of clustering
	ret <- vector("list",3)
	ret[[1]] <- m.o
	ret[[2]] <- toclustN1
	ret[[3]] <- idsN1
	names(ret) <- c("mclust","filtered.data","filtered.ids")
	return(ret)
}

clusterEM_PostProc <- function(m.o,gaussianCutoff,uncertaintyMax,toclust,ids) {

	choose.k <- length(unique(m.o$classification)) #remove duplicates
	if (0 %in% unique(m.o$classification)) {
		choose.k <- choose.k - 1
	}
	m <- vector("list", choose.k)
	s.inv <- vector("list", choose.k)
	ps <- vector("list", choose.k)
	p.gauss <- function(x, m, s.inv) { return(exp( -0.5*t(x - m)%*%s.inv%*%(x - m)) ) }
	de <- c(0)
	for (i in 1:choose.k) {
		m[[i]] <- m.o$parameters$mean[,i]
		s <- m.o$parameters$variance$sigma[,,i]
		s.inv[[i]] <- solve(s)
		ps[[i]] <- apply(toclust,1,p.gauss, m[[i]], s.inv[[i]])
		de[i] <- det(m.o$parameters$variance$sigma[,,i])
	}
	#x <- grid1( 100, range = range(toclust[,1]))
	#y <- grid1( 100, range = range(toclust[,2]))
	#xy <- grid2(x,y)
	if (0 %in% unique(m.o$classification)) {
		m.o$parameters$pro <- m.o$parameters$pro[-length(m.o$parameters$pro)]
	}
	#xyDens <- dens(modelName = m.o$modelName, data = xy, parameters = m.o$parameters)
	#xyDens <- matrix(xyDens, nrow = length(x), ncol = length(y))
	
	#m.dens <- c(0)
	#for (iii in 1:choose.k) {
	#	m.dens[iii] <- xyDens[which(x >= m[[iii]][1])[1],which(y >= m[[iii]][2])[1]]
	#}
	
	count <- table(m.o$classification)

	b.de <- boxplot(de,plot=F)$stats
	outl <- b.de[5]
	k.rem <- c(0)

	cutoff <- rep(gaussianCutoff,choose.k)

	tmp <- matrix(ncol=length(ps), nrow=length(ps[[1]]))
	for (i in 1:length(ps)) {
		tmp[,i] <- ps[[i]] > cutoff[i]
	}
	check.cond <- as.logical(apply(tmp, 1, sum))
	condition <- check.cond & !(m.o$classification %in% k.rem)
	if (m.o$G > 1) {
		uncertainty <- 1-apply(m.o$z,1,max)
		condition <- condition & (uncertainty < uncertaintyMax)
	}
	g.data <- toclust[condition, ]
	cln <- data.frame(id=ids[condition], bin=m.o$classification[condition])

	ret <- vector("list",5)
	ret[[1]] <- cln
	ret[[2]] <- de
	ret[[3]] <- count
	ret[[4]] <- c(0) # DENSITY REMOVED
	ret[[5]] <- g.data
	names(ret) <- c("classification","det","count","density","data")
	return(ret)

}


clusterCleanEM_L2 <- function(toclust, ids, pnoise) {
	library(mclust)
	library(covRobust)

	if (pnoise > 0) {
		xcve1 <- cov.nnve(toclust,pnoise=pnoise,k=12)
		noise1 <- xcve1$classification == 0
		toclustN1 <- toclust[!noise1,]
		idsN1 <- ids[noise1]
		idsG1 <- ids[!noise1]
	
		xcve2 <- cov.nnve(toclustN1,pnoise=pnoise,k=8)
		noise2 <- xcve2$classification == 0
		toclustN2 <- toclustN1[!noise2,]
		idsN2 <- idsG1[noise2]
		idsG2 <- idsG1[!noise2]
	
		idsN.all <- unique(as.factor(c(as.character(idsN1),as.character(idsN2))))
	
		noise <- ids %in% idsN.all
	} else {
		toclustN2 <- toclust
		idsG2 <- ids
	}

	m.o <- Mclust(toclustN2,modelNames="VVV")
	
	ret <- vector("list",3)
	ret[[1]] <- m.o
	ret[[2]] <- toclustN2
	ret[[3]] <- idsG2
	names(ret) <- c("mclust","filtered.data","filtered.ids")
	return(ret)
}

clusterNoiseEM_L2 <- function(toclust, ids, pnoise) {

	xcve1 <- cov.nnve(toclust,pnoise=pnoise,k=12)
	noise1 <- xcve1$classification == 0
	toclustN1 <- toclust[!noise1,]
	idsN1 <- ids[noise1]
	idsG1 <- ids[!noise1]

	xcve2 <- cov.nnve(toclustN1,pnoise=pnoise,k=8)
	noise2 <- xcve2$classification == 0
	toclustN2 <- toclustN1[!noise2,]
	idsN2 <- idsG1[noise2]
	idsG2 <- idsG1[!noise2]
	
	idsN.all <- unique(as.factor(c(as.character(idsN1),as.character(idsN2))))
	noise <- ids %in% idsN.all

	toclustNbic <- mclustBIC(toclust,modelNames="VVV",initialization = list(noise = noise))
	toclustNsummary <- summary(toclustNbic, toclust)
	gmmN <- toclustNsummary
		
	m.o <- gmmN
	
	ret <- vector("list",3)
	ret[[1]] <- m.o
	ret[[2]] <- toclust
	ret[[3]] <- ids
	names(ret) <- c("mclust","data","ids")
	return(ret)
}

