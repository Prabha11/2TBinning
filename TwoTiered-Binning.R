rm(list=ls())


library(Matrix)
library(mclust)
library(covRobust)

source('R_func/Two-Tiered_Binning_Func.R')

# Data set name (please ensure that this is the same as the name specified
# in the generate_features.sh script)
path <- 'sample_data/'
name <- 'simBG'


# Clustering parameters
p_scaling_GC 	<- 0.51154290
p_scaling_OFDEG	<- 0.07964383	# Note: this corresponds to a 4-mer OFDEG basis
p_OFDEG_basis 	<- 4
p_min_length 	<- 2000
p_max_length 	<- 1e99
p_max_Ns 		<- 0

# Tier 1
p_T1_noise_est 	<- 0.1
p_T1_noise_K 	<- 5
p_T1_epsilon_m 	<- 0.5
p_T1_epsilon_u 	<- 0.2

# Tier 2
p_T2_noise_est 	<- 0.4
p_T2_noise_K 	<- 20
p_T2_epsilon_m 	<- 0.2
p_T2_epsilon_u 	<- 0.2

nfile1 <- paste(path,name,'/',name,'_view1.OFDEG',sep="")
tfile <- paste(path,name,'/',name,'_view2.n4',sep="")
output_file <- paste(path,name,'/',name,"_output",sep="")

(cat("yes", fill=TRUE))
# Tier 1 - CLUSTERING
d1f <- read.table(file=nfile1, sep=",", header=T)
# Filtering outliers (you may need to adjust these values according to your data set; outliers will cause the mclust functions to fail)

plot(ofdeg~gc,d1f,xlab="GC content (%)",pch=20,col=densCols(cbind(d1f$gc,d1f$ofdeg)),ylab="OFDEG",cex.lab=0.9,cex.axis=0.8, sub="Initial datapoints")
abline(h=-0.095)
abline(h=-0.06)
d1 <- subset(d1f,d1f$ofdeg > -0.095 & d1f$ofdeg < -0.06)

plot(ofdeg~gc,d1,xlab="GC content (%)",pch=20,col=densCols(cbind(d1$gc,d1$ofdeg)),ylab="OFDEG",cex.lab=0.9,cex.axis=0.8, sub="OFDEG filtered datapoints")

d <- subset(d1,length >= p_min_length & length <= p_max_length)

nrow(d)
d <- subset(d, d$n <= p_max_Ns)
nrow(d)
sum(d$length)
d$ofdeg <- abs(d$ofdeg)
ds <- data.frame(id=d$id,gc=scale(d$gc,center=F,scale=p_scaling_GC ),ofdeg=scale(d$ofdeg,center=F,scale=p_scaling_OFDEG))
mt <- remap(ds,p_OFDEG_basis)
toclust <- cbind(mt$man,mt$dist)
ids <- mt$id

plot(toclust,xlab="U",pch=20,col=densCols(toclust),ylab="V",cex.lab=0.9,cex.axis=0.8)

clustL1 <- clusterCleanEM_L1(toclust,ids,pnoise=p_T1_noise_est,k=p_T1_noise_K)
clustL1.post <- clusterEM_PostProc(clustL1$mclust,gaussianCutoff=p_T1_epsilon_m,uncertaintyMax=p_T1_epsilon_u,clustL1$filtered.data,clustL1$filtered.ids)
vz <- clustL1.post$classification
lengths <- data.frame(id=d$id,length=d$length)
gcs <- data.frame(id=d$id,gc=d$gc)

# In most cases mclust will provide the correct number of clusters based on the BIC, but in some cases it is worthwhile to check that you are not overfitting
plot(clustL1$mclust,clustL1$filtered.data,what="BIC")
plot(toclust,xlab="U",pch=20,col="lightgray",ylab="V",cex.lab=0.9,cex.axis=0.8)
points(clustL1.post$data,pch=20,col=blues9[5])

# Tier 1 - SUMMARY
(perc.noise <- (nrow(d) - nrow(clustL1$filtered.data))/nrow(d))
(count.noise <- nrow(d) - nrow(clustL1$filtered.data))
m.L1 <- merge(data.frame(id=clustL1$filtered.ids,bin=clustL1$mclust$classification), lengths, by.x="id",by.y="id")
g.L1 <- merge(data.frame(id=clustL1$filtered.ids,bin=clustL1$mclust$classification), gcs, by.x="id",by.y="id")
(L1.count.bps <- tapply(m.L1$length,m.L1$bin,sum))
(L1.count <- tapply(m.L1$length,m.L1$bin,length))
(L1.mean.gc <- tapply(g.L1$gc,g.L1$bin,mean))
(uncertainty.L1 <- apply(clustL1$mclust$z,2,mean))

(perc.L1 <- nrow(vz)/nrow(d))
(count.L1 <- nrow(vz))
m.L1.post <- merge(clustL1.post$classification, lengths, by.x="id",by.y="id")
(L1.post.bps <- tapply(m.L1.post$length,m.L1.post$bin,sum))
(L1.post.count <- tapply(m.L1.post$length,m.L1.post$bin,length))
(filtered <- (L1.count + 19 - L1.count)/(L1.count))
(L1.pre <- data.frame(id=names(L1.count),count=L1.count,bps=L1.count.bps,gc=L1.mean.gc,uncertainty=uncertainty.L1))
(L1.post <- data.frame(id=names(L1.post.bps),count=L1.post.count,bps=L1.post.bps,filtered=round(filtered*L1.count,2),filtered.perc=round(100*filtered,2)))
(cat("plotted", fill=TRUE))
	
(nfile2 <- tfile)
(d2 <- read.table(file=nfile2, sep=",", header=T))
(mm <- merge(d2,vz,by.x="id",by.y="id"))
(v2 <- data.frame(id=mm$id,v1=mm$bin,v2=numeric(nrow(mm))))
((v1.classes <- unique(mm$bin)))
(features <- 4:(ncol(d2)-1))
	
(cln2 <- vector("list",length(v1.classes)))
(uncertainty.L2 <- vector("list",length(v1.classes)))
(noise.L2 <- vector("list",length(v1.classes)))
(noise.L2.count <- vector("list",length(v1.classes)))
(mclust.L2 <- vector("list",length(v1.classes)))

	
(for (i in 1:length(v1.classes)) {
	(v1.class <- v1.classes[i])
	(ii <- v1.class)
	(toclust.subset <- mm[mm$bin == v1.class,])
	(toclust <- toclust.subset[,features])
	(pca <- prcomp(toclust,scale=T))
	(perc <- round(pca$sdev/sum(pca$sdev)*100,2))
	
	# You may only need to the first two principal components, depending on your data set
	choices <- c(1,2,3)
	#choices <- c(1,2)
	
	(toclust.pca <- pca$x[,choices])
	(toclust <- toclust.pca)
	(ids <- toclust.subset[,1])
	(C_C.L2 <- clusterCleanEM_L1(toclust, ids, pnoise=p_T2_noise_est,k=p_T2_noise_K))
	(mclust.L2[[ii]] <- tapply(C_C.L2$filtered.ids,C_C.L2$mclust$classification,length))
	(noise.L2[[ii]] <- 1-nrow(C_C.L2$filtered.data)/nrow(toclust))
	(noise.L2.count[[ii]] <- nrow(toclust)-nrow(C_C.L2$filtered.data))
	if (C_C.L2$mclust$G > 1) {
		(uncertainty.L2[[ii]] <- apply(C_C.L2$mclust$z,2,mean))
	} else {
		(uncertainty.L2[[ii]] <- 0)
	}
	(cln2[[ii]] <- C_C.L2.post <- clusterEM_PostProc(C_C.L2$mclust,gaussianCutoff=p_T2_epsilon_m,uncertaintyMax=p_T2_epsilon_u,C_C.L2$filtered.data,C_C.L2$filtered.ids))
})
(clx.l2 <- cln2[[1]]$classification)
(det.all <- cbind(rep(1,length(cln2[[1]]$det)),cln2[[1]]$det))
(count.all <- cbind(rep(1,length(cln2[[1]]$count)),cln2[[1]]$count))
(clx.l2$bin <- as.numeric(paste(1,clx.l2$bin,sep="")))
(for (j in 2:length(cln2)) {
	(tmp.clx <- cln2[[j]]$classification)
	(tmp.clx$bin <- as.numeric(paste(j,tmp.clx$bin,sep="")))
	(clx.l2 <- rbind(clx.l2,tmp.clx))
	(det.all <- rbind(det.all,cbind(rep(j,length(cln2[[j]]$det)),cln2[[j]]$det)))
	(count.all <- rbind(count.all,cbind(rep(j,length(cln2[[j]]$count)),cln2[[j]]$count)))
})
(L.2 <- unlist(lapply(table(det.all[,1]),function(x) {return(seq(1,x))})))
(L.12 <- as.numeric(paste(det.all[,1],L.2,sep="")))
(count.each <- vector("list",length(cln2)))
(count.bps <- vector("list",length(cln2)))
(mean.gc <- vector("list",length(cln2)))
(for (cc in 1:length(cln2)) {
	count.each[[cc]] <- tapply(cln2[[cc]]$classification$id,cln2[[cc]]$classification$bin,length)
	m.L2 <- merge(cln2[[cc]]$classification, lengths, by.x="id",by.y="id")
	count.bps[[cc]] <- (tapply(m.L2$length,m.L2$bin,sum))
	g.L2 <- merge(cln2[[cc]]$classification, gcs, by.x="id",by.y="id")
	mean.gc[[cc]] <- (tapply(g.L2$gc,g.L2$bin,mean))
})
(m.L2.all <- merge(clx.l2, lengths, by.x="id",by.y="id"))
(tapply(m.L2.all$length,m.L2.all$bin,sum))
(perc.L2 <- nrow(clx.l2)/nrow(d))
(count.L2 <- nrow(clx.l2))


# Exclude bins? (i.e. due to an insufficient number of sequences)
(count.cxb <- table(clx.l2$bin))
(exclude <- names(count.cxb)[which(count.cxb < 30)])
(clx.l2.cond <- rep(TRUE,nrow(clx.l2)))
(for (i in 1:length(clx.l2.cond)) {
	if (clx.l2$bin[i] %in% exclude) {
		clx.l2.cond[i] <- FALSE
	}
})

# Tier 2 - SUMMARY
(filtered.L2 <- (unlist(mclust.L2)-unlist(count.each)))
(filtered.L2.perc <- round(100*filtered.L2/unlist(mclust.L2),2))
(uncertainty.L2.all <- round(unlist(uncertainty.L2),2))
(noise.L2.all <- round(100*unlist(noise.L2),2))
(noise.L2.count.all <- unlist(noise.L2.count))
(noise.L2.summary <- data.frame(id=v1.classes,noise.perc=noise.L2.all,noise=noise.L2.count.all))
(L2.bins <- data.frame(id=L.12,bps=unlist(count.bps),count=unlist(count.each),gc=round(100*unlist(mean.gc),2),filtered=filtered.L2,filtered.perc=filtered.L2.perc,uncertainty=uncertainty.L2.all))



# Output - summary
(write.table(L2.bins,file=paste(output_file,".L2_summary",sep=""),sep="\t",row.names=F,col.names=F,quote=F))
# Output - binning results (sequnece.id, bin)
(write.table(clx.l2[clx.l2.cond,],file=paste(output_file,'.L2_BINS',sep=""),sep="\t", row.names=F,col.names=F,quote=F))
(cat("yes", fill=TRUE))

# CHECK
(ans <- read.table(file="sample_data/simBG/sim.contig.ans",sep="\t",header=F))
(names(ans) <- c("id","taxon"))
(mx2 <- merge(clx.l2[clx.l2.cond,],ans,by.x="id",by.y="id"))
(taxon.r <- mx2$taxon[drop=TRUE])
(ac2 <- xtabs(~taxon.r+bin,mx2))
(ac2)

