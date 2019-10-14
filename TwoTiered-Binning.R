rm(list=ls())
# setwd("C:\\Users\\PRABHA\\Downloads\\2TBinning_v1.2.3.tar\\2TBinning_v1.2.3\\2TBinning")

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

# Tier 1 - CLUSTERING
d1f <- read.table(file=nfile1, sep=",", header=T)
# Filtering outliers (you may need to adjust these values according to your data set; outliers will cause the mclust functions to fail)

plot(ofdeg~gc,d1f,xlab="GC content (%)",pch=20,col=densCols(cbind(d1f$gc,d1f$ofdeg)),ylab="OFDEG",cex.lab=0.9,cex.axis=0.8)
abline(h=-0.095)
abline(h=-0.06)
d1 <- subset(d1f,d1f$ofdeg > -0.095 & d1f$ofdeg < -0.06) # filtering out some contigues as outliers.

plot(ofdeg~gc,d1,xlab="GC content (%)",pch=20,col=densCols(cbind(d1$gc,d1$ofdeg)),ylab="OFDEG",cex.lab=0.9,cex.axis=0.8)

d <- subset(d1,length >= p_min_length & length <= p_max_length) # filtering out again
write.table(d,file=paste(output_file,".filtered_contiges_1st_phase",sep=""),sep=",",row.names=F,col.names=F,quote=F)

print("first filter result number of rows")
nrow(d)
d <- subset(d, d$n <= p_max_Ns)
print("second filter result number of rows")
nrow(d)
print("total length of all sequences")
sum(d$length) # get the total length of all the sequnce
d$ofdeg <- abs(d$ofdeg) # getting absolute value of OFDEG
ds <- data.frame(id=d$id,gc=scale(d$gc,center=F,scale=p_scaling_GC ),ofdeg=scale(d$ofdeg,center=F,scale=p_scaling_OFDEG)) # scaled data
mt <- remap(ds,p_OFDEG_basis) # function is defined in R_func folder
#print (mt)
toclust <- cbind(mt$man,mt$dist)
ids <- mt$id
plot(toclust,xlab="U",pch=20,col=densCols(toclust),ylab="V",cex.lab=0.9,cex.axis=0.8)
clustL1 <- clusterCleanEM_L1(toclust,ids,pnoise=p_T1_noise_est,k=p_T1_noise_K) 
clustL1.post <- clusterEM_PostProc(clustL1$mclust,gaussianCutoff=p_T1_epsilon_m,uncertaintyMax=p_T1_epsilon_u,clustL1$filtered.data,clustL1$filtered.ids)
vz <- clustL1.post$classification
lengths <- data.frame(id=d$id,length=d$length) # create two dem array  of id, length
gcs <- data.frame(id=d$id,gc=d$gc) # create to dem array of id, gc

# In most cases mclust will provide the correct number of clusters based on the BIC, but in some cases it is worthwhile to check that you are not overfitting
plot(clustL1$mclust,clustL1$filtered.data,what="BIC")

plot(toclust,xlab="U",pch=20,col="lightgray",ylab="V",cex.lab=0.9,cex.axis=0.8) # all data send to clustering
points(clustL1.post$data,pch=20,col=blues9[5]) #  after post processing


# Tier 1 - SUMMARY
print("Tier 1 - SUMMARY")

(perc.noise <- (nrow(d) - nrow(clustL1$filtered.data))/nrow(d)) # calculate the percentage noise
(count.noise <- nrow(d) - nrow(clustL1$filtered.data))
m.L1 <- merge(data.frame(id=clustL1$filtered.ids,bin=clustL1$mclust$classification), lengths, by.x="id",by.y="id")
g.L1 <- merge(data.frame(id=clustL1$filtered.ids,bin=clustL1$mclust$classification), gcs, by.x="id",by.y="id")
(L1.count.bps <- tapply(m.L1$length,m.L1$bin,sum)) # tapply function use to break up a vector into groups defined by some classifying factor, compute a function on the subsets, and return the results in a convenient form
(L1.count <- tapply(m.L1$length,m.L1$bin,length))
(L1.mean.gc <- tapply(g.L1$gc,g.L1$bin,mean))
(uncertainty.L1 <- apply(clustL1$mclust$z,2,mean))


perc.L1 <- nrow(vz)/nrow(d)
count.L1 <- nrow(vz)
m.L1.post <- merge(clustL1.post$classification, lengths, by.x="id",by.y="id")
L1.post.bps <- tapply(m.L1.post$length,m.L1.post$bin,sum)
(L1.post.count <- tapply(m.L1.post$length,m.L1.post$bin,length))
(filtered <-(L1.count - L1.post.count)/(L1.count))
(L1.pre <- data.frame(id=names(L1.count),count=L1.count,bps=L1.count.bps,gc=L1.mean.gc,uncertainty=uncertainty.L1))
(L1.post <- data.frame(id=names(L1.post.bps),count=L1.post.count,bps=L1.post.bps,filtered=round(filtered*L1.count,2),filtered.perc=round(100*filtered,2)))

print("end tier one")

nfile2 <- tfile
d2 <- read.table(file=nfile2, sep=",", header=T)
mm <- merge(d2,vz,by.x="id",by.y="id") # results are contiges TTf's which are filtered after post process in tier 1
v2 <- data.frame(id=mm$id,v1=mm$bin,v2=numeric(nrow(mm))) #is the bin number
v1.classes <- unique(mm$bin) # bin categories
features <- 4:(ncol(d2)-1) #TNF features list

cln2 <- vector("list",length(v1.classes))
uncertainty.L2 <- vector("list",length(v1.classes))
noise.L2 <- vector("list",length(v1.classes))
noise.L2.count <- vector("list",length(v1.classes))
mclust.L2 <- vector("list",length(v1.classes))
	
for (i in 1:length(v1.classes)) { # loop through each bin
	v1.class <- v1.classes[i]
	ii <- v1.class
	toclust.subset <- mm[mm$bin == v1.class,] # selects contigs of perticular tier 1 bin
	toclust <- toclust.subset[,features] # selects the TNF set
	pca <- prcomp(toclust,scale=T) # return pCA with std deviations and rotaions
	perc <- round(pca$sdev/sum(pca$sdev)*100,2)
	
	# You may only need to the first two principal components, depending on your data set
	choices <- c(1,2,3)
	#choices <- c(1,2)

	toclust.pca <- pca$x[,choices] # select the first 3 priciple comp
	toclust <- toclust.pca
	ids <- toclust.subset[,1] # ids

	C_C.L2 <- clusterCleanEM_L1(toclust, ids, pnoise=p_T2_noise_est,k=p_T2_noise_K)
	mclust.L2[[ii]] <- tapply(C_C.L2$filtered.ids,C_C.L2$mclust$classification,length) # contain bins in tier 2 with number of contigs
	noise.L2[[ii]] <- 1-nrow(C_C.L2$filtered.data)/nrow(toclust)
	noise.L2.count[[ii]] <- nrow(toclust)-nrow(C_C.L2$filtered.data) # coise contains contig count

	if (C_C.L2$mclust$G > 1) {
		(uncertainty.L2[[ii]] <- apply(C_C.L2$mclust$z,2,mean))
	} else {
		(uncertainty.L2[[ii]] <- 0)
	}
	cln2[[ii]] <- C_C.L2.post <- clusterEM_PostProc(C_C.L2$mclust,gaussianCutoff=p_T2_epsilon_m,uncertaintyMax=p_T2_epsilon_u,C_C.L2$filtered.data,C_C.L2$filtered.ids)
}

clx.l2 <- cln2[[1]]$classification
det.all <- cbind(rep(1,length(cln2[[1]]$det)),cln2[[1]]$det)
count.all <- cbind(rep(1,length(cln2[[1]]$count)),cln2[[1]]$count)
clx.l2$bin <- as.numeric(paste(1,clx.l2$bin,sep="")) # bin number

for (j in 2:length(cln2)) {
	tmp.clx <- cln2[[j]]$classification
	tmp.clx$bin <- as.numeric(paste(j,tmp.clx$bin,sep=""))
	clx.l2 <- rbind(clx.l2,tmp.clx) # aggregate clusters 
	det.all <- rbind(det.all,cbind(rep(j,length(cln2[[j]]$det)),cln2[[j]]$det))
	count.all <- rbind(count.all,cbind(rep(j,length(cln2[[j]]$count)),cln2[[j]]$count)) # summery of bin assignmnets in each tier
}

L.2 <- unlist(lapply(table(det.all[,1]),function(x) {return(seq(1,x))}))
L.12 <- as.numeric(paste(det.all[,1],L.2,sep=""))
count.each <- vector("list",length(cln2))
count.bps <- vector("list",length(cln2))
mean.gc <- vector("list",length(cln2))

for (cc in 1:length(cln2)) {
	count.each[[cc]] <- tapply(cln2[[cc]]$classification$id,cln2[[cc]]$classification$bin,length) #bin assignmnet
	m.L2 <- merge(cln2[[cc]]$classification, lengths, by.x="id",by.y="id") # merge  id's
	count.bps[[cc]] <- (tapply(m.L2$length,m.L2$bin,sum)) # get the sum of each bin's len
	g.L2 <- merge(cln2[[cc]]$classification, gcs, by.x="id",by.y="id")
	mean.gc[[cc]] <- (tapply(g.L2$gc,g.L2$bin,mean)) # get GC content mean of each bins
}

m.L2.all <- merge(clx.l2, lengths, by.x="id",by.y="id")
tapply(m.L2.all$length,m.L2.all$bin,sum)
perc.L2 <- nrow(clx.l2)/nrow(d)
count.L2 <- nrow(clx.l2)


# Exclude bins? (i.e. due to an insufficient number of sequences)
count.cxb <- table(clx.l2$bin)
print(count.cxb)
exclude <- names(count.cxb)[which(count.cxb < 10)] # ginore bin counts less than this
clx.l2.cond <- rep(TRUE,nrow(clx.l2))
for (i in 1:length(clx.l2.cond)) {
	if (clx.l2$bin[i] %in% exclude) {
		clx.l2.cond[i] <- FALSE  # attribute cond - > False if ignored
	}
}

# Tier 2 - SUMMARY
filtered.L2 <- (unlist(mclust.L2)-unlist(count.each))
filtered.L2.perc <- round(100*filtered.L2/unlist(mclust.L2),2) # count of each bin after processing
uncertainty.L2.all <- round(unlist(uncertainty.L2),2)
noise.L2.all <- round(100*unlist(noise.L2),2)
noise.L2.count.all <- unlist(noise.L2.count)
noise.L2.summary <- data.frame(id=v1.classes,noise.perc=noise.L2.all,noise=noise.L2.count.all)
(L2.bins <- data.frame(bin_id=L.12,bps=unlist(count.bps),count=unlist(count.each),mean_gc=round(100*unlist(mean.gc),2),filtered=filtered.L2,filtered.perc=filtered.L2.perc,uncertainty=uncertainty.L2.all))

# Output - summary
write.table(L2.bins,file=paste(output_file,".L2_summary",sep=""),sep=",",row.names=F,col.names=F,quote=F)

# Output - binning results (sequnece.id, bin)
write.table(clx.l2[clx.l2.cond,],file=paste(output_file,'.L2_BINS',sep=""),sep="\t", row.names=F,col.names=F,quote=F) # write only clx.l2.cond true attributes
binnedContigs <- clx.l2[clx.l2.cond,]
#print(binnedContigs)

# FIND UNBINNED CONTIGS
unbinnedContigs <- subset(d1f, binnedContigs$id != d1f$id) # got unbinned sequences in 2T binning method
write.table(unbinnedContigs,file=paste(output_file,".unbinned_contiges",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
#print(unbinnedContigs)

# CHECK
ans <- read.table(file="sample_data/simBG/sim.contig.ans",sep="\t",header=F)
names(ans) <- c("id","taxon")
mx2 <- merge(clx.l2[clx.l2.cond,],ans,by.x="id",by.y="id")
taxon.r <- mx2$taxon[drop=TRUE]
(ac2 <- xtabs(~taxon.r+bin,mx2))

(cat("Accuracy check finished", fill=TRUE))

# Creating binned points csv
binned_points_with_features <- merge(mx2, d1f, by.x="id",by.y="id")
write.csv(binned_points_with_features, file = "binned_points.csv")

