if("pheatmap" %in% rownames(installed.packages()) == FALSE) {install.packages("pheatmap")}; require(pheatmap)  
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}; require(scales)
if("ape" %in% rownames(installed.packages()) == FALSE) {install.packages("ape")}; require(ape)
if("pegas" %in% rownames(installed.packages()) == FALSE) {install.packages("pegas")}; require(pegas)
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}; require(RColorBrewer)
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}; require(dplyr)
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}; require(ggplot2)

setwd("/directory/")

id <- "scaffold3"
start <- 1
end <- 9443038

########### input data ################
fasta_file_M <- paste("/directory/file/fasta")
stats_file_M <- paste("/directory/file/saguaro.stats.txt")
fasta_file_V <- paste("/directory/file/fasta")
stats_file_V <- paste("/directory/file/saguaro.stats.txt")

chromosome <- paste("Chromosome", id, sep=" ")

stats_M <- read.table(stats_file_M, header=FALSE)
fasta_M <- ape::read.dna(fasta_file_M, format='fasta')
stats_V <- read.table(stats_file_V, header=FALSE)
fasta_V <- ape::read.dna(fasta_file_V, format='fasta')

header <- c("Chromosome", "Start", "End", "cactus", "pi1", "pi2", "dxy", "HBK", "HSM")

names(stats_M) <- header
names(stats_V) <- header

# correct for zero based start
stats_V$Start <- stats_V$Start + start;stats_V$End <- stats_V$End + start
end <- 9443038
stats_M$Start <- stats_M$Start + start;stats_M$End <- stats_M$End + start

# for plotting remove coordinates of zero coverage region in P. nyererei reference genome
stats_M_nocov <- stats_M[with(stats_M,  !((Start %in% 404391:431327))), ];stats_M_nocov <- stats_M_nocov[with(stats_M_nocov,  !((End %in% 404391:431327))), ]
stats_V_nocov <- stats_V[with(stats_V,  !((Start %in% 404391:431327))), ];stats_V_nocov <- stats_V_nocov[with(stats_V_nocov,  !((End %in% 404391:431327))), ]

stats_M_nocov$new_start <- c();
stats_M_nocov$new_end <- c()
stats_V_nocov$new_start <- c();
stats_V_nocov$new_end <- c()

stats_M_nocov$new_start[stats_M_nocov$Start > 431327] <- stats_M_nocov$Start[stats_M_nocov$Start > 431327]-26936;stats_M_nocov$new_end[stats_M_nocov$End > 431327] <- stats_M_nocov$End[stats_M_nocov$End > 431327]-26936
stats_M_nocov$new_start[stats_M_nocov$Start < 431327] <- stats_M_nocov$Start[stats_M_nocov$Start < 431327];stats_M_nocov$new_end[stats_M_nocov$End < 431327] <- stats_M_nocov$End[stats_M_nocov$End < 431327]
stats_V_nocov$new_start[stats_V_nocov$Start > 431327] <- stats_V_nocov$Start[stats_V_nocov$Start > 431327]-26936;stats_V_nocov$new_end[stats_V_nocov$End > 431327] <- stats_V_nocov$End[stats_V_nocov$End > 431327]-26936
stats_V_nocov$new_start[stats_V_nocov$Start < 431327] <- stats_V_nocov$Start[stats_V_nocov$Start < 431327];stats_V_nocov$new_end[stats_V_nocov$End < 431327] <- stats_V_nocov$End[stats_V_nocov$End < 431327]

# set negative FST values to zero
stats_M_nocov$HSM[stats_M_nocov$HSM < 0] <- 0
stats_V_nocov$HSM[stats_V_nocov$HSM < 0] <- 0

########## read data from TWISST results ##################
weights_file_V <- "/directory/weights.tsv"
window_data_file_V <- "/directory/data.tsv"
topos_phyl_V = read.tree("/directory/topologies.trees")
topos_V = readLines("/directory/topologies.trees")
topos_phyl_names_V = read.tree("/directory/topologies.with.species.names.trees")

weights_V = read.table(weights_file_V, header = T)
weights_V <- weights_V / apply(weights_V, 1, sum) #normalise rows so weights sum to 1
topoNames_V = names(weights_V) #retrieve the names of the topologies
window_data_V = read.table(window_data_file_V, header = T)
good_rows = which(is.na(apply(weights_V,1,sum)) == F) #exclude any rows where data is missing
weights_V <- weights_V[good_rows,]
window_data = window_data_V[good_rows,]
data.df.V <- cbind(weights_V, window_data_V)

# for plotting remove coordinates of zero coverage region in P. nyererei reference genome
data.df.V.nocov <- data.df.V[with(data.df.V,  !((start %in% 404391:431327))), ];data.df.V.nocov <- data.df.V.nocov[with(data.df.V.nocov,  !((end %in% 404391:431327))), ]
data.df.V.nocov$new_start <- c();data.df.V.nocov$new_end <- c()
data.df.V.nocov$new_start[data.df.V.nocov$start > 431327] <- data.df.V.nocov$start[data.df.V.nocov$start > 431327]-26936;data.df.V.nocov$new_end[data.df.V.nocov$end > 431327] <- data.df.V.nocov$end[data.df.V.nocov$end > 431327]-26936
data.df.V.nocov$new_start[data.df.V.nocov$start < 431327] <- data.df.V.nocov$start[data.df.V.nocov$start < 431327];data.df.V.nocov$new_end[data.df.V.nocov$end < 431327] <- data.df.V.nocov$end[data.df.V.nocov$end < 431327]
# cut end
values_to_keep <- data.df.V.nocov$new_start[data.df.V.nocov$new_start <= 512000]
length(values_to_keep)
rows_to_keep<- data.df.V.nocov[1:length(values_to_keep),]
data.df.V.nocov <- rows_to_keep

data.df.V.nocov <- data.df.V.nocov[, -c(107:111)]
data.df.V.str <- data.df.V.nocov %>% select(scaffold, new_start, new_end, topo1, topo2, topo3, topo16, topo17, topo18, topo76, topo77, topo78)
data.df.V.nst <- data.df.V.nocov[, -c(76:78)]
data.df.V.nst <- data.df.V.nst[, -c(16:18)]
data.df.V.nst <- data.df.V.nst[, -c(1:3)]

data.df.V.str$weight_sum <- NA
data.df.V.str$weight_sum <- rowSums(data.df.V.str[4:12])

weights_file_M <- "/directory/weights.tsv"
#coordinates file for each window
window_data_file_M <- "/directory/data.tsv"
#files with topologies
topos_phyl_M = read.tree("/directory/topologies.trees")
topos_M = readLines("/directory/topologies.trees")
topos_phyl_names_M = read.tree("/directory/topologies.with.species.names.trees")
########## read data ##################
weights_M = read.table(weights_file_M, header = T)
#normalise rows so weights sum to 1
weights_M <- weights_M / apply(weights_M, 1, sum)
#retrieve the names of the topologies
topoNames_M = names(weights_M)
window_data_M = read.table(window_data_file_M, header = T)
#exclude any rows where data is missing
good_rows = which(is.na(apply(weights_M,1,sum)) == F)
weights <- weights_M[good_rows,]
window_data = window_data_M[good_rows,]
data.df.M <- cbind(weights_M, window_data_M)

values_to_keep <- stats_M_nocov$new_start[stats_M_nocov$new_start <= 512000]
length(values_to_keep)
rows_to_keep<- stats_M_nocov[1:length(values_to_keep),]
stats_M_nocov <- rows_to_keep

data.df.M.nocov <- data.df.M[with(data.df.M,  !((start %in% 404391:431327))), ];data.df.M.nocov <- data.df.M.nocov[with(data.df.M.nocov,  !((end %in% 404391:431327))), ]
data.df.M.nocov$new_start <- c();data.df.M.nocov$new_end <- c()
data.df.M.nocov$new_start[data.df.M.nocov$start > 431327] <- data.df.M.nocov$start[data.df.M.nocov$start > 431327]-26936;data.df.M.nocov$new_end[data.df.M.nocov$end > 431327] <- data.df.M.nocov$end[data.df.M.nocov$end > 431327]-26936
data.df.M.nocov$new_start[data.df.M.nocov$start < 431327] <- data.df.M.nocov$start[data.df.M.nocov$start < 431327];data.df.M.nocov$new_end[data.df.M.nocov$end < 431327] <- data.df.M.nocov$end[data.df.M.nocov$end < 431327]

data.df.M.nocov <- data.df.M.nocov[, -c(107:111)]
data.df.M.str <- data.df.M.nocov %>% select(scaffold, new_start, new_end, topo1, topo2, topo3, topo16, topo17, topo18, topo76, topo77, topo78)
data.df.M.nst <- data.df.M.nocov[, -c(76:78)]
data.df.M.nst <- data.df.M.nst[, -c(16:18)]
data.df.M.nst <- data.df.M.nst[, -c(1:3)]

data.df.M.str$weight_sum <- NA
data.df.M.str$weight_sum <- rowSums(data.df.M.str[4:12])

# apply smoothening function for FST regions
HSM_V_01 <- loess(HSM ~ new_start, data=stats_V_nocov, span=0.004)
V_smoothed01 <- predict(HSM_V_01)

HSM_M_01 <- loess(HSM ~ new_start, data=stats_M_nocov, span=0.03)
M_smoothed01 <- predict(HSM_M_01)

pdf(file = paste0("Plot_TWISST_saguaro_FST_smoothened_both_populations.pdf"), width = 12, height = 6)
par(mar=c(5,4,4,5)+.1)
plot(NA,NA,xlim=c(312000, 512000),ylim=c(-1,1),axes=F,bty="n",xlab="",ylab="")
axis(side=1, at=c(312000, 412000, 512000), label=c(312,412,512))
axis(side=2, at=c(0.5, 1), labels = T, col.axis="darkorange")
axis(side=2, at=c(0, 0.5, 1), labels = F, col="darkorange")
axis(side=2, at=c(0, -0.5, -1), labels = c(0, 0.5, 1), col.axis="skyblue")
axis(side=2, at=c(0, -0.5, -1), labels = F, col="skyblue")
axis(4, at=c(0, 0.5, 1))
axis(4, at=c(0, -0.5, -1), labels = F)
axis(4, at=c(-0.5, -1), labels = c(0.5, 1))
mtext(side=1, text="Position [bp] on scaffold 3", cex=0, line=2, outer=FALSE)
mtext(side=2, text="Weight", cex=1, line=2.3, at=0, adj = 0)
mtext(side=4, text="FST", cex=1, line=2.3, at=0, adj = 0.5)
par(xpd=T)
# annotation of genes on scaffold
Agrp2_Exon_01 <- rect(438849-26936, 1, 438957-26936, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
Agrp2_Exon_02 <- rect(443411-26936, 1, 443485-26936, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
Agrp2_Exon_03 <- rect(448241-26936, 1, 448408-26936, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="agrp2", cex=0.5, line=0.4, outer=FALSE, at=443411-26936, col = "black")
line <- rect(438849-26936, 1.024, 448241-26936, 1.025, density = NULL, angle = 45, col = "black", border = "grey40")
atp6V0d2 <- rect(452121-26936, 1, 455894-26936, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="atp6V0d2", cex=0.5, line=0.4, outer=FALSE, at=427071, col = "black")
eif4a2	<- rect(313963, 1, 	318639, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="eif4a2", cex=0.5, line=0.4, outer=FALSE, at=(313963+318639)/2, col = "black")
setd9	<- rect(321003, 1, 	321854, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="setd9", cex=0.5, line=0.4, outer=FALSE, at=(321003+321854)/2, col = "black")
dcun1d1	<- rect(322145, 1, 	327779, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="dcun1d1", cex=0.5, line=0.4, outer=FALSE, at=(322145+327779)/2, col = "black")
ENSPNYG00000001687	<- rect(328305, 1, 	336628, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="ENSPNYG00000001687", cex=0.5, line=0.4, outer=FALSE, at=(328305+336628)/2, col = "black")
ENSPNYG00000001655	<- rect(344147, 1, 	346030, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="ENSPNYG00000001655", cex=0.5, line=0.4, outer=FALSE, at=(344147+346030)/2, col = "black")
ripk2	<- rect(362537, 1, 	373966, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="ripk2", cex=0.5, line=0.4, outer=FALSE, at=(362537+373966)/2, col = "black")
zgc	<- rect(475374, 1, 	479416, 1.05, density = NULL, angle = 45, col = "grey40", border = "grey40")
mtext(side=3, text="zgc:153115", cex=0.5, line=0.4, outer=FALSE, at=(475374+479416)/2, col = "black")

sep <- rect(312000, 1, 512000, 1, density = NULL, angle = 45, col = "black", border = "black")
par(xpd=F)
for (row in 1:nrow(data.df.V.str)){
  rect(data.df.V.str$new_start[row], 0, data.df.V.str$new_end[row], data.df.V.str$weight_sum[row],col="orange",border=NA)
}
par(new=TRUE)
for (row in 1:nrow(data.df.M.str)){
  rect(data.df.M.str$new_start[row], 0, data.df.M.str$new_end[row], -data.df.M.str$weight_sum[row],col="skyblue",border=NA)
}
par(new=TRUE)
plot(((stats_V_nocov$new_start+stats_V_nocov$new_end)/2), stats_V_nocov$HSM, type="p", col=alpha("black",alpha=0.5), pch = 16, xlab=NA, ylab=NA, xlim=c(312000, 512000), ylim=c(-1, 1), axes=F)
par(new=TRUE)
plot(stats_V_nocov$new_start, stats_V_nocov$HSM, type="n", xlab=NA, ylab=NA, xlim=c(312000, 512000), ylim=c(-1, 1), axes=F)
lines(V_smoothed01, x=stats_V_nocov$new_start, col="black", lwd=2)

par(new=TRUE)
plot(((stats_M_nocov$new_start+stats_M_nocov$new_end)/2), -(stats_M_nocov$HSM), type="p",  col=alpha("black",alpha=0.5), pch = 16, xlab=NA, ylab=NA, xlim=c(312000, 512000), ylim=c(-1, 1), axes=F)
par(new=TRUE)
plot(stats_M_nocov$new_start, -(stats_M_nocov$HSM), type="n", xlab=NA, ylab=NA, xlim=c(312000, 512000), ylim=c(-1, 1), axes=F)
lines(-(M_smoothed01), x=stats_M_nocov$new_start, col="black", lwd=2)
abline(h=0, col = "black") # seperator
dev.off()

########################################################################################################################################################
### Fig. 2B candidate regions (FST>0.75) position
########################################################################################################################################################

setwd("/directory/")

id <- "scaffold3"
start <- 1
end <- 9443038
min<-438171
max<-448838
stats_file_M <- paste("/Volumes/ORYX/Backup/Baits/fasta/Pnye/new/stats_LM_new_haplotypes_scaffold3_1_9443038.txt")
stats_file_V <- paste("/Volumes/ORYX/Backup/Baits/fasta/Pnye/new/stats_LV_new_haplotypes_scaffold3_1_9443038.txt")

stats_M <- read.table(stats_file_M, header=FALSE)
stats_V <- read.table(stats_file_V, header=FALSE)
header <- c("Chromosome", "Start", "End", "cactus", "pi1", "pi2", "dxy", "HBK", "HSM")
names(stats_M) <- header
names(stats_V) <- header
stats_V$Start <- stats_V$Start + start;stats_V$End <- stats_V$End + start
stats_M$Start <- stats_M$Start + start;stats_M$End <- stats_M$End + start
# set negative FST to zero
stats_V$HSM[stats_V$HSM < 0] <- 0
stats_M$HSM[stats_M$HSM < 0] <- 0

stats_V <- data.frame(append(stats_V, c(Lake='LV'), after=11))
stats_M <- data.frame(append(stats_M, c(Lake='LM'), after=11))
stats_all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(stats_V, stats_M))
stats_all<-stats_all[rev(order(stats_all$HSM)),] # Order by name
stats_all$ID <- c(1:nrow(stats_all))

stats_all<-stats_all[stats_all$End>min & stats_all$Start<max,]
stats_all$Start[stats_all$Start<min]<-min
stats_all$End[stats_all$End>max]<-max

Start <- cbind(stats_all$Start,stats_all$HSM)
End <- cbind(stats_all$End,stats_all$HSM)
#ID <- cbind(stats_all$ID,stats_all$HSM)
Lake <- cbind(as.character(stats_all$Lake),stats_all$ID)

stats_start <- cbind(Start, Lake) 
stats_end <- cbind(End, Lake)

###############  ##########################################################################################
stats <- as.data.frame(rbind(stats_start, stats_end))
stats <- as.data.frame(stats)
header <- c("Position", "FST", "Lake","ID")
names(stats) <- header
stats$FST <-as.numeric(as.character(stats[,2]))
stats$Position <-as.numeric(as.character(stats[,1]))
stats$ID <-as.numeric(as.character(stats[,4]))

stats<-stats[order(stats$ID,stats$Position),] # Order by name

pdf(file = paste0("LVLM_weights_FST_438849_448241_v2.pdf"), width = 12, height = 6)
par(mar=c(5,4,4,5)+.1)
plot(NA,NA,xlim=c(438849, 448241),ylim=c(-1,1),axes=F,bty="n",xlab="",ylab="")
axis(side=1, at=c(438849, 448241))
axis(side=2, at=c(0, 0.5, 1), labels = T)
axis(side=2, at=c(0, -0.5, -1), labels = c(0.,0.5, 1))
mtext(side=1, text="Position on scaffold 3", cex=0, line=2, outer=FALSE)
mtext(side=2, text="FST", cex=1, line=2.3, at=0.25, adj = 0)
mtext(side=2, text="Weight", cex=1, line=2.3, at=-0.75, adj = 0)
Agrp2_Exon_01 <- rect(438849, 1, 438957, 1.05, density = NULL, angle = 45, col = "black", border = "black")
Agrp2_Exon_02 <- rect(443411, 1, 443485, 1.05, density = NULL, angle = 45, col = "black", border = "black")
Agrp2_Exon_03 <- rect(448409, 1, 449024, 1.05, density = NULL, angle = 45, col = "black", border = "black")
UTR_iso1 <- rect(438775, 1, 438848, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
UTR_iso2.1 <- rect(438835, 1, 438848, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
UTR_iso2.2 <- rect(438668, 1, 438711, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
#EnhA <- rect(442320, 1, 443410, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")

mtext(side=3, text="agrp2", cex=1, line=0.75, outer=FALSE, at=443411, col = "black")
line <- rect(438800, 1.024, 448241, 1.025, density = NULL, angle = 45, col = "black", border = NA)
par(new=TRUE)
for (row in 1:nrow(stats)){
  rect(stats_all$Start[stats_all$Lake=="LM"][row], 0, stats_all$End[stats_all$Lake=="LM"][row], stats_all$HSM[stats_all$Lake=="LM"][row],col=alpha("skyblue",alpha=0.5),border=NA)
}
par(new=TRUE)
for (row in 1:nrow(stats)){
  rect(stats_all$Start[stats_all$Lake=="LV"][row], 0, stats_all$End[stats_all$Lake=="LV"][row], stats_all$HSM[stats_all$Lake=="LV"][row],col=alpha("darkorange",alpha=0.5),border=NA)
}

line <- rect(438800, 0, 448250, 0.01, density = NULL, angle = 45, col = "black", border = NA)
par(new=TRUE)
for (row in 1:nrow(data.df.M.str)){
  rect((data.df.M.str$new_start[row]+26936), 0, (data.df.M.str$new_end+26936)[row], -data.df.M.str$weight_sum[row],col=alpha("skyblue",alpha=0.5),border=NA)
}
par(new=TRUE)
for (row in 1:nrow(data.df.V.str)){
  rect((data.df.V.str$new_start+26936)[row], 0, (data.df.V.str$new_end+26936)[row], -data.df.V.str$weight_sum[row],col=alpha("darkorange",alpha=0.5),border=NA)
}
dev.off()

########################################################################################################################################################
########################################################################################################################################################
setwd("/Volumes/ORYX/Backup/Baits/fasta/Pnye/new/2020_resubmission")
header <- c("Chromosome", "Start", "End", "cactus", "pi1", "pi2", "dxy", "HBK", "HSM")
min<-438171
max<-448838

LM<-read.csv("/Volumes/ORYX/Backup/Baits/fasta/Pnye/new/2020_resubmission/stats_20200718_LM2_haplotypes_scaffold3_1_9443038.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(LM)<-header
LM$HSM[LM$HSM<0]<-0
LM<-LM[LM$End>min & LM$Start<max,]
LM$Start[LM$Start<min]<-min
LM$End[LM$End>max]<-max
LV<-read.csv("/Volumes/ORYX/Backup/Baits/fasta/Pnye/new/stats_LV_new_haplotypes_scaffold3_1_9443038.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(LV)<-header
LV$HSM[LV$HSM<0]<-0
LV<-LV[LV$End>min & LV$Start<max,]
LV$Start[LV$Start<min]<-min
LV$End[LV$End>max]<-max

d1=data.frame(x1=LM$Start,x2=LM$End,y1=rep(1.5,length(LM$End)),y2=rep(2.5,length(LM$End)))
d2=data.frame(x1=LV$Start,x2=LV$End,y1=rep(3,length(LV$End)),y2=rep(4,length(LV$End)))
d3=data.frame(x1=c(438849,443411,448241),x2=c(438957,443485,448408),y1=rep(0,3),y2=rep(1,3))

LV$col <- NA
LM$col <- NA

LM$HSM <- round(LM$HSM, digits = 2)
LV$HSM <- round(LV$HSM, digits = 2)

for (row in 1:nrow(LM)){
    if (LM$HSM[row] > 0.75)
      {
      LM$col[row] <- as.character("#87CEEB")
      }
    else if (LM$HSM[row] > 0.5){
    LM$col[row] <- as.character((alpha("#87CEEB",alpha=0.75)))
    }
    else if (LM$HSM[row] > 0.25){
    LM$col[row] <- as.character((alpha("#87CEEB",alpha=0.50)))
    }
    else if (LM$HSM[row] > 0.125){
    LM$col[row] <- as.character((alpha("#87CEEB",alpha=0.25)))
    }
    if (LM$HSM[row] < 0.125){
    LM$col[row] <- as.character((alpha("#87CEEB",alpha=0.01)))
  }
}

for (row in 1:nrow(LV)){
  if (LV$HSM[row] > 0.75){
    LV$col[row] <- as.character("#FF8C00")
  }
  else if (LV$HSM[row] > 0.5){
    LV$col[row] <- as.character(alpha("#FF8C00",alpha=0.75))
  }
  else if (LV$HSM[row] > 0.25){
    LV$col[row] <- as.character((alpha("#FF8C00",alpha=0.50)))
  }
  else if (LV$HSM[row] > 0.125){
    LV$col[row] <- as.character((alpha("#FF8C00",alpha=0.25)))
  }
  else if (LV$HSM[row] < 0.125){
    LV$col[row] <- as.character((alpha("#FF8C00",alpha=0.01)))
  }
}

  
pdf(file = paste0("LVLM_FST_438849_448241_v5.pdf"), width = 12, height = 6)
par(mar=c(5,4,4,5)+.1)
plot(NA,NA,xlim=c(438500, 449500),ylim=c(-1,1),axes=F,bty="n",xlab="",ylab="")
axis(side=1, at=c(438500, 438775, 444000, 449024, 449500), label=T)
axis(side=1, at=c(438500, 439500, 440500, 441500, 442500,443500,444500,445500,446500,447500,448500,449500), label=F)
mtext(side=1, text="Position on scaffold 3", cex=0, line=2, outer=FALSE)
Agrp2_Exon_01 <- rect(438849, 1, 438957, 1.05, density = NULL, angle = 45, col = "black", border = "black")
Agrp2_Exon_02 <- rect(443411, 1, 443485, 1.05, density = NULL, angle = 45, col = "black", border = "black")
Agrp2_Exon_03 <- rect(448409, 1, 449024, 1.05, density = NULL, angle = 45, col = "black", border = "black")
UTR_iso1 <- rect(438775, 1, 438848, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
UTR_iso2.1 <- rect(438835, 1, 438848, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
UTR_iso2.2 <- rect(438668, 1, 438711, 1.05, density = NULL, angle = 45, col = "darkgrey", border = "darkgrey")
mtext(side=3, text="agrp2", cex=1, line=0.75, outer=FALSE, at=443411, col = "black")
line <- rect(438800, 1.024, 448241, 1.025, density = NULL, angle = 45, col = "black", border = NA)
for (row in 1:nrow(LV)){
  rect(LV$Start[row], 0, LV$End[row], 1,col=LV$col[row],border=NA)
}
for (row in 1:nrow(LM)){
  rect(LM$Start[row], 0, LM$End[row], -1,col=LM$col[row],border=NA)
}
dev.off()


# get scaffold mean FST
Av_FST_M <- stats_M_nocov$HSM
Av_FST_M <- mean(Av_FST_M)
