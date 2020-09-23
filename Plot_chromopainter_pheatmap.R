### Heatmap for Agrp2 haplotype paper, Sabine Urban
# 1.July.2019

### Load Pheatmap and other packages
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}; require(stringr)  
if("pheatmap" %in% rownames(installed.packages()) == FALSE) {install.packages("pheatmap")}; require(pheatmap)  
if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}; require(plyr)  
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}; require(RColorBrewer)

setwd("//") # Set working directory
### Load input, haplotypes are in separate tables (one file per haplotype with headers: "position donor1 donor2 donor3 donor4")
# Uses chromopainter output (*copyprobsperlocus.out)
donor_hap1<-"Asta2_xxx_HAP1_LV_nst" # Prefix of all files that I want to read
csv.dir<-"//"
files<-list.files(csv.dir)[str_sub(list.files(csv.dir),1,21)==donor_hap1] # 1 to 21 is the position/length of the prefix
input <- paste0(csv.dir,files)

data.list<-list()
dataset1<-""

for (i in 1:length(input)) {
  dataset <- read.csv(input[i], sep=" ")
  Recipient=str_sub(as.character(input[i]),123,126)
  header <- c("pos", paste0("H. sauvagei","_",Recipient,"_","hap1"), paste0("A.stappersii","_",Recipient,"_","hap1"), paste0("H.vittatus","_",Recipient,"_","hap1"), paste0("H.gracilior","_",Recipient,"_","hap1")) # change header to more meaningful names
  names(dataset) <- header
  dataset <- dataset[with(dataset,  ((pos %in% 438000:449000))), ]
  dataset1 <- merge(dataset1,dataset)
}

dataset1 <- subset( dataset1, select = -x )
##########################################################################################
donor_hap2<-"Asta2_xxx_HAP2_LV_nst"
csv.dir<-"//"
files<-list.files(csv.dir)[str_sub(list.files(csv.dir),1,21)==donor_hap2]
input <- paste0(csv.dir,files)

data.list<-list()

dataset2<-""

for (i in 1:length(input)) {
  dataset <- read.csv(input[i], sep=" ")
  Recipient=str_sub(as.character(input[i]),123,126)
  header <- c("pos", paste0("H. sauvagei","_",Recipient,"_","hap2"), paste0("A.stappersii","_",Recipient,"_","hap2"), paste0("H.vittatus","_",Recipient,"_","hap2"), paste0("H.gracilior","_",Recipient,"_","hap2"))
  names(dataset) <- header
  dataset <- dataset[with(dataset,  ((pos %in% 438000:449000))), ]
  dataset2 <- merge(dataset2,dataset)
}

dataset2 <- subset( dataset2, select = -x )
##########################################################################################

merge<-cbind(dataset1,dataset2[,-1]) #combine both tables and remove first column from second input (redundant)
merge.red <- merge[with(merge,  ((pos %in% 438000:449000))), ]  # Take a subset of the data (here region of interest)... otherwise also too much for pheatmap
merge.red<-merge.red[ , order(names(merge.red))] # Order by name
merge.red <- data.frame(append(merge.red, c(" "=NA), after=6))
merge.red <- data.frame(append(merge.red, c(" "=NA), after=13))
merge.red <- data.frame(append(merge.red, c(" "=NA), after=20))
rownames(merge.red)<-merge.red[,28] #Move position information into the rownames
merge.red<-merge.red[,-28] # remove the position column
#merge.red<-merge.red[nrow(merge.red):1,] # Change order ... pheatmap will invert it
merge.red<-t(merge.red) #transpose the matrix to have position on x axis

### define annotation
annotation.col<-list(Annotation = c(Enhancer="grey",Exonic="black",' ' = "white"))

### define where enhancer and gene is
Enh<-rep(" ",ncol(merge.red))
Enhancer.Start<-442320
Enhancer.End<-443410
Exon1.Start<-438849
Exon1.End<-438957
Exon2.Start<-443411
Exon2.End<-443485
Exon3.Start<-448241
Exon3.End<-448408
Enh[which(colnames(merge.red)>=Enhancer.Start & colnames(merge.red)<=Enhancer.End)]<-"Enhancer"
Enh[which(colnames(merge.red)>=Exon1.Start & colnames(merge.red)<=Exon1.Start)]<-"Exonic"
Enh[which(colnames(merge.red)>=Exon2.Start & colnames(merge.red)<=Exon2.End)]<-"Exonic"
Enh[which(colnames(merge.red)>=Exon3.Start & colnames(merge.red)<=Exon3.End)]<-"Exonic"

annotation<-data.frame(Annotation = factor(Enh))
rownames(annotation)<-colnames(merge.red)
palette_V<-brewer.pal(9,"Oranges")

pdf(file = paste0("Asta2_heatmap_nst_recipient.pdf"), width = 12, height = 6)
pheatmap(merge.red,
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = T,
         cellwidth=5,
         cellheight=5,
         fontsize=5,
         #width= 120,
         #height=90,
         annotation_col = annotation,
         annotation_colors = annotation.col,
         na_col = "white",
         color=palette_V,
         border_color = NA)
dev.off()
