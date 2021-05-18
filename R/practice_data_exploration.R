# packages
library(vcfR)
library(adegenet)
library(parallel)
library(ape)
library(tidyverse)
library(readxl)
library(report)

sink("./docs/package_versions.txt")
report_packages()
sink(NULL)

# Set file paths
vars_path <- "./data/bowtie_final_variants.vcf"
refseq_path <- "./data/kudzu_01Nov2018_Xb0hT.fasta"
# vars_path <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/Kudzu/bowtie_output/bowtie_final_variants.vcf"
# refseq_path <- "/uufs/chpc.utah.edu/common/home/u6033249/Data/Kudzu/kudzu_01Nov2018_Xb0hT_oneline.fasta.gz"


# import metadata and subset to correct samples ####

# read metadata excel sheet
meta <- read_xlsx("./data/Man3_Table1_5aug18.xlsx")
# clean Enumber column
meta$E_number <- meta$E_number %>% str_split("/") %>% map_chr(1) %>% str_trim()
meta


# read in DNA files ####
vcf <- read.vcfR(vars_path)
dna <- read.dna(refseq_path,format = "fasta")

# convert to genlight object
x <- vcfR2genlight(vcf)

# clean up names
indNames(x) <- indNames(x) %>% 
  str_remove_all("/scratch/kingspeak/serial/u6033249/GBS_Files/BOWTIE2/alignments/") %>% 
  str_remove_all(".fq.gz.sam.bam.sorted.bam")

# subset to only Enumbers in metadata
x <- x[indNames(x) %in% meta$E_number,]

# add lat lon country etc
x@other[["Lat"]] <- meta$Latitude
x@other[["Lon"]] <- meta$Longitude
x@other[["Variety"]] <- meta$Variety_genomics
x@other[["Country"]] <- meta$Country
x@other[["State_Province"]] <- meta$`State / Province`

saveRDS(x,"./data/genlight_object.RDS")

nInd(x) # number of individuals
nLoc(x) # number of loci
indNames(x) # names of the individuals
chromosome(x) # names of chromosomes (contigs) for each SNP
position(x) # positions of each SNP within the contig they were found in




# remove rare SNPs (detected fewer than 1000 times)
x2 <- x[,as.matrix(x) %>% colSums() >= 1000]
position(x2)





# REEEEALLY need to clean up variants somehow to get them to a reasonable count

# table of contig SNP counts (not all contigs same length of course!)
# locNames(x) %>% str_split("_") %>% map_chr(2) %>% str_split(";") %>% map_chr(1) %>% table() %>% plot()
# 
# pop(x) # grouping factor .... need to add this (population metadata?)
# other(x) # misc info stored as list
# NA.posi(x) # positions of missing values for each individual
# NA.posi(x) %>% unlist()

# chromosome(x) # names of chromosomes (contigs) for each SNP
# position(x) %>% plot# positions of each SNP within the contig they were found in

# position(x) %>% sqrt() %>% density(bw=10) %>% plot()
# 
# plot(density(position(x),bw=10))

# what if we took the reference genome and treated it like one BIIIG contig (remove all internal headers and line breaks)


snpposi.plot(position(x),genome.size = max(position(x))) +
  theme_minimal() +
  labs(caption = "Does this reflect the condition of the ref genome? (lots of short contigs)")
ggsave("$HOME/bowtie_snp_position_plot.png")
# snpposi.test(position(x),genome.size = max(position(x)))

# png("./output/glPlot.png")
# glPlot(x,posi = "topleft")
# dev.off()


# PCA
pca1 <- glPca(x2,nf = 6,parallel = TRUE)
pca1 <- readRDS("./output/bowtie_pca1.RDS")

scatter(pca1)
colorplot(pca1$scores,pca1$scores,transp = TRUE,cex=4,xlab="PC1",ylab="PC2")


# Tree
tre <- nj(dist(as.matrix(x2)))
# tre <- readRDS("./output/bowtie_tree.RDS")
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=1.5)
# abline(h=0,v=0, col="grey")
# add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)
plot(tre, typ="fan", show.tip=FALSE,cex=.5)
tiplabels(pch=20, col=myCol, cex=3)

km4 <- kmeans(pca1$scores,centers = 4)
km <- kmeans(pca1$scores,centers = 6)

pop(x2) <- km$cluster

# # DAPC
dapc4 <- dapc(x2,n.pca = 2,n.da = 1)
saveRDS(dapc4,"./output/dapc4.RDS")
?dapc

pred <- predict(dapc4)
# 
# 
# 8
# 2
scatter(dapc4,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,txt.leg=paste("group", 1:2), col=c("red","blue"))
compoplot(dapc4, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)
# 
loadingplot(dapc4$var.contr, thres=1e-3)
# 
# 
# 
# 
# x.dist <- dist(x)
# plot(x)
# 
# ?create.chromR
# 
# chrom <- create.chromR(vcf, name = "Pueraria_montana", seq=dna,verbose = TRUE)
# plot(chrom)
# chromoqc(chrom, dp.alpha = 50)
# chrom@var.info
# object.size(vcf)
# 
# 
# 
# 


