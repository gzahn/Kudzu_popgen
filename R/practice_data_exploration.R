# packages
library(vcfR)
library(adegenet); packageVersion("adegenet")
library(parallel)
library(ape)
library(tidyverse)

# read in files
vcf <- read.vcfR("./data/bowtie_final_variants.vcf")
dna <- read.dna("./data/kudzu_01Nov2018_Xb0hT.fasta",format = "fasta")

readLines("./data/mappedreadcounts.txt") %>% as.numeric() %>% plot()


x <- vcfR2genlight(vcf)
ploidy(x) %>% table()

# in a genlight object: rows=individuals, cols=SNPs
x[1,2] %>% as.matrix()

nInd(x) # number of individuals
nLoc(x) # number of loci
indNames(x) # names of the individuals

# clean up names
indNames(x) <- indNames(x) %>% 
  str_remove_all("/scratch/kingspeak/serial/u6033249/GBS_Files/BOWTIE2/alignments/") %>% 
  str_remove_all(".fq.gz.sam.bam.sorted.bam")

locNames(x) # names of loci (contigs from draft genome)


# REEEEALLY need to clean up variants somehow to get them to a reasonable count

# table of contig SNP counts (not all contigs same length of course!)
locNames(x) %>% str_split("_") %>% map_chr(2) %>% str_split(";") %>% map_chr(1) %>% table() %>% plot()

pop(x) # grouping factor .... need to add this (population metadata?)
other(x) # misc info stored as list
NA.posi(x) # positions of missing values for each individual
NA.posi(x) %>% unlist()

chromosome(x) # names of chromosomes (contigs) for each SNP
position(x) # positions of each SNP within the contig they were found in



plot(density(position(x),bw=10))

snpposi.plot(position(x),genome.size = max(position(x))) +
  theme_minimal() +
  labs(caption = "Does this reflect the condition of the ref genome? (lots of short contigs)")
ggsave("output/snp_position_plot.png")
# snpposi.test(position(x),genome.size = max(position(x)))

png("./output/glPlot.png")
glPlot(x,posi = "topleft")
dev.off()


# make small subset of genlight for testing
x <- x[1:10,1:1000]

# PCA
pca1 <- glPca(x,nf = 3,parallel = TRUE)
saveRDS(pca1,"./output/pca1.RDS")
pca1 <- readRDS("./output/pca1.RDS")
pca1$scores
scatter(pca1)
beepr::beep(sound=8)

colorplot(pca1$scores,pca1$scores,transp = TRUE,cex=4,xlab="PC1",ylab="PC2")
abline(h=0,v=1.5,col="grey")

# Tree
tre <- nj(dist(as.matrix(x)))

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)
dev.off()
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)


# DAPC
dapc1 <- dapc(x,pop = c(rep("pop1",5),rep("pop2",5)),n.pca = 5,n.da = 1)


8
2
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,txt.leg=paste("group", 1:2), col=c("red","blue"))
compoplot(dapc1, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)

loadingplot(dapc1$var.contr, thres=1e-3)




x.dist <- dist(x)
plot(x)

?create.chromR

chrom <- create.chromR(vcf, name = "Pueraria_montana", seq=dna,verbose = TRUE)
plot(chrom)
chromoqc(chrom, dp.alpha = 50)
chrom@var.info
object.size(vcf)

