
# ALFDA-Affinity based Local Fisher Discriminant Analysis with network weights


#### Xinghu Qin & Peilin Jia, Beijing Institute of Genomics, Chinese Academy of Sciences & China National Center for Bioinformation.


ALFDA is based on a weighted graph with edge weights representing the genetical relatedness between two individuals. 

This package implements ALFDA algorithm. ALFDA is a nonlinear method that can capture both local and global structure for multimodal data. ALFDA searches the closest edges on the graph and then preserves the weights between and within populations with singular value decomposition, producing highly stable and reproduceable features.   

## Install the package from github:

```{R}
library(devtools)

install_github("xinghuq/ALFDA")

library("ALFDA")
```

## Basic examples using the simulated data


```{r}
### examples 

n_individuals <- 100
n_loci <- 1000
pop=rep(1:10,each=10)
pop=factor(pop,levels=unique(pop))
# Simulate genotype matrix
genotype <- matrix(sample(c(0, 1, 2), n_individuals * n_loci, replace = TRUE, prob = c(0.7, 0.2, 0.1)),
                   nrow = n_individuals, ncol = n_loci)

#In this code, we first set the number of individuals and loci using the variables n_individuals and n_loci, respectively.

#We then simulate the genotype matrix using the matrix() function. The sample() function randomly samples from the values 0, 1, and 2, with probabilities of 0.7, 0.2, and 0.1, respectively.

#The replace = TRUE argument indicates that sampling is done with replacement, and the nrow and ncol arguments specify the number of rows (individuals) and columns (loci) in the matrix.

#The resulting genotype matrix has 100 rows (individuals) and 1000 columns (loci), with each cell containing a randomly sampled genotype value of 0, 1, or 2.

allele_freq=apply(genotype, 2, function(x){x/2})

Rd=ALFDA(allele_freq,y=pop, r=3, kaf = 10, sigma = 0.5, knn = 6, reg = 0.001))


```

### Ancient genome examples

```{R}
## Obtaining the genomic data from the European Nucleotide Archive under accession number PRJEB47891 and preprocess them to get vcf or gds format

genofile <- snpgdsOpen("brit.gds")
read.gdsn(index.gdsn(genofile, "sample.id"))

read.gdsn(index.gdsn(genofile, "snp.id"))


#SNP_posit_UKM=read.table("brit.bim",header=FALSE)
### MAF 0.01
snpset <- snpgdsLDpruning(genofile, sample.id=as.character(UKM_geo_ID$sample.id),autosome.only=TRUE,remove.monosnp=TRUE,ld.threshold=0.2,maf = 0.01,missing.rate=0.5)
# Get all selected snp id 826 samples, 130,577 SNPs using maf 0.01, and 99,154 markers with MAF=0.05
snpset.id <- unlist(unname(snpset))
save(snpset.id,file="prundingLD_0.2_MAF_0.01_missingrate0.5_rm_monosnp_snpID.RData")


## LD pruning,MAF >0.01, individual missing.rate=0.01,remove monsnp
genofile <- snpgdsOpen("brit.gds")
read.gdsn(index.gdsn(genofile, "sample.id"))
snpset <- snpgdsLDpruning(genofile, sample.id=as.character(UKM_geo_ID$sample.id),autosome.only=TRUE,remove.monosnp=TRUE,ld.threshold=0.2,maf = 0.01,missing.rate=0.5)
snpset.id <- unlist(unname(snpset))
##### get sample ID, SNP ID
g <- snpgdsGetGeno(genofile, snp.id = snpset.id, with.id=TRUE)
colnames(g$genotype)=as.vector(g$snp.id)
rownames(g$genotype)=as.vector(g$sample.id)

snpgdsClose(genofile)


###

UKM_AFLDA_pop0.5=AFLDA(as.matrix(g),y=UKM_sample_new$Country_sub,r=10,kaf=10,knn=2,metric = "plain")
write.csv(as.data.frame(UKM_AFLDA_pop0.5$Z),file="rd_UKM_AFLDA_popk10sigam0.5.csv")


### visualize the genetic sturcture

library(plotly)
cols1=rainbow(length(unique(UKM_sample_new$Country_sub)))
p_2D <- plot_ly(as.data.frame(UKM_AFLDA_pop0.5$Z), x =UKM_AFLDA_pop0.5$Z[,1], y =UKM_AFLDA_pop0.5$Z[,1],  color =UKM_sample_new$Country_sub, colors=cols1,symbol=UKM_sample_new$Country_sub,symbols = unique(popinf$pop)) %>% 
  layout(autosize = TRUE)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X 1'),
                      yaxis = list(title = 'X 2')))

```
### Citations
Qin.X, Jia.,P. 2023. New machine learning method identifies subtle fine-scale genetic stratification in diverse populations. Submitted.

#### Contact
Email: qinxinghu@gmail.com
 

