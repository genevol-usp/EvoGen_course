
The commands below are an illustration of how to impute genotypes using BEAGLE. This script do not work with the current files.

Here we will briefly show how to impute data from genotype likelihoods using **BEAGLE** software.
These commands are not supposed to be run on the example datasets.
They are just to illustrate how imputation can be done (lines by T. Korneliussen).

First we need to prepare the information needed for imputation, and we will use ANGSD for this.
With ANGSD we will call variable sites (-min_pval) based on allele frequencies (-doMaf) computed from the reference allele (-ref) and the inferred alternative allele (-doMajorMinor). 
Based on these sites we will output the genotype likelihoods based on the SAMtools model (-GL 1) in the BEAGLE software format (-doGlf 2).

This is the command to generate such files:
```
ngsTools/angsd/angsd -bam bams.list -doSNP 1 -minLRT 6.6  -doMaf 2 -doMajorMinor 4 -r 1: -doGlf 2  -GL 1 -out angsd  -ref ref/hg19.fa.gz
```
You can look at the BEAGLE genotype likelihood format:
```
less -S angsd.beagle.gz
```
What are these 3 columns we have for each individuals?

Now we run the imputation based on the genotype likelihoods:
```
java -jar beagle.jar like=angsd.beagle.gz out=imputation
```

---------

The 33 individuals sequenced as part of the 1000 genomes project was also part of the HapMap project. 
Therefore, we have accurate information about many of the genotypes from these individuals. 

You can extract genotype information from hapmap website (e.g. hapmap_CEU_r23a_filteredHg19Chr1.tar.gz).
Unzip these files with known genotypes.
```
tar -xf hapmap_CEU_r23a_filteredHg19Chr1.tar.gz
```
Return to the previous folder and open R.

You need R package `snpMatrix` to run this example
```
# for linux
wget http://popgen.dk/albrecht/open/snpMatrix_1.14.6.tar.gz
# for mac
curl http://popgen.dk/albrecht/open/snpMatrix_1.14.6.tar.gz
# install the package from the terminal
R CMD INSTALL snpMatrix_1.14.6.tar.gz
```

Open R.
```
library(snpMatrix) # from Anders Albrechtsen

#read in the genotype data. Change the path to where you unpacked the genotypes
pl<-read.plink("data/exercise4/hapmap_CEU_r23a_filteredHg19Chr1")
bim<-read.table("data/exercise4/hapmap_CEU_r23a_filteredHg19Chr1.bim",as.is=T)

#read in the called genotypes from angsd, samtools and GATK just as before
gatk.pos <- read.table("gatk.012.pos")[,2]
sam.pos <-  read.table("sam.012.pos")[,2]
angsd <-  read.table("angsd.geno.gz") #you might have to remove the .gz
ta<-table(c(angsd[,2],gatk.pos,sam.pos))
pos <- as.numeric(names(ta[ta==3]))
angsd<-as.data.frame(t(angsd[angsd[,2] %in% pos,][,-c(1,2)]))
gatk<-read.table("gatk.012")[,-1][gatk.pos %in% pos]
sam<-read.table("sam.012")[,-1][sam.pos %in% pos]

#find the overlab with the hapmap data
table(keep<-pos%in%bim[,4])
int<-which(bim[,4]%in%pos)

indNames<-rownames(pl) # individual names in HapMap
indNamesBam<-sapply(strsplit(sub("small","",basename(scan("bams.list"
,what="theFck"))),".m"),function(x)x[1]) #individual names in the 33 1000genome

#convert the genotype data into integers
genoHap<-as.integer(pl[,int])
genoHap<-matrix(3-genoHap,nrow=60,ncol=length(int))
rownames(genoHap)<-indNames
genoHap[genoHap==3]<-NA

# match the strands
mafs<-read.table("angsd.mafs.gz",as.is=T,head=T) #you might have the remove the .gz
mafs<-mafs[mafs$position%in%pos,]
refStrand<-bim[int,5]==mafs$ref[keep]
refStrand<-bim[int,5]==mafs$ref[keep]
genoHap[rep(refStrand,each=60)]<-2-genoHap[rep(refStrand,each=60)]

#estimate the mean concordance rate
sapply(list(sam=sam,gatk=gatk,angsd=angsd),function(x)
mean(genoHap[indNamesBam,]==as.matrix(x[,keep]),na.rm=T))

#the above calculations assumed that the missing genotypes are discordance.
# Estimate the concordance rate without missing data
angsdNA<-as.matrix(angsd)
angsdNA[angsd==-1]<-NA
gatkNA<-as.matrix(gatk)
gatkNA[gatkNA==-1]<-NA
# do not type \ symbol
sapply(list(sam=sam,gatk=gatk,angsd=angsd,gatkNA=gatkNA,angsdNA=angsdNA),function(x) \
mean(genoHap[indNamesBam,]==as.matrix(x[,keep]),na.rm=T))
```
Do not close R.

Which of the genotyping methods performs the best in the above analysis?
How could the comparison of methods be improved?
We now include the imputation results in the comparisons

```
#read in the genotype probabilities
gprobsDat<-read.table("imputation.angsd.beagle.gz.gprobs.gz",head=T,as.is=T)
impuPos<-as.integer(sub("1_","",gprobsDat[,1])) #position

#keep only overlapping positions and convert to a matrix
gprob<-as.matrix(gprobsDat[impuPos %in% pos[keep],-c(1:3)])

#call the genotype with the highest probability
funn<-function(x)
  apply(gprob[,(x-1)*3+1:3],1,which.max)-1

imputa <- t(sapply(1:33, funn))

#get the mean concordance rate
mean(genoHap[indNamesBam,]==imputa,na.rm=T)
```

Why does this method perform so much better?
How could you increase the accuracy even further?
When will this method not increase performance?
How do you think the performance difference correlates with sample allele frequency?





