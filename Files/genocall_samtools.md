
The following commands may not be fully compatible with more recent versions of SAMtools and therefore should be considered as pure indication.

## SAMtools

We now use SAMtools to call genotypes on the human dataset.
Firstly, we generate a VCF file.
```
samtools-1.2/samtools mpileup -f input/human/hg19_chr1.fa -b input/human/bams.list -u -r 1: | bcftools-1.2/bcftools view - > output/human.vcf
```
We can use **vcftools** to extract genotypes
```
vcftools_0.1.12b/bin/vcftools --vcf output/human.vcf --012 --out output/human.sam
```
Let us have a look at the output files
```
less -S output/human.sam.012
less output/human.sam.012.pos
```
Again, genotypes are coded as 0, 1, and 2.

------

We want to compare results obtained with ANGSD (using a prior or not) and SAMtools.
If not stored, recompute genotypes with ANGSD with and without (or uniform) prior:
```
./angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa.gz $OPT1 $OPT2 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01  -r 1: -doGeno 2 -doPost 2 -out output/human.prior -postCutoff 0
./angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa.gz $OPT1 $OPT2 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01  -r 1: -doGeno 2 -doPost 1 -out output/human.unif -postCutoff 0
./angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa.gz $OPT1 $OPT2 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01  -r 1: -doGeno 2 -doPost 2 -out output/human.prior2 -postCutoff 0.8
# in the latter case we set low-confidence assignments as missing data
```

Let us use R to quickly parse these files. Open R:
```
# read files;
 angsd.prior <-  read.table("output/human.prior.geno.gz");
 angsd.unif <-  read.table("output/human.unif.geno.gz");
 sam.pos <-  read.table("output/human.sam.012.pos")[,2];
# get overlapping sites;
 tab <- table(c(angsd.prior[,2],angsd.unif[,2],sam.pos));
 pos <- as.numeric(names(tab[tab==3]));
# extract genotypes
 angsd.prior<-as.data.frame(t(angsd.prior[angsd.prior[,2] %in% pos,][,-c(1,2)]))
 angsd.unif<-as.data.frame(t(angsd.unif[angsd.unif[,2] %in% pos,][,-c(1,2)]))
 sam<-read.table("output/human.sam.012")[,-1][sam.pos %in% pos]
# how many genotypes differ?
 # between SAMtools and ANGSD with prior
 sum(sam!=angsd.prior)
 # between SAMtools and ANGSD
 sum(sam!=angsd.unif)
 # between ANGSD with and without prior
 sum(angsd.prior!=angsd.unif)
# what is the proportion of genotypes that differ?
 # between SAMtools and ANGSD with prior
 sum(sam!=angsd.prior)/prod(dim(sam))
 # between SAMtools and ANGSD
 sum(sam!=angsd.unif)/prod(dim(sam))
 # between ANGSD with and without prior
 sum(angsd.prior!=angsd.unif)/prod(dim(sam))
# and if we do not consider missing sites (the one with low probability)?
 angsd.prior2 <-  read.table("output/human.prior2.geno.gz");
 angsd.prior2 <-as.data.frame(t(angsd.prior2[angsd.prior2[,2] %in% pos,][,-c(1,2)]))
 angsd.prior2[angsd.prior2==-1]<-NA
 # between SAMtools and ANGSD with prior
 sum(sam!=angsd.prior2, na.rm=T)/length(which(!is.na(angsd.prior2)))
 # between ANGSD with and without prior
 sum(angsd.prior2!=angsd.unif, na.rm=T)/length(which(!is.na(angsd.prior2)))
# what is the proportion of genotyeps with data? 
 length(which(!is.na(angsd.prior2)))/prod(dim(sam))
```

As a general recommendation, in case of low depth data, using a prior based on allele frequencies has been shown to well perform in genotype calling.




