
**WORKFLOW**:
... FILTERED DATA > GENOTYPE AND SNP CALLING > POPULATION GENETICS (SFS)

Another important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described). 
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site. Finally, an estimate of the SFS is computed.

Sequence data -> Genotype likelihoods -> Posterior probabilities of SAF -> SFS

These steps can be accomplished in ANGSD using `-doSaf 1` option and the program `realSFS`, as following:
```
ngsTools/angsd/angsd -b input/lyca/bams.list -anc input/lyca/referenceseq.fasta -sites input/lyca/sites.angsd.bed -minMapQ 10 -minQ 10 -GL 1 -doSaf 1 -out output/lyca
```
What are the output files?
```
ls output/lyca.*
# output/lyca.arg		output/lyca.saf		output/lyca.saf.pos.gz
```

File saf.pos summarizes the genomic coordinates of all sites analyzed. 
```
gunzip -c output/lyca.saf.pos.gz | head
```
File .saf is a binary file, so we use a quick R script to have a look at its contents:
```
Rscript scripts/printSAF.R output/lyca.saf 20 99085 0
```

What do these values represent?
Let us convert them in proper probabilities.
```
Rscript scripts/printSAF.R output/lyca.saf 20 99085 1
```

Now we can estimate the SFS:
```
ngsTools/angsd/misc/realSFS output/lyca.saf 40 -P 2 -nSites 1000000 2> /dev/null > output/lyca.sfs 
cat output/lyca.sfs
```

Again `.sfs` file is in log-scale. We can plot it.
```
Rscript -e 'sfs=as.numeric(scan("output/lyca.sfs",quiet=T)); sfs=exp(sfs)/sum(exp(sfs));pdf(file="output/lyca.sfs.pdf"); par(mfrow=c(1,2));barplot(sfs,names=0:40,main="SFS"); sfs[1]=NA;sfs=sfs/sum(sfs,na.rm=T); barplot(sfs[2:41],names=1:40,main="SFS",sub="only variable sites");dev.off()';
```
It this does not work type:
```
Rscript scripts/DoPlotSFS.R
```
or open R and type:
```
sfs=as.numeric(scan("output/lyca.sfs",quiet=T));
sfs=exp(sfs)/sum(exp(sfs));
pdf(file="output/lyca.sfs.pdf");
par(mfrow=c(1,2));
barplot(sfs,names=0:40,main="SFS");
sfs[1]=NA;
sfs=sfs/sum(sfs,na.rm=T);
barplot(sfs[2:41],names=1:40,main="SFS",sub="only variable sites");
dev.off()
```

```
open output/lyca.sfs.pdf # on mac
```
Why is it a bit bumpy?

The method used here does not rely on called genotypes and it is suitable for low-depth data.
Let us explore the different behaviour using different options. 

We can estimate the SFS using different genotype likelihoods or genotype priors.

```
ngsTools/angsd/angsd -b input/lyca/bams.list -ref input/lyca/referenceseq.fasta -sites input/lyca/sites.angsd.bed -minMapQ 10 -minQ 10 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -doGeno 2 -doPost 1 -out output/lyca.angsd.v11
ngsTools/angsd/angsd -b input/lyca/bams.list -ref input/lyca/referenceseq.fasta -sites input/lyca/sites.angsd.bed -minMapQ 10 -minQ 10 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -doGeno 2 -doPost 2 -out output/lyca.angsd.v12
ngsTools/angsd/angsd -b input/lyca/bams.list -ref input/lyca/referenceseq.fasta -sites input/lyca/sites.angsd.bed -minMapQ 10 -minQ 10 -GL 2 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -doGeno 2 -doPost 1 -out output/lyca.angsd.v21
ngsTools/angsd/angsd -b input/lyca/bams.list -ref input/lyca/referenceseq.fasta -sites input/lyca/sites.angsd.bed -minMapQ 10 -minQ 10 -GL 2 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -doGeno 2 -doPost 2 -out output/lyca.angsd.v22
```

Let us plot the results using R. Open it and type
```
# function to get the SFS
getSFS <-function(x){
  a<-read.table(x)[,-c(1:2)]
  ns=ncol(a)
  a[a==-1] <- NA;
  tot=as.numeric(table(rowSums(a,na.rm=T)))
  tot/sum(tot)
}

# read file
sfs=rbind(getSFS("output/lyca.angsd.v11.geno.gz"), getSFS("output/lyca.angsd.v21.geno.gz"), getSFS("output/lyca.angsd.v12.geno.gz"), getSFS("output/lyca.angsd.v22.geno.gz"))

# let us merge the previous Maximum Likelihood Estimate of SFS
mle=as.numeric(scan("output/lyca.sfs",quiet=T));
mle=exp(mle[2:41])/sum(exp(mle[2:41]));
sfs=rbind(sfs, mle)

# plot
pdf(file="output/lyca.many.sfs.pdf")
barplot(sfs, beside=T, legend=c("SAMtools, HWE-prior","GATK, HWE-prior", "SAMtools, no prior", "GATK, no prior", "MLE"))
dev.off()
pdf(file="output/lyca.many.zoom.sfs.pdf")
barplot(sfs[,1:7], beside=T, legend=c("SAMtools, HWE-prior","GATK, HWE-prior", "SAMtools, no prior", "GATK, no prior", "MLE"))
dev.off()
```

Open the plots
```
open output/lyca.many.sfs.pdf output/lyca.many.zoom.sfs.pdf # on mac
# evince output/lyca.many.sfs.pdf output/lyca.many.zoom.sfs.pdf 
```
Do you get similar or different estimates of SFS?


**EXERCISE**
Compute the SFS on the same dataset after calling SNPs and genotypes, using either ANGSD or SAMtools, by using different thresolds for SNP calling.
Compare the results.
A possible **SOLUTION** for this exercise can be found [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/solutions.txt).
Alternatively, estimate the SFS on your dataset or on the example human dataset.


**ADDITIONAL MATERIAL**
In case of low-depth sequencing data, summary statistics can be computed by keeping statistical uncertainty into account.
We provide some examples on how to compute summary statistics from low-depth sequencing data using [ANGSD](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/lowcov.md) and [ngsTools](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).



