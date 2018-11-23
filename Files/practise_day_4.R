
## This script shows an example of how to detect natural positive selection in SNP data

# Data set is comprised of ~ 200,000 SNPs from the metaboChip in HapMap populations: YRI (Africa), CEU (Europe), CHB (East Asian)

# Data is in the widely-used plink format (http://pngu.mgh.harvard.edu/~purcell/plink/), specifically in allele frequency per cluster format

# if you want to have a look at the input file:
# less -S ../Data/hapmap.frq.strat

# ----

# one of the most powerful approaches to detect selection is to compare genetic variation between 2 or more populations, assuming that the control populations are neutrally evolving for that specific gene

# FST is a very commonly used metric to identify changes in allele (or haplotype) frequencies between populations

# Let us assume that we want to detect selection in Europeans (CEU).
# Therefore compute FST between CEU and a reference/control population (e.g. Africans YRI).

# we can compute per-site FST values using a method-of-moments estimator
# less -S ../Scripts/plink2fst.R
#Rscript ../Scripts/plink2fst.R

# this generates a file Results/hapmap.fst
# it may take a while so you can copy it from ../Data

# cp ../Data/hapmap.fst Results/.
# less -S Results/hapmap.fst

# -----

# it is convenient to plot the FST values across the genome to identify potential outliers

# produce a Manhattan plots
#Rscript ../Scripts/manPlotFST.R Results/hapmap.fst Results/hapmap.fst.jpg

# open Results/hapmap.fst.jpg

# -----

# now open R

# high values of FST indicate a change in allele frequencies between CEU and YRI
# however, we cannot know whether either YRI or CEU is under selection for these high-FST values

# one strategy would be to calculate an additional statistic and look, for instance, whether loci with high FST are associated with a decrease in nucleotide diversity in YRI or CEU.

# an alternative is to compute a measure of genetic differentiation that uses a 3rd population

# let us compute the Population Branch Statistic (PBS)

# read the FST values
fst=read.table("Results/hapmap.fst", stringsAsFact=F, header=T)

par(mfrow = c(1, 3))
hist(fst$MAF.YRI)
hist(fst$MAF.CEU)
hist(fst$MAF.CHB)

sfs_ceu <- fst$MAF.CEU
sfs_ceu[which(sfs_ceu>0.5)] <- 1 - sfs_ceu[which(sfs_ceu>0.5)]  

par(mfrow = c(1,1))
hist(sfs_ceu)

# plot the results (PBS values) across the genome
# PBS(CEU) = ( T(YRI.CEU) + T(CEU.CHB) - T(YRI.CHB) ) / 2

# T could be computed as -log(1-FST)

# QUESTION:
# calculate PBS values for all SNPs towards testing selection in CEU

# ANSWER:
# YRI.CEU + CEU.CHB - YRI.CHB
pbs= ( (-log(1-fst$FST.YRI.CEU)) + (-log(1-fst$FST.CEU.CHB)) - (-log(1-fst$FST.YRI.CHB)) ) / 2
pbs[which(pbs<0)]=0 # set to 0 negative values, for convenience

# a Manhattan-like plot
cols=rep("grey", nrow(fst)); cols[which( (fst$chrom %% 2) == 1)]="lightgrey"

# -----

# how would you identify outliers for PBS?

# first, calculate and plot empirical thresholds (e.g. 99th and 99.9th)
pbs_th<-quantile(pbs, seq(0.99,0.999,0.001), na.rm=T)[c(1,10)]

plot(x=fst$cpos, y=pbs, col=cols, frame=F, xlab="", xaxt="n", ylab="PBS", main="CEU", pch=16)
abline(h=pbs_th, lty=2)

# -----

# check which SNP is the top hit
fst[which.max(pbs),]

# QUESTION:
# Where is this SNP located? In which gene? What are the allele frequencies in human populations? Is the derived or ancestral allele at high frequency in Europeans? Does it show any other signature of selection?

# ANSWER:
# https://genome-euro.ucsc.edu
# http://haplotter.uchicago.edu/
# http://hgdp.uchicago.edu/cgi-bin/gbrowse/HGDP/

# QUESTION:
# Is it an already known locus under positive selection in Europeans?
# If so, what is its phenotypic effect? What is the most likely selective pressure?

# ANSWER:
# http://www.ncbi.nlm.nih.gov/pubmed


# -----

# Let us consider the second best candidate SNP based on PBS analysis

fst[which.max(pbs[-which.max(pbs)]),]

# again, let's check its annotation using, for instance, the UCSC Genome Browser
# rs: rs482000
# gene: SLC35F3

# now we want to assign a p-value to this value of PBS by using neutral simulations

# we will use the program "ms" (be sure to have it installed) to simulate allele frequencies in YRI, CEU and CHB under a model of neutral evolution

# first thing, we need to decide a demographic model for our human populations, and assume it is the one proposed in this paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695

source("../Scripts/functions.R")

# define the directory where you have installed "ms"
ms_dir<-"/home/BIO5789/utils/bin/ms"

ms.command <- paste(ms_dir, "326 10000 -s 1 -I 3 118 120 88 -n 1 1.68 -n 2 1.12 -n 3 1.12 -eg 0 2 72 -eg 0 3 96 -ma x 2.42 1.52 2.42 x 7.73 1.52 7.73 x -ej 0.029 3 2 -en 0.029 2 0.29 -ej 0.19 2 1 -en 0.30 1 1 | gzip > Results/ms.txt.gz")

system(ms.command, intern=F)

# this will create a file with 10,000 neutral simulations

# read the output
sim.chroms=readMs("Results/ms.txt.gz" , 326)$hap

# below please find enclosed a possible (slow) solution to compute PBS for each simulated SNP, which will be recorded in an array called sim.pbs

nreps=length(sim.chroms)

# initialise
sim.pbs<-rep(NA, nreps)
for (i in 1:nreps) {

	pops=list()
        pops[[1]]=sim.chroms[[i]][1:118] # YRI
        pops[[2]]=sim.chroms[[i]][119:238] # CEU
        pops[[3]]=sim.chroms[[i]][239:326] # CHB

	# compute FST
        sim.fst=chroms2fst(pops)

	# compute PBS
        sim.pbs[i]=( (-log(1-sim.fst[1])) + (-log(1-sim.fst[2])) - (-log(1-sim.fst[3])) ) / 2

	# assign p-value
        if ((i %% 100)==0) cat(i, "\t", length(which(sim.pbs>=1.42))/i, "\n")

}
sim.pbs[which(sim.pbs<0)]=0

# ---

# QUESTION:
# from these simulations, compute a p-value to test whether this SNP is (likely to be) under positive selection in Europeans

# ANSWER:
length(which(sim.pbs>=1.42)) / length(sim.pbs)
# or
length(which(sim.pbs>=1.42)) / length(which(!is.na(sim.pbs)))

# plot
hist(sim.pbs, main="Simulations under neutrality", xlab="PBS", breaks=20)
abline(v=1.42, lty=2)


# EXERCISE
# Run at least 100,000 simulations and assess the significance of this PBS value by jointly considering its DAF (derived allele frequency) and PBS value. In other words, record the significance intervals on a DAF-PBS Cartesian plot and check whether this SNP is indeed an outlier


