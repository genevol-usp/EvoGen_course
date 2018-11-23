
We will use now ngsTools, a set of programs to estimate several quantities of interest in population genetics, especially in the case of multiple populations.
Its implementation is described [here](http://www.ncbi.nlm.nih.gov/pubmed/24458950).
We will use the butterfly dataset, since it contains several populations.
Our first aim is to compute some measures of genetic differentiation, as described [here](http://www.ncbi.nlm.nih.gov/pubmed/23979584).

------------

We need to compute probabilities of allele frequencies for 2 populations we are considering now separately.
We will consider 'sin' and 'vic' populations, with 7 and 6 individuals each.
```
ngsTools/angsd/angsd -b input/lyca/bams_sin.list -anc input/lyca/referenceseq.fasta -GL 1 -doSaf 1 -out output/lyca.sin -sites input/lyca/sites.angsd.bed
ngsTools/angsd/angsd -b input/lyca/bams_vic.list -anc input/lyca/referenceseq.fasta -GL 1 -doSaf 1 -out output/lyca.vic -sites input/lyca/sites.angsd.bed
```

Let us check their number of sites:
```
gunzip -c output/lyca.sin.saf.pos.gz | wc -l
#  99001
gunzip -c output/lyca.vic.saf.pos.gz | wc -l
#  98907
```
Of the original 100,000 sites, some positions are skipped by the basic internal filtering done by ANGSD.
Therefore, we first need to extract common sites in the two files, and recreate new .saf files.
This can be done with these commands (from F.G. Vieira):
```
gunzip -c output/lyca.sin.saf.pos.gz > output/A.saf.pos
gunzip -c output/lyca.vic.saf.pos.gz > output/B.saf.pos
awk 'FNR==NR {x[$1"_"$2]=NR; next} x[$1"_"$2] {print x[$1"_"$2]; print FNR > "/dev/stderr"}' output/A.saf.pos output/B.saf.pos >output/A.pos 2>output/B.pos
rm output/A.saf.pos output/B.saf.pos
```
Now we have the indexes of common sites, indeed they have the same number of sites:
```
wc -l output/A.pos 
#  98862
wc -l output/B.pos 
#  98862
```

Finally we create new .saf files with only sites in common:
```
ngsTools/ngsUtils/GetSubSfs -infile output/lyca.sin.saf -posfile output/A.pos -nind 7 -nsites 99001 -len 98862 -outfile output/lyca.sin.fix.saf
ngsTools/ngsUtils/GetSubSfs -infile output/lyca.vic.saf -posfile output/B.pos -nind 6 -nsites 98907 -len 98862 -outfile output/lyca.vic.fix.saf
```
A quick trick to check that everything went fine is to retrieve the dimension of each file (using "ls -l"), then divide this number by 8 and then by the double of the number of individuals plus 1.
You should get the final number of sites (e.g. `ls -l output/lyca.sin.fix.saf` gives 11863440 which is 11863440/8/15=98862).

-------------

Genotype likelihoods > Probabilities of Sample Allele Frequency > 2D-SFS

Now we can estimate a 2D-SFS which will be used as prior in our calculation of genetic differentiation:
```
ngsTools/ngsPopGen/ngs2dSFS -postfiles output/lyca.sin.fix.saf output/lyca.vic.fix.saf -outfile output/lyca.2dsfs -relative 1 -nind 7 6 -maxlike 1 -nsites 98862 -block_size 50000 -islog 1
cat output/lyca.2dsfs
```
We can even plot it using R (you need ggplot2 installed)
```
Rscript ngsTools/ngsPopGen/scripts/plot2dSFS.R output/lyca.2dsfs output/lyca.2dsfs.pdf sin vic
open output/lyca.2dsfs.pdf
```

A better way would be to use posterior probabilities and thus estimate the per-species SFS to be used as prior.

---------

... > 2D-SFS > FST <br>

Finally, we are able to estimate per-site FST:
```
ngsTools/ngsPopGen/ngsFST -postfiles output/lyca.sin.fix.saf output/lyca.vic.fix.saf -priorfile output/lyca.2dsfs -outfile output/lyca.fst -nind 7 6 -nsites 99862 -block_size 50000 -islog 1
```

Look at the results:
```
less -S output/lyca.fst 
```
The 4th column is the per-site FST (negative values mean 0) while the 5th is the probability of being variable.
Can you see which sites are more likely to be SNPs?

You can plot these values using the R script `ngsTools/ngsPopGen/scripts/plotFST.R`.
It is convenient to filter out sites that are clearly not variable (using -t option, but try not to be strict), as these may add some noise.
Indeed, try to set this values to 0 (so consider all sites) and to 0.9 (very strict), and you should get a range a values that relates to the contribution of SNP calling.

As an additional note, we can estimate FST using a different 2D-SFS as a prior, by setting -maxlike 0, although this should be used only in case of very few sites (or to calculate a local SFS) and to avoid patchy spectra.
On a global scale and for very low coverage, FST using these 2 different spectra as priors may change quite a bit.

---------

... > SFS > SUMMARY STATISTICS

Similarly, we can also compute some basic statistics from these data, like (the expectation of) the number of segregating sites and fixed differences between populations, as briefly described [here](http://www.ncbi.nlm.nih.gov/pubmed/23979584) and [here](http://www.ncbi.nlm.nih.gov/pubmed/24260275).
This can be achieved by the following commands (first we need to incorporate 1D-SFS as prior in our .saf files):
```
ngsTools/angsd/misc/realSFS output/lyca.sin.fix.saf 14 -P 2 -nSites 1000000 2> /dev/null > output/lyca.sin.fix.sfs
ngsTools/angsd/misc/realSFS output/lyca.vic.fix.saf 12 -P 2 -nSites 1000000 2> /dev/null > output/lyca.vic.fix.sfs

ngsTools/angsd/angsd -b input/lyca/bams_sin.list -anc input/lyca/referenceseq.fasta -GL 1 -doSaf 1 -pest output/lyca.sin.fix.sfs -out output/lyca.sin.pp
ngsTools/angsd/angsd -b input/lyca/bams_vic.list -anc input/lyca/referenceseq.fasta -GL 1 -doSaf 1 -pest output/lyca.vic.fix.sfs -out output/lyca.vic.pp

ngsTools/ngsUtils/GetSubSfs -infile output/lyca.sin.pp.saf -posfile output/A.pos -nind 7 -nsites 99001 -len 98862 -outfile output/lyca.sin.pp.fix.saf
ngsTools/ngsUtils/GetSubSfs -infile output/lyca.vic.pp.saf -posfile output/B.pos -nind 6 -nsites 98907 -len 98862 -outfile output/lyca.vic.pp.fix.saf
```

Then you can run the program and plot the results (in sliding windows of 10kbp)
```
ngsTools/ngsPopGen/ngsStat -npop 2 -postfiles output/lyca.sin.pp.fix.saf output/lyca.vic.pp.fix.saf -nind 7 6 -nsites 98862 -iswin 1 -outfile output/lyca.stats -islog 1 -block_size 10000

## This is an example of output on windows of 10kbp
# head -n 20 output/lyca.stats
# 1     10000   88.460136       28.846587       64.962117       23.485323       4.566297
# 10001 20000   96.353437       31.153763       69.837429       25.690679       3.255709
# 20001 30000   76.649889       21.896089       65.718856       23.579390       0.336617
# 30001 40000   77.682144       22.883620       53.643992       18.386499       2.058954
# 40001 50000   99.323242       29.755223       72.716092       25.539969       0.658504

less -S output/lyca.stats
```

---------

... > SFS > PCA

Lastly, we can investigate the genetic structure of our samples via Principal Components Analysis (PCA).
ngsTools implements an estimation of the covariance matrix, as described [here](http://www.ncbi.nlm.nih.gov/pubmed/23979584).

First, we need to compute genotype posterior probabilities, as previously explained.
```
ngsTools/angsd/angsd -bam input/lyca/bams.list -ref input/lyca/referenceseq.fasta -GL 1 -doMaf 2 -doMajorMinor 1 -SNP_pval 0.05 -doGeno 32 -doPost 2 -out output/lyca -sites input/lyca/sites.angsd.bed
gunzip output/lyca.geno.gz
gunzip -c output/lyca.mafs | wc -l
#   1514 # -1 you get the nr of sites
```
Please note that here, in this scenario, we perform a (relaxed) SNP calling.

Now we can estimate the covariance matrix, decompose it, and plot the first 2 components:
```
# compute the covariance
ngsTools/ngsPopGen/ngsCovar -probfile output/lyca.geno -outfile output/lyca.covar -nind 20 -nsites 1513 -call 0 -norm 0

# prepare a suitable annotation file
awk -v OFS="\t" '$1=$1' input/lyca/lyca.clst > input/lyca/lyca.clst2

# plot
Rscript ngsTools/ngsPopGen/scripts/plotPCA.R -i output/lyca.covar -c 1-2 -a input/lyca/lyca.clst2 -o output/lyca.pca.pdf
open output/lyca.pca.pdf # on mac 
```



