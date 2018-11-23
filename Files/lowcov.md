
In case of low-depth sequencing data, it is recommended to avoid calling SNPs and genotypes, since bad inferences at this level may inflate all downstream analyses.
With this data, summary statistics can be computed by keeping statistical uncertainty into account.
We will show some examples on how to compute summary statistics from low-depth sequencing data using ANGSD and ngsTools.

### ANGSD

We first need to compute the per-site sample allele frequency posterior probabilities (.saf files) and the SFS. 
Methods used here are described [here](http://www.biomedcentral.com/1471-2105/14/289).
We use the human dataset as an illustration.

Genotype likelihoods > SAF probabilities > SFS

```
OPT1='-minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1'
OPT2='-minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1'

ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -anc input/human/hg19_chr1.fa $OPT1 $OPT2 -GL 1 -doMaf 2 -doMajorMinor 5 -doSaf 1 -r 1: -out output/human

ngsTools/angsd/misc/realSFS output/human.saf 66 -P 2 2> /dev/null > output/human.sfs
```
We have already discussed how to estimate the SFS using ANGSD.

SAF probabilities + SFS > NUCLEOTIDE DIVERSITY

Now we can compute several indexes of nucleotide diversity per site (note that we give the estimated SFS as input):
```
ngsTools/angsd/angsd -bam input/human/bams.list -anc input/human/hg19_chr1.fa -ref input/human/hg19_chr1.fa $OPT1 $OPT2 -r 1: -out output/human -doThetas 1 -doSaf 1 -GL 1 -doMaf 2 -doMajorMinor 5 -pest output/human.sfs
```

Let us have a look at the output file:
```
gunzip -c output/human.thetas.gz | head

 Chromo	Pos	Watterson	Pairwise	thetaSingleton	thetaH	thetaL
 1	13999902	-8.118760	-8.174015	-7.906803	-8.704033	-8.404312
 1	13999903	-8.118760	-8.174015	-7.906803	-8.704033	-8.404312
 1	13999904	-8.118760	-8.174015	-7.906803	-8.704033	-8.404312
 1	13999905	-8.201048	-8.317305	-7.914413	-9.035470	-8.613257
 1	13999906	-8.201047	-8.317303	-7.914413	-9.035457	-8.613251
```
Values are per site and in log scale.

Often, one is interested in performing a sliding window analysis of some summary statistics, like Tajima D.
In ANGSD, this can be achieved using the following commands:
```
ngsTools/angsd/misc/thetaStat make_bed output/human.thetas.gz # for file indexing, optional
ngsTools/angsd/misc/thetaStat do_stat output/human.thetas.gz -nChr 66 -win 10000 -step 2000
```
Here, we computed statistics on windows of 10kbp and a step of 2kbp.

Results are stored here:
```
less -S output/human.thetas.gz.pestPG
```
and columns are:
 - (indexStart,indexStop)(posStart,posStop)(regStart,regStop) chrom window_center; <br>
 - 5 estimators of theta: Watterson, pairwise, Fu & Li, Fay H, L; <br>
 - 5 neutrality test statistics: Tajima D, Fu&Li F, Fu&Li D, Fay H, Zeng E. <br>
 - The final column is the effetive number of sites with data in the window. <br>

Therefore to get the values of Tajima D in these windows we can simply extract the corresponding columns
```
cut -f 2,3,9 output/human.thetas.gz.pestPG
```

**EXERCISE**
Perform a sliding window analysis on the butterfly dataset, by computing Tajima D and Fu & Li metrics.
Use the reference sequence to polarize alleles.
Plot the results. 
A possible **SOLUTION** for this exercise is given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/solutions.txt).
Alternatively, use your own data set.




