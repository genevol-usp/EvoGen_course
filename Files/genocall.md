
**WORKFLOW**:
FILTERED DATA -> GENOTYPE CALLING

Here we will explore several ways to call genotypes from sequencing data.
We will use a subset of human sequencing data from the [1000 Genomes Project](http://www.1000genomes.org/). 
This dataset comprises 33 individuals of European descent.
Our goal here is to call genotypes for each individual at each site.
We will use ANGSD (and SAMtools as additional material), and compare results using different options.

As a preliminary step (if not already provided) we need to create a text file with a list of BAM files, and created indexes for such BAM files.
You should have already run this, but in case here is the command line:
```
gunzip -c input/human/hg19_chr1.fa.gz > input/human/hg19_chr1.fa
ls input/human/*bam > input/human/bams.list
for i in input/human/smallNA*.bam; do samtools-1.2/samtools index $i; done
samtools-1.2/samtools faidx input/human/hg19_chr1.fa
```

### ANGSD

[ANGSD](http://popgen.dk/angsd/index.php/Main_Page#Overview) is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes.
ANGSD is not part of ngsTools, but the latter automaticlaly links to the latest version of ANGSD.

First, let us see all things ANGSD can do:
```
ngsTools/angsd/angsd

	-> angsd version: 0.615	 build(Apr  2 2015 09:56:30)
	-> Please use the website "http://www.popgen.dk/angsd" as reference
	-> Use -nThreads or -P for number of threads allocated to the program

Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-doHWE		Est inbreeding per site
	-doGeno		Call genotypes
	-doFasta	Generate a fasta for a BAM file
	-doAbbababa	Perform an ABBA-BABA test
	-sites		Analyse specific sites (can force major/minor)
	-doSaf		Estimate the SFS and/or neutrality tests genotype calling
	-doHetPlas	Estimate hetplasmy by calculating a pooled haploid frequency

	Below are options that can be usefull
	-bam		Options relating to bam reading
	-doMajorMinor	Infer the major/minor using different approaches
	-ref/-anc	Read reference or ancestral genome
	many others

For information of specific options type: 
	./angsd METHODNAME eg 
		./angsd -GL
		./angsd -doMaf
		./angsd -doAsso etc
		./angsd sites for information about indexing -sites files
Examples:
	Estimate MAF for bam files in 'list'
		'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```

We now see how to use ANGSD to call genotypes.
The specific option is `-doGeno`:
```
ngsTools/angsd/angsd -doGeno
	...
	-doGeno	0
	1: write major and minor
	2: write the called genotype encoded as -1,0,1,2, -1=not called otherwise counts of derived
	4: write the called genotype directly: eg AA,AC etc 
	8: write the posterior probability of all possible genotypes
	16: write the posterior probability of called genotype
	32: write the posterior probability of called genotype as binary
	-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
	-postCutoff=0.333333 (Only genotype to missing if below this threshold)
	-geno_minDepth=-1	(-1 indicates no cutof)
	-geno_maxDepth=-1	(-1 indicates no cutof)
	...
```

We also want to first perform a SNP calling so we assign genotypes only for putatively variable sites.
We may also want to explore some basic filtering that can be done by ANGSD, as explained [here](http://popgen.dk/angsd/index.php/Filters).
Briefly, these are the options for filtering we may want to use:

Parameter | Meaning
--------- | -------
-minInd 10 | use only sites with data from at least N individuals
-setMinDepth 20 | minimum total depth
-setMaxDepth 400 | minimum total depth

Other options can be see with:
```
ngsTools/angsd/angsd -bam
...
-remove_bads	1	Discard 'bad' reads, (flag >=255) 
-uniqueOnly	0	Discards reads that does not map uniquely
-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
-minMapQ	0	Discard reads with mapping quality below
-minQ		13	Discard bases with base quality below
-trim		0	Number of based to discard at both ends of the reads
-only_proper_pairs	1	Only use reads where the mate could be mapped
-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
-baq		0	adjust qscores around indels (as SAMtools), supply -ref
```

Therefore, our command can be:
```
ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1 -minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01  -r 1: -doGeno 2 -doCounts 1 -doPost 1 -out output/human
```
Here we used the reference allele as major, and analysing only the first chromosome (data is a subset anyway).

Here, we also used another option, `doPost`:
```
ngsTools/angsd/angsd -doPost
	...
	-doPost	0	(Calculate posterior prob ...)
	1: Using frequency as prior
	2: Using uniform prior
	...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium. 
We will see later what to do when the assumption of HWE is not valid.

Let us look at the output files:
```
# output/human.arg	output/human.geno.gz	output/human.mafs.gz
gunzip -c output/human.geno.gz | wc -l
# 138
gunzip -c output/human.geno.gz | less -S
```
Genotypes are coded as 0,1,2, as the number of alternative alleles. 

If we want to print the major and minor allele then we set `-doGeno 3`:
```
ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1 -minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01  -r 1: -doGeno 3 -doCounts 1 -doPost 1 -out output/human
gunzip -c output/human.geno.gz | less -S
```

Values coded as -1 represent missing data. Indeed, genotypes with a posterior probability less than a specified threshold are set as missing. 
You can vary this cutoff by setting the option -postCutoff.
For instance, if we set `-postCutoff 0`:
```
ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1 -minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -r 1: -doGeno 2 -doPost 2 -out output/human -postCutoff 0
```
we get
```
gunzip -c output/human.geno.gz | grep -1 - | wc -l
#       0
```
sites with missing genotypes, while if we are more stringent and use 0.95 as cutoff
```
ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1 -minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1 -GL 1 -doMaf 2 -doMajorMinor 4 -SNP_pval 0.01 -r 1: -doGeno 2 -doPost 2 -out output/human -postCutoff 0.95
```
we get
```
gunzip -c output/human.geno.gz | grep -1 - | wc -l
```
sites with at least one individual with missing genotype (all!).
Why is that?

Let us look at the depth distribution for this dataset, by using `-doCounts 1 -doDepth 1`:
```
ngsTools/angsd/angsd -bam input/human/bams.list -ref input/human/hg19_chr1.fa -minInd 10 -setMinDepth 20 -remove_bads 1 -uniqueOnly 1 -minMapQ 10 -minQ 10 -trim 5 -only_proper_pairs 1 -C 50 -baq 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -r 1: -out output/human
# script to plot distribution of depth
Rscript scripts/plotDepthANGSD.R output/human.depthGlobal output/human.depthSample output/human.depth.pdf
open output/human.depth.pdf # on mac
```
Indeed, the mean depth per sample is around 4, therefore genotypes cannot be assigned with very high confidence.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application. 
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity). 
We will show later how to accurately estimate summary statistics with low-depth data.


### Inbred species

**ADDITIONAL MATERIAL**
For some studies (domesticated or self-pollinated species), it is important to consider any deviation from HWE when performing genotype or SNP calling.
We provide some command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/inbreeding.md)) to achieve this using ngsTools.


### SAMtools

**ADDITIONAL MATERIAL**
We also provide command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/genocall_samtools.md)) to call genotypes using SAMtools, and to compare results with ANGSD.


### BEAGLE

**ADDITIONAL MATERIAL**
You can also use BEAGLE to increase accuracy of your genotype calling.
Several examples using BEAGLE to impute data are given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/imputation.md).

### FreeBayes

Freebayes is another tool for SNP and Genotype calling, available [here](https://github.com/ekg/freebayes).
It is especially suitable for indels and CNVs detection.

### GATK

Alternatively, one can use GATK, which runs slower and requires more steps. Here are commands to generate a VCF file from the previous examples:
```
# GATK does not support .gz compressed references
gunzip ref/hg19.fa.gz
# GATK requires a dictionary for the reference
# remember to change the path to the program
java -jar CreateSequenceDictionary.jar R= ref/hg19.fa O=ref/hg19.dict
# run gatk
./gatk  -R ref/hg19.fa -T UnifiedGenotyper -I bams.list -L 1 -nt 4 -o gatk.vcf
# zip the file again for safety
gzip ref/hg19.fa
```

Look at the generated VCF file (gatk.vcf). How many SNPs does the program predict?











