
The instructions below may not be fully operating using the latest version of SAMtools.

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES

Here we will use a custom Perl script to filter sites, although there are many available tools for that.
The scope here is to show the main steps for data filtering of sites, rather than giving a fixed pipeline.
Our aim here is also to illustrate how filtering indeed changes the data and how important it is to keep track of such changes.

-----------------------

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES (mpileup , VCF/BCF)

Main quality control filters can be grouped as:

1. Depth
* Minimum site read depth <br>
* Maximum site read depth <br>
* Even depth across individuals <br>
* Minimum number of alternate alleles <br>

2. Bias and other quality aspects
* Minimum RMS mapping quality for SNPs <br>
* Mapping quality bias <br>
* Strand bias <br>
* Allele bias in potential heterozygotes <br>
* Base quality bias <br>
* Distance from end of read bias <br>

3. Hardy-Weinberg Equilibrium
* Exact test of HWE <br>
* Excess of heterozygotes <br>
* Excess of homozygotes <br>

4. Biallelic site filter

5. Mutation type removal filters

Read depth cutoffs can be chosen by first investigating the empirical distribution.
Again, we may want to remove sites with excessively high or low depth compared to the distribution.

Strand bias measures the imbalance for the number of reads at each strand. The end read bias test looks for unusually high/low number of alternate alleles at the of the reads, which could be due to bad trimming.
Allele bias in potential heterozygotes tests for potential barcode/allele swapping.

Let us see the options for this script (implemented by T. Linderoth and F.G. Vieira and others).
This script reads a VCF file and filter SNPs based on a set of rules.
It is similar to SAMtools `varfilter.pl` script but it extends some of the filters.
Please note that this script is only for SNP filtering and will ignore INDELs (which we have already filtered out using SAMtools).

```
perl scripts/SNPcleaner.pl --help

Usage:
SNPcleaner.pl [OPTIONS] <infile.vcf>
or
cat <infile.vcf> | SNPcleaner.pl [OPTIONS]
############ OPTIONS ############

input files:

--pop_info|-P   FILE    input file with population information (each line should list `sample_name      pop_ID`) for HWE filtering
--anc|-A        FILE    ancestral-state fasta file (with FAI in same directory)
--pileup|-G     FILE    pileup format file (all sites and individuals in input vcf must be in pileup)
--exons|-X      FILE    BED format file of exonic regions (sorted from lowest to highest numbered contig)

coverage filters:

--minDepth|-d   INT     minimum site read depth [2]
--maxDepth|-D   INT     maximum site read depth [1000000]
--minIndiv|-k   INT     minimum number of individuals with at least [-u INT]X coverage (requires SNPcleaner -u and mpileup -D) [1]
--minIndiv_cov|-u       INT     minimum individual coverage threshold used for -k (requires SNPcleaner -k and mpileup -D) [0]
--minalt|-a     INT     minimum number of alternate alleles per site [0]

bias and other quality-aspect filters:

--minRMSmap|-Q  INT     minimum RMS mapping quality for SNPs [10]
--mapqual|-f    FLOAT   min p-value for map quality bias [0]
--strand_ind|-S FLOAT   min p-value for strand bias from combining p-values across individuals [0.0001]
--strand_site|-s        FLOAT   min p-value for strand bias determined from read counts summed over all individuals at the site [0.0001]
--allele_bias|-T        FLOAT   min p-value for allele bias in potential heterozygotes [0.001]
--hetero_llh|-R FLOAT   skip called homozygotes with heterozygote likelihood less than -R FLOAT for allele bias filter [0.05]
--basequal|-b   FLOAT   min p-value for base quality bias [1e-100]
--endbias|-e    FLOAT   min p-value for biased distance of alternate bases from ends of reads (indication of misalignment) [0.0001]

Hardy-Weinberg equilibrium filters:

--hwe|-h                FLOAT   min p-value for exact test of HWE (two-tailed) [0]
--hetero_excess|-H      FLOAT   min p-value for exact test of excess heterozygotes [0]
--hetero_deficit|-L     FLOAT   min p-value for exact test of deficient heterozygotes [0]
--inbreed_coef|-F       FLOAT   inbreeding coefficient value [0]
--rmv_nonHWEexons|-g    filter-out exons containing at least one SNP out of HWE (requires -X)

mutation type filters:

--mutation_rmv|-M       STR     mutation type(s) to remove (ex: `-M CT_GA` means remove C<=>T and G<=>A)
--one_dir|-w    remove mutations defined by -M in one direction (ex: `-M CT_GA` means remove C=>T and G=>A) (requires -A)
--alt_excess|-J FLOAT   min p-value for excess substitutions defined by -M in sites called as nonvariable (requires -A if -w) [1e-06]
--error|-E      FLOAT   sequencing error rate (required by -J for filtering mutation types) [0.01]

general filters:

--nonvar|-v     process nonvariant sites (in addition to variants)
--keep_nonbinary|-2     keep non-biallelic sites
--exclude_contigs|-r    FILE    list of contigs/chromosomes to exclude (each line should list a contig name exactly as it appears in the input vcf file)
--rmv_nonexonic|-t      filter-out non-exonic sites (requires -X)

output:

--bed|-B        FILE    name of dumped BED format file for sites that pass all filters
--failed_sites|-p       FILE    name of dumped file containing sites that failed at least one filter (bziped)
--vcfout|-o     FILE    name of dumped vcf file containing sites that passed filters
--ind_depth|-I  FILE    dumped file with individual mean depth

##########

Notes:

Some of the filters rely on annotations generated by SAMtools/BCFtools.
To use the eveness-of-coverage filters (options -k and -u), -D must be used with satmools mpileup.
If option -s or -S is set to 0, that particular strand bias filter is not performed.
It is recomended to use mpileup -I to ignore indels.
Characters in front of filtered sites (dumped with option -p) indicate filters that the site failed to pass.

```
So there are a lot of options.
Input file is a VCF.
Some of these options require the mpileup file, which we have already produced.
Moreover, we may need to give an ancestral sequence (indexed, `input/lyca/referenceseq.fasta`) for some mutation-type filters (e.g. for ancient DNA).


-----

**WORKFLOW**:
... > FILTERING SITES > DEPTH

First we want to remove sites with extremely low or high depth.
Let us get the empirical distribution.

```
cat output/lyca.raw.vcf | perl scripts/SNPcleaner.pl -v -a 0 -k 10 -u 2 -d 1 -S 0 -s 0 -G output/lyca.mpileup -M NN_NN -o output/lyca.filt.vcf 2> /dev/null
# it should print how many sites were processed and the individual mean depths
```

Let us get the distribution of per-site depth values (script by F.G Vieira):
```
cat output/lyca.filt.vcf | tr ";" "\t" | awk 'BEGIN{print "chr_pos\tdepth"} !/#/{sub("DP=","",$8); if(rand()<=1) print $1"_"$2"\t"$8}' > output/lyca.sdepth
Rscript scripts/plotDepth.R output/lyca.sdepth output/lyca.sdepth.pdf
open output/lyca.sdepth.pdf # on mac
```
You can investigate the distribution of per-site depth:
```
     0%    2.5%      5%    7.5%     10%   12.5%     15%   17.5%     20%   22.5%
  42.00   43.00   53.00   54.00   61.00   63.00   63.00   72.00   77.00   84.00
    25%   27.5%     30%   32.5%     35%   37.5%     40%   42.5%     45%   47.5%
 116.00  151.00  154.00  157.00  163.00  176.00  178.00  180.00  183.00  183.00
    50%   52.5%     55%   57.5%     60%   62.5%     65%   67.5%     70%   72.5%
 191.00  196.00  236.65  250.00  274.00  303.00  313.00  327.00  334.00  371.00
    75%   77.5%     80%   82.5%     85%   87.5%     90%   92.5%     95%   97.5%
 397.00  399.00  571.00  655.00  667.00  716.25  748.00  873.00 1067.00 1131.00
   100%
1576.00
```
and decide which cutoffs are more suitable for your scope.


---------------

**WORKFLOW**:
... > FILTERING SITES > BIAS

We are now ready to perform the first round of filtering.
Here some parameters:

Parameter | Value | Description
--------- | ----- | -----------
-G | output/lyca.mpileup | mpileup file <br>
-P | input/lyca/clst.txt | this gives the correspondence sample-pop for HWE filtering <br>
-p | output/lyca.bad.bz2 | list of filtered sites <br>
-d | 50 | minimum site depth <br>
-D | 900 | maximu site depth <br>
-k,-u | 10,2 | at least 10 samples must have 2 reads <br>
-a | 0 | minimum number of alternate alleles per site <br>
-Q | 20 | minimum RMS mapping quality for SNPs <br>
-f | 0.001 | min p-value for map quality bias <br>
-S | 0.001 | min p-value for strand bias from combining p-values across individuals <br>
-T | 0.001 | min p-value for allele bias in potential heterozygotes <br>
-b | 0.001 | min p-value for base quality bias <br>
-e | 0.001 | min p-value for biased distance of alternate bases from ends of reads <br>
-h | 0.001 | min p-value for exact test of HWE <br>
-v | process nonvariant sites <br>

Let us write down the command:
```
OPT1='-a 0 -k 10 -u 2 -d 50 -D 900 -Q 20'
OPT2='-f 0.001 -S 0.001 -s 0 -T 0.001 -b 0.001 -e 0.001 -h 0.001 -v'
cat output/lyca.raw.vcf | perl scripts/SNPcleaner.pl $OPT1 $OPT2 -G output/lyca.mpileup -P input/lyca/clst.txt -M NN_NN -o output/lyca.filt.vcf -p output/lyca.bad.bz2
```
How many sites did we filter? How many sites were retained?
We can inspect the removed sites to better tune our filtering parameters.
```
bunzip2 output/lyca.bad.bz2
less -S output/lyca.bad
```

Here some examples of filtered sites:

* low depth:
```
kd      ref_contig      50160   .       A       .       23.7    .       DP=1;;AC1=25;FQ=-23.7   PL:DP:SP        0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0
```

* N in the reference
```
Nkd     ref_contig      81199   .       N       .       5.42    .       DP=1;AF1=1;AC1=40;DP4=0,0,1,0;MQ=50;FQ=-23.7    PL:DP:SP        0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   26:1:0  0:0:0   0:0:0   0:0:0   0:0:0   0:0:0   0:0:0
```

* HWE test:
```
kdh(p=0.000108018128040375;0)   ref_contig      51494   .       C       A       8.93    .       \
DP=14;VDB=5.960000e-02;RPB=8.745357e-01;AF1=0.6656;AC1=27;DP4=1,0,2,0;MQ=37;FQ=3.43;PV4=1,0.14,0.11,1   GT:PL:DP:SP:GQ  0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:20,3,0:1:0:4        0/1:26,3,0:1:0:4        0/0:0,3,34:1:0:4        0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3 0/1:0,0,0:0:0:3
```
Indeed here all individuals are heterozygotes.

As you can see we have many sites that failed `S` test (strand bias) so maybe we have been too strict.
For other tests, it seems that all sites passed so maybe we have been too relaxed.
Let us check the proportion of sites filetered out by each test.
```
Rscript ../scripts/plotFilters.R output/lyca.bad output/lyca.bad.pdf 7000
open output/lyca.bad.pdf # on mac
```
We can rerun it by adjusting our parameters from the considerations.
```
OPT1='-a 0 -k 7 -u 1 -d 20 -D 900 -Q 20'
OPT2='-f 0.001 -S 1e-100 -s 0 -T 0.05 -b 0.05 -e 0.05 -h 0.05 -v'
cat output/lyca.raw.vcf | perl scripts/SNPcleaner.pl $OPT1 $OPT2 -B output/lyca.good.bed -G output/lyca.mpileup -P input/lyca/clst.txt -M NN_NN -o output/lyca.filt.vcf -p output/lyca.bad2.bz2
bunzip2 output/lyca.bad2.bz2
```

Let us check sites that failed due to quality bias:

* base quality bias
```
b       ref_contig      50209   .       C       .       24.8    .       DP=39;RPB=1.633743e+00;AF1=0.02792;AC1=1;DP4=34,0,1,0;MQ=44;FQ=-24.7;PV4=1,0.024,0.031,1        PL:DP:SP        0:2:0   0:1:0   0:5:0   0:0:0   0:2:0   0:4:0   0:1:0   0:1:0   0:2:0   0:1:0   0:2:0   0:1:0   0:2:0   15:4:0  0:1:0   0:0:0   0:1:0   0:0:0   0:1:0   0:4:0
grep 50209 output/lyca.mpileup
ref_contig      50209   C       2       ..      GH      1       .       H       5       .....   HGIIE   0       *       *       2       ..      GI      4       ....    IHHI    1       ..      HG      1       .       I       2       ..      6G      1       .       I       2       ..      GG      4       ...T    G5D<    1       .       D       0       *       *       ....    EHHI
```

* map quality bias
```
f       ref_contig      50222   .       T       .       8.03    .       DP=39;RPB=1.528942e+00;AF1=0.04621;AC1=1;DP4=33,0,1,0;MQ=44;FQ=-7.86;PV4=1,1,3.9e-05,1  PL:DP:SP        0:2:0   0:1:0   0:5:0   0:0:0   0:2:0   0:5:0   0:1:0   0:0:0   0:2:0   0:1:0   0:2:0   0:1:0   0:1:0   0:4:0   0:1:0   0:0:0   36:1:0  0:0:0   0:1:0   0:4:0
grep 50222 output/lyca.mpileup
ref_contig      50222   T       2       ..      GI      1       .       G       5       .....   HIIEH   0       *       *       2       ..      >I      5       .....   IEHGI   1       ..      G:      1       .       G       2       ..      HG      1       .       F       1       .       I       4       ....    >H@E    1       .       G       0       *       *       ....    EDFI
```

* strand bias
```
DS      ref_contig      55250   .       C       .       82.8    .       DP=567;AF1=0;AC1=0;DP4=551,0,0,0;MQ=46;FQ=-82.7 PL:DP:SP        0:31:0  0:36:0  0:34:0  0:22:0  0:42:0  0:27:0  0:33:0  0:33:0  0:15:0  0:18:0  0:28:0  0:22:0  0:27:0  0:14:0  0:18:0  0:21:0  0:33:0  0:25:0  0:37:0  0:35:0
mac:try matteo$ grep 55250 output/lyca.mpileup
ref_contig      55250   C       31      ............................... IGI7EIHHDI@E?GBIGHIIIIIIIDIHHGI 36      ....................................    HGGIEHEAHIEIGHGIBH>GIHFGGBH=HHGDIHHH    34      ..................................      IIIHIHHIID=GHIIHI@GICI;EID;HAIIHHE      22      ......................  BHDGIDIIHBDIGHGHGGHH?C  42      ..........................................      @=IHHIHIHIIDHIIIDDHHIGIGDIHHIHGII7IGDEI@II      27      ...........................     IEIHDII5BHIIHIIDHHHGIAHIIIH     33      .................................       <IGIIIIEHGGHIE>IHIHGGGGIDIGIGIHHD       33      .................................       HBIHIIGIIIIIIHIHGH?III?>BIH@D=IIH       15      ............... GHGHIIGIIGIHBGG 18      ..................      >HHIGGDEHHIIHIHGII      28      ............................    HGHGEHIHIDHIHIGIIIIH?III=HII    22      ......................  IIGIGDBGIHBIIHIIBIIGFI  27      ...........................     8II<?DHIGDIHHGIHFIHIIHADHIH     14      ..............  IFGHDGDGHHHHDH  18      ..................      IIDBIDGBGCIH>IDGII21    .....................   IGFHIH<IHIEIIGIDIHGIH   33      .................................       HGIBIIGFHIGHHIHDIDHHDIHHIHDHIHHH;       25      .........................       DIBHIIBDEDDH:BADIGIBDH:II       37      .....................................   IHIIGHHIGGIGHIAHHHGGIIHIHIBIH7IGHHFDI   35      ...................................     ?FIIBIEIIIHG?IIIHHIIGIAEIGHDDHIHHIH
```

-------------

-------------

**WORKFLOW**:
... > FILTERING SITES > CHECK

Now we check the resulting Site Frequency Spectrum (SFS).
We will see later how to accurately estimate it from sequencing reads.
So far, we will quickly compute it from the resulting VCF file.

```
perl -ne 'print "$1\n" if /AF1=([^,;]+)/' output/lyca.filt.vcf > output/lyca.vcf.freqs
Rscript scripts/plotSFS.R output/lyca.vcf.freqs output/lyca.vcf.sfs.pdf
open output/lyca.vcf.sfs.pdf # on mac
```
You can also get it [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/lyca.vcf.afs.pdf).

The output for this filtering step is a BED file (along with the filtered VCF file) with the positions that passed the quality filters.
```
head output/lyca.good.bed
```
This file can be used for all downstream analyses, along with the final VCF with only the valid sites.



