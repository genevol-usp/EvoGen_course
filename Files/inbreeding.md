
Please note that all command lines below may not work with more recent versions of ANGSD, but are provided as a general reference.
For more information on how to call genotypes in case of inbred species please refer to [ngsF](https://github.com/fgvieira/ngsF) web page.

### Inbreeding

When studying domesticated or self-pollinated species it is extremely important to estimate inbreeding coefficients and include such estimates into genotype and SNP calling.
ANGSD and ngsTools are suitable for this purpose as they have built-in functionalities to compute per-site and per-individual inbreeding coefficients, and incorporate such estimates when calling SNPs and genotypes and calculating the SFS.
Notably, methods implemented are suitable for low-coverage data as again they do not rely on called genotypes.

As an illustration, we will apply these methods to simulated data which will incorporate a significant level of inbreeding.
We will use ngsTools to simulate such data. 
Please note that more sophisticated programs to simulated NGS data are available, e.g. [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/).
ngsTools is useful when simulating multiple populations with a specified degree of genetic differentiation, variable coverage among samples, and deviation from HWE.

First, let us simulate 10 individuals with a mean per-sample inbreeding coefficient of 0.5.
We simulated 1000 sites variable (in the population, not in the sample), with 6X as mean coverage
```
ngsTools/ngsSim/ngsSim -npop 1 -nsites 1000 -pvar 1 -depth 6 -errate 0.01 -F 0.5 -seed 1 -outfiles output/simul
```

There are a lot of files in output. For a complete description read [here](https://github.com/mfumagalli/ngsSim/tree/master/examples).
We are interested in the `.glf.gz` files which are the per-site per-sample genotype likelihoods.
We first need to convert these values into an accetable format for ngsF (the program in ngsTools dedicated to the calculation of inbreeding coefficients).
```
ngsTools/angsd/angsd -sim1 output/simul.glf.gz -nInd 10 -doGlf 3 -doMaf 2 -SNP_pval 0.001 -out output/simul.geno -doMajorMinor 1
```

Now we can estimate inbreeding coefficients (fast way)
```
N_SITES=$((`gunzip -c output/simul.geno.mafs.gz | wc -l`-1))
gunzip -c output/simul.geno.glf.gz | ./ngsTools/ngsF/ngsF -n_ind 10 -n_sites $N_SITES -glf - -verbose 2 -min_epsilon 0.01 -chunk_size 100 -approx_EM -init_values r -max_iters 100 -out output/simul.indF
```
Please note that this may not work on all machines, especially of mac.
If this is the case please type:
```
cp input/simul.indF* output/.
```

We can finally incorporate these estimate inbreeding coefficients in our calculation of allele frequencies, genotypes and summary statistics.
```
ngsTools/angsd/angsd -sim1 output/simul.glf.gz -nInd 10 -doGeno 3 -doPost 1 -doMaf 2 -out output/simul.indF -doMajorMinor 1 -indF output/simul.indF
gunzip -c output/simul.mafs.gz | head
```







