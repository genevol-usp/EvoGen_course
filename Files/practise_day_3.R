
## These commands are meant to be run in R
## Make sure you have created a folder called Results in the current directory (if not do that with "mkdir Results")
## This script calls the program "ms", so be sure you have installed it
## Make sure you are in the Work/ folder.

# open R

# load all R functions we need
source("../Scripts/functions.R")

## 1) READ AND STORE OBSERVED GENETIC DATA (2D-SFS)

# the file "../Data/polar.brown.sfs" includes the joint (2D) site frequency spectrum (SFS) between polar bears (on the rows) and brown bears (on the columns)
# if you want to see this file type "cat ../Data/polar.brown.sfs" in your terminal

# read this file and store the the joint SFS of polar VS brown bears into a matrix
polar.brown.sfs<-as.matrix(read.table("../Data/polar.brown.sfs", stringsAsFactors=FALSE, header=FALSE))
# this matrix represents the joint (2D) site frequency spectrum (SFS)

# -----

# QUESTION: what is the sample size for this dataset?

# ANSWER: this is the unfolded spectrum which means that each population has 2n+1 entries in its spectrum, with n being the number of individuals
polar.nr_ind<-(nrow(polar.brown.sfs)-1)/2
brown.nr_ind<-(ncol(polar.brown.sfs)-1)/2
# on the other hand the number of chromosomes can be retrieved as
polar.nr_chrom<-nrow(polar.brown.sfs)-1
brown.nr_chrom<-ncol(polar.brown.sfs)-1

# -----

# QUESTION: how many sites (and SNPs) are there in this 2D-SFS?

# ANSWER: the number of analysed sites is simply the sum of all entries in the SFS
nr_sites<-sum(polar.brown.sfs)
# whereas the number of polymorphic sites is equal to the number of sites minus the count of sites with joint allele frequency (0,0) or (2*polar.nr_ind; 2*brown.nr_ind)
nr_snps<-as.numeric(nr_sites-polar.brown.sfs[1,1]-polar.brown.sfs[(2*polar.nr_ind+1), (2*brown.nr_ind+1)])

# -----

# for convenience, let's set to NA entries in the matrix where the sites is not a SNP
polar.brown.sfs[1,1]=NA
polar.brown.sfs[(2*polar.nr_ind+1), (2*brown.nr_ind+1)]=NA

# plot the spectrum
plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="2D-SFS")

# this is the unfolded spectrum BUT it has been generated using a reference allele to polarise the spectrum; this means that the polarisation is arbitrary and we need to use the folded spectrum instead
polar.brown.sfs<-fold2DSFS(polar.brown.sfs)

# plot the folded spectrum
plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="folded 2D-SFS")

# as an illustration, compute the FST value for the comparison polar vs brown
doFST(polar.brown.sfs)

# -----

## 2) SIMULATE GENETIC DATA UNDER DIFFERENT VALUES OF THE PARAMETER TO BE ESTIMATED

# define how many simulations we want to perform (ideally a lot)
nr_simul<-100

# define the prior distribution of our parameter to be estimated (divergence time)
# use a uniform prior bounded at realistic range values
tdiv_min<-300000 # 300k years ago
tdiv_max<-700000 # 700k years ago

# convert time in coalescent units
# we need generation time and the reference effective population size we used for writing our model in ms
# T(years)=T(coal)*4*Nref
gen_time <- 8.423
ref_pop_size <- 68000
tdiv_min_coal <- tdiv_min/(gen_time*ref_pop_size*4)
tdiv_max_coal <- tdiv_max/(gen_time*ref_pop_size*4)

# -----

# initialise output files
cat("", file="Results/divergence_distance.txt")

debug<-FALSE # set this to TRUE if you want to plot intermediate results

# set the directory wher you installed "ms" software
ms_dir<-"/home/BIO5789/utils/bin/ms" # this is my specific case, yours could be different

# iterate across all simulations
for (i in 1:nr_simul) # run this "for loop" nr_simul times
{

	# print
        cat("\n",i, "/", nr_simul, "\n")

	# pick a divergence time randomly from the prior distribution defined above
        tdiv_random_coal<-runif(1, min=tdiv_min_coal, max=tdiv_max_coal)
	# convert this value into years-units
	tdiv_random<-round(tdiv_random_coal*gen_time*4*ref_pop_size)

	# record the sampled value of divergence time in a file
	cat(tdiv_random, "\t", file="Results/divergence_distance.txt", append=TRUE)
	cat("\t", tdiv_random, "\n")

	# simulate genetic data of polar and brown bears under this sampled divergence time
	# use "ms" to simulate data

	# initalise ms output file (create an empty file)
	cat("", file="ms.txt")        

	# QUESTION: how many SNPs do we need to simulate?

	# ANSWER: this should match the observed value, stored in "nr_snps" 

	ms.command <- paste(ms_dir, "50", nr_snps, "-s 1 -I 2 36 14 -n 1 1 -n 2 6.8 -en 0.02269 1 0.07353 -en 0.05281 2 0.2941 -em 0.06459 1 2 4.896 -em 0.1392 1 2 0 -en 0.1392 1 0.2941 -ej", tdiv_random_coal, "2 1 -en 0.3924 1 1.809 > ms.txt")
	# run ms
	system(ms.command)

	# read ms output file and compute the 2D-SFS
	simulated.sfs<-fromMStoSFSwith1site("ms.txt", nr_snps, polar.nr_chrom, brown.nr_chrom)

	# fold it
	simulated.sfs<-fold2DSFS(simulated.sfs)
	simulated.sfs[1,1]=NA

	# if you want to plot it
	# plot2DSFS(simulated.sfs, xlab="Polar", ylab="Brown", main="Simulated 2D-SFS")
	# and calculate FST
	cat("\t", doFST(simulated.sfs), "\n")

	# compute the Eucledian distance between observed and simulated SFS (our summary statistics)
	# sqrt((OBS-SIM)^2)
	eucl_dist<-sum(sqrt((polar.brown.sfs-simulated.sfs)^2), na.rm=T)
	# store it in a file
	cat(eucl_dist, "\n", file="Results/divergence_distance.txt", append=T)

	cat("\t", eucl_dist, "\n") # print in std output

	# if you want to inspect the comparison
	if (debug) {

		matrix.dist<-sqrt((polar.brown.sfs-simulated.sfs)^2)
		# residuals
		matrix.resid<-(simulated.sfs-polar.brown.sfs)/(polar.brown.sfs+1)

		par(mfrow=c(2,2))
		plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="2D-SFS")
		plot2DSFS(simulated.sfs, xlab="Polar", ylab="Brown", main="Simulated 2D-SFS")
		plot2DSFS(matrix.dist, xlab="Polar", ylab="Brown", main="Eucledian distance")
		image(x=seq(0,polar.nr_chrom), y=seq(0,brown.nr_chrom), z=matrix.resid, xlab="Polar", ylab="Brown", main="Std. Residuals (SIM-OBS)/OBS")

	}

} 


## 3) ESTIMATION OF DIVERGENCE TIMES

# either stop your results or use my pre-computed file ("cp ../Data/divergence_distance.txt Results/.")

# read results from simulation
results <- read.table("Results/divergence_distance.txt", stringsAsFactors=F, head=F, sep="\t")
colnames(results) <- c("divergence","distance")
nr_sim <-nrow(results)

# check that we explored the whole prior range uniformily
hist(results$divergence, xlim=c(300000, 700000), xlab="Divergence time", main="Prior distribution")

# retain the 5% top simulations with the closest distance
rankings <- rank(results$distance, ties.method="min")/nr_sim
top <- which(rankings<0.05)

# plot the posterior distribution
hist(results$divergence[top], xlim=c(300000, 700000), xlab="Divergence time", main="Posterior distribution")
# explore the median and mode
summary(results$divergence[top])

# compute the 95th bayesian credible interval (BCI)
quantile(results$divergence[top], seq(0,1,0.025))[c(2,40)]

# compute bayes factors (posterior vs prior) in bins

prior_hist <- hist(results$divergence, breaks=seq(300000,700000,20000), plot=F)
prior_hist$counts <- prior_hist$counts/sum(prior_hist$counts)

post_hist <- hist(results$divergence[top], breaks=seq(300000,700000,20000), plot=F)
post_hist$counts <- post_hist$counts/sum(post_hist$counts)

bf_time <- post_hist$counts/prior_hist$counts
x_time <- post_hist$mid

plot(x=x_time, y=bf_time, main="Divergence time", frame=F, ty="b", xlab="years ago", xlim=c(300000,700000), ylab="Bayes Factor", pch=16, col="grey", cex=1.5, cex.axis=1.2, cex.lab=1.2)

# QUESTION:
# Are we getting a good fit to the observed data? Do we need to add more parameters or change the prior range?


## 4) SUGGESTED EXERCISES

# a) inspect the 2D-SFS under the estimated divergence time (simulate genetic data using "ms" assuming a divergence time equal to the median of the posterior probability)

# b) use other possible summary statistics (for instance FST) and vary the prior range

# c) use "abc" R package to compute distances and estimations (and specify the type of ABC algorithm to be applied, for instance "loclinear")


## 5) POSSIBLE SOLUTIONS

# a) 

tdiv.est <- 328800 / (gen_time*ref_pop_size*4)

cat("", file="ms.txt")
ms.command <- paste(ms_dir, "50", nr_snps, "-s 1 -I 2 36 14 -n 1 1 -n 2 6.8 -en 0.02269 1 0.07353 -en 0.05281 2 0.2941 -em 0.06459 1 2 4.896 -em 0.1392 1 2 0 -en 0.1392 1 0.2941 -ej", tdiv.est, "2 1 -en 0.3924 1 1.809 > ms.txt")
system(ms.command)

est.sfs<-fromMStoSFSwith1site("ms.txt", nr_snps, polar.nr_chrom, brown.nr_chrom)
est.sfs<-fold2DSFS(est.sfs)
est.sfs[1,1]=NA

doFST(est.sfs)

matrix.dist<-sqrt((polar.brown.sfs-est.sfs)^2)
matrix.resid<-(est.sfs-polar.brown.sfs)/(polar.brown.sfs+1)

par(mfrow=c(2,2))
plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="2D-SFS")
plot2DSFS(est.sfs, xlab="Polar", ylab="Brown", main="Simulated 2D-SFS under the estimated divergence time")
plot2DSFS(matrix.dist, xlab="Polar", ylab="Brown", main="Eucledian distance")
image(x=seq(0,polar.nr_chrom), y=seq(0,brown.nr_chrom), z=matrix.resid, xlab="Polar", ylab="Brown", main="Std. Residuals (SIM-OBS)/OBS")


## CREDITS
# Few lines has been written by S.D. Gopal.












