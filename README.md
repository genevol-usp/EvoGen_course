# EvoGen_course

Evolutionary Genomics course - Data Analysis module

13th-17th April 2015

## Material

To download all the material in this web page use [git](http://git-scm.com/):

	git clone https://github.com/mfumagalli/EvoGen_course

Enter the folder and create a subfolder "Work" where you will be running the exercises:

	cd EvoGen_course
	mkdir Work
	cd Work
	mkdir Results # you results will be saved here

Instructions to download remaining input files and to obatin necessary programs for the practical sessions are given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/install.md).

Feel free to contact me (mfumagalli82@gmail.com) if you have problems.

## Agenda

### Day 1

 *	Lecture: From raw NGS data to genotypes ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides/))

 *	Paper discussion: ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides/)) and ([paper](https://github.com/mfumagalli/EvoGen_course/tree/master/Papers/))

 *	Practical session: 
	+	[data filtering](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/filtering.md)
	+       [genotype calling](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/genocall.md)

 *	Additional material:

	+       genotype calling using:
		-	[SAMtools](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/genocall_samtools.md)
		-	[BEAGLE](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/imputation.md)

### Day 2

 *	Lecture: SNP calling and advanced methods for evolutionary inferences from NGS data ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides))

 *	Paper discussion: ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides)) and ([papers](https://github.com/mfumagalli/EvoGen_course/tree/master/Papers/))

 *	Practical session: 

	+       [SNP calling](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/snpcall.md)
	+	advanced methods to estimate [SFS](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/sfs.md)

 *	Additional material: 

	+       SNP calling using [SAMtools](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/snpcall_samtools.md)
	+	estimation of [inbreeding](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/inbreeding.md) coefficients
	+	advanced methods to calculate summary statistics using:
		-	[ANGSD](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/lowcov.md)
		-	[ngsTools](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/lowcov_ngstools.md)
		-	[ngsTools/ngsDist](https://github.com/mfumagalli/ngsTools/tree/master/TUTORIAL.md)


### Day 3

 *	Lecture: Population structure and demographic inferences ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides/))

 *	Paper discussion: ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides)) and ([paper](https://github.com/mfumagalli/EvoGen_course/tree/master/Papers))

 *	Practical exercises:

	+	estimating demographic parameters through simulations ([R script](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/practise_day_3.R))
	+	inference of population splits and mixture events using TreeMix ([bash script](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/practise_day_3_extra.txt))

### Day 4

 *	Lecture: Detecting natural selection ([slides](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides))

 *	Paper discussion: slides ([pdf](https://github.com/mfumagalli/EvoGen_course/tree/master/Slides)) and ([paper](https://github.com/mfumagalli/EvoGen_course/tree/master/Papers))

 *	Practical exercises:

	+	selection scan in the human genome and assessment of significance ([R script](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/practise_day_4.R))
	+	summary statistics using [ngsTools/ngsDist](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)

### Day 5

Evaluation.



## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en).
Sequencing data on butterfiles has been kindly provided by [Zach Gompert](https://gompertlab.wordpress.com/).
ANGSD has been developed by Thorfinn Korneliussen and Anders Albrechtsen. 
Filipe G. Vieira implemented the inbreeding estimation. 
Tyler Linderoth wrote most of data filtering scripts together with [Sonal Singhai](https://systemsbiology.columbia.edu/people/sonal-singhal), [Ke Bi](http://scholar.google.ca/citations?user=ymcwERQAAAAJ), and Filipe G. Vieira.



