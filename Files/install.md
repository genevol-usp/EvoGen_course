
## Requirements

Please note that this tutorial has been designed and tested on unix and mac.

Mac users need the developer tools, as [here](https://developer.apple.com/xcode/) or installed using [homebrew](http://brew.sh/).

If you are using Windows please consider using a virtual machine.
Some quick guidelines can be found [here](http://www.howtogeek.com/170870/5-ways-to-run-linux-software-on-windows).
You can use VMware Player and VirtualBox, and then you can install a light linux distribution like Xubuntu or Lubuntu.
Alternatively, you can use 'cygwin' to recreate a command line environment like in unix/mac.

You should have installed at least SAMtools, FastQC, ANGSD, ngsTools (see details below) as these will be the main programs used.
No knowledge of command line scripting is required.

## Software

You need the following software to run most of these examples:

 - For days 1 and 2:

 	+ [SAMtools](http://samtools.sourceforge.net/), [bcftools](https://samtools.github.io/bcftools/bcftools.html) and [VCFtools](http://vcftools.sourceforge.net/)

 	+ [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

 	+ [Sickle](https://github.com/najoshi/sickle) and [Scythe](https://github.com/vsbuffalo/scythe)

 	+ [ANGSD](http://popgen.dk/wiki/index.php/ANGSD)

 	+ [ngsTools](https://github.com/mfumagalli/ngsTools)

 - For day 3:

	+ [ms](http://home.uchicago.edu/rhudson1/source.html)

	+ [TreeMix](https://code.google.com/p/treemix/) [optional]

 - For day 4:

	+ [ngsDist](https://github.com/fgvieira/ngsDist) [optional]

For parsing and plotting you will also need [R](http://www.r-project.org/) and 'Rscript' installed.
You will also need some R packages to run all examples, and these are: 'ggplot2', 'VennDiagram', 'optparse', 'fields'.
[Perl](http://www.perl.org/) and [Python](https://www.python.org/) may also be used.

Please read individual documentation for proper installation of each software. As a general reference, we provide below instructions on how to download and install main programs used.

Additional instructions to install optional programs can be found [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/install_additional.md) but it is very likely that you will not need these during this course, so feel free to ignore it. It is provided for the sake of completeness.

## Data

We will use some example data sets:

 - Raw reads data of **chipmunks** are from [this](http://www.ncbi.nlm.nih.gov/pubmed/24118668) paper, and are kindly provided by Tyler Linderoth and Ke Bi.
In this paper, Authors sequenced around 4Mbp from early 20th century samples. We will use a very small subset of the original dataset.

 - BAM files of **butterflies** are a small sample subset of the original study and we will analyse around 1 Mbp.
More information about the original aligned files can be obtained at the original paper [here](http://www.ncbi.nlm.nih.gov/pubmed/22759293).
In this study, authors analysed a butterfly Next-Generation Sequencing dataset sequenced using Illumina GAII technology.
The dataset consists of a total of 381 samples of Lycaeides idas, Lycaeides melissa and an hybrid collection (Jackson Hole).
DNA resequencing was conducted on custom reduced genomic complexity libraries (RAD-seq).
We will analysed only a subset of 20 samples.

 - Other BAM files that will be used are a subset of **human** sequencing data from the 1000 Genomes Project.
More information on this project can be found [here](http://www.1000genomes.org/).
This dataset comprises 33 individuals of European descent.

 - Genotype likelihoods from inbred samples (optional) will be generated on-the-run.

 - Genetic data for demographic and selection inferences are provided in the [Data](https://github.com/mfumagalli/EvoGen_course/tree/master/Data) folder.

All these files can be downloaded [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/input.tar.gz).
User is `student` and password is `ANU2014`.

Untar these files by typing:
```
tar -xvf input.tar.gz
rm input.tar.gz
```
This file is rather big (~550M)
You can use a smaller version (~100M) if you prefer, but you will not be able to run examples on human dataset.
If this is the case, you can download it [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/input_light.tar.gz).
```
tar -xvf input_light.tar.gz
rm input_light.tar.gz
mv input_light/ input/ # if you run this you will not be able to use human examples
```

As a backup, some pre-generated output files can be downloaded [here](https://www.dropbox.com/s/fw7ec89hy97pluq/myoutput.tar.gz?dl=0):
```
tar -xvf myoutput.tar.gz
rm myoutput.tar.gz
```


## Preparation

Assuming you are inside the `EvoGen_course/Work` folder previously created, you can create a folder where we will put all material (and results) for the first 2 days.
```
mkdir bioinfo
cd bioinfo
mkdir output
mv ../input .
```

----------

To install the latest version of ANGSD and ngsTools the best way is to first install git, for instance as described [here](http://git-scm.com/downloads).

```
git clone --recursive https://github.com/mfumagalli/ngsTools.git
cd ngsTools
make clean
make
```
You may see a lot of warnings, but as long as they are not errors (and compilation fails), please ignore them.

To view a list of possible options, and check that everything worked for instance for ANGSD, simply type:
```
ngsTools/angsd/angsd
```

Then go back to your previous directory by typing `cd ..`.

------------

SAMtools and BCFtools can be downloaded [here](http://www.htslib.org/download/). Please use version 1.2. Then type:
```
tar xvjf samtools-1.2.tar.bz2
rm samtools-1.2.tar.bz2
cd samtools-1.2
make
cd ..
tar xvjf bcftools-1.2.tar.bz2
rm bcftools-1.2.tar.bz2
cd bcftools-1.2
make
cd ..
```


To check that everything went fine, type the following commands and check whether help messages are printed:
```
samtools-1.2/samtools
bcftools-1.2/bcftools
```

-------------

VCFtools can be obtained [here](http://sourceforge.net/projects/vcftools/files/). Version 0.1.12b is fine. Then type:
```
tar xzvf vcftools_0.1.12b.tar.gz
rm vcftools_0.1.12b.tar.gz
cd vcftools_0.1.12b
make
cd ..
./vcftools_0.1.12b/bin/vcftools
```

-------

For data filtering, we will use:

- [Scythe](https://github.com/vsbuffalo/scythe)

```
git clone https://github.com/vsbuffalo/scythe.git
cd scythe
make all
cd ..
./scythe/scythe
```

- [Sickle](https://github.com/najoshi/sickle)

```
git clone https://github.com/najoshi/sickle.git
cd sickle
make
cd ..
./sickle/sickle
```

You may want to install FastQC as well, available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc).

------

For days 3 and 4 you will need 'ms' and 'TreeMix'.
Instructions to download them cna be found [here](http://home.uchicago.edu/rhudson1/source/mksamples.html) and [here](https://code.google.com/p/treemix/).

------

Lastly, we will also assume that you have perl (for parsing) and R (with Rscript, for plotting) installed in your /usr/bin directory.  Also, python may be used.
To test if this is the case, type:
```
perl --help
R --help
Rscript --help
python --help
````
and check that this prints out help messages. Otherwise install them following instructions reported [here](http://www.perl.org/get.html) and [here](http://www.r-project.org/). You will also need some R packages to run all examples, and these are: 'ggplot2', 'VennDiagram', 'optparse'. To install a package in R, open R and type:
```
install.packages("ggplot2")
install.packages("VennDiagram")
install.packages("optparse")
```
You will be asked to select a mirror, and eventually to first install some missing dependencies. Close R by typing 'q()' and type 'n' if you are asked whether to save the current workspace.

Also [awk](http://en.wikipedia.org/wiki/AWK) will be sometimes used to parse files, in some optional examples.

------------

You will also need some scripts available [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/scripts.tar.gz). Untar these files by typing:
```
tar -xvf scripts.tar.gz
rm scripts.tar.gz
# make sure the scripts folder is inside bioinfo/ folder
```


