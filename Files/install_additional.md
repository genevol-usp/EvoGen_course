
These are various instructions to download and install some optional programs.
It is very likely you will not need to use such programs during this course.

* To run some additional optional exercises you will also need [Picard](http://picard.sourceforge.net/)
Download the latest zipped version [here](http://sourceforge.net/projects/picard/files) and unzip it.
You need [Java](http://www.java.com/en/) to run Picard tool.

* For the solution of an exercise we will use [cutadapt](https://code.google.com/p/cutadapt/).
You can use the git repository to install it from [here](https://github.com/marcelm/cutadapt).
You need Python to run it.

* [ngsDist](https://github.com/fgvieira/ngsDist) is mentioned at some point.

* In general, files will be un/compressed using gunzip, bzip2 and bunzip2, some details can be found [here](http://osxdaily.com/2012/05/29/create-extract-bz2-mac-os-x/) and [here](http://linux.about.com/library/cmd/blcmdl1_bunzip2.htm).

* For some optional exercises, a Perl script will use Statistics::Distributions package, available for download [here](http://search.cpan.org/~mikek/Statistics-Distributions-1.02/Distributions.pm).
The easiest way to install it is to run:
```
sudo cpan Statistics::Distributions
```
If it fails, download the .tar.gz file, then unzipped it and install it:
```
tar -xvzf Statistics-Distributions-1.02.tar.gz
cd Statistics-Distributions-1.02
perl Makefile.PL
make test
cd..
```
You need to make sure that the package has been correctly added your Perl directory.

Likewise you may need to install this additional package [IO-Compress](http://search.cpan.org/~pmqs/IO-Compress-2.064/lib/IO/Compress/Bzip2.pm):
```
sudo cpan IO::Compress::Bzip2
# if it does not work try: sudo cpan force install IO::Compress::Bzip2
```
or manually:
```
tar -xvzf IO-Compress-2.064.tar.gz
cd IO-Compress-2.064
perl Makefile.PL
cd ..
```
and Getopt which should be already installed, otherwise see [here](http://search.cpan.org/~jhi/perl-5.8.1/lib/Getopt/Std.pm).
Before being worried that something failed, read below for a quick test to check that this perl script indeed works.

To test that you have all required packages in Perl type:
```
perl scripts/SNPcleaner.pl --help
```
and see if it prints out something or errors.

* Optionally, you can download [GATK](http://www.broadinstitute.org/gatk/) and [FreeBayes](https://github.com/ekg/freebayes), as they will be briefly mentioned and discussed as supplementary information.


