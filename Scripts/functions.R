
library("grid")
library("maps")
library("spam")
library("fields")


plot2DSFS<-function(sfs, xlab="", ylab="", main="") {

	# this function plots a 2D-SFS in log10 scale; it requires "fields" package
	nchroms_pop1<-nrow(sfs)-1
        nchroms_pop2<-ncol(sfs)-1

	sfs=log10(sfs)
	sfs[which(is.na(sfs))]=0
	sfs[which(sfs==(-Inf))]=0	

	brk=seq(min(sfs),max(sfs)+0.3,0.3)
	lab.brk=round(10^(seq(min(sfs),max(sfs)+0.3,0.3)))
	lab.brk[1]=0
	cols=rainbow(length(brk)-1)
	cols[1]="white"

	image.plot(x=seq(0,nchroms_pop1), y=seq(0,nchroms_pop2), z=sfs, breaks=brk, col=cols, lab.breaks=lab.brk, xlab=xlab, ylab=ylab, main=main)




}


fold2DSFS<-function(unfolded) {

	# this function takes as input an unfolded spcetrum and returns a folded one

        folded<-unfolded

        sample_size_pop1<-(nrow(unfolded)-1)/2
        sample_size_pop2<-(ncol(unfolded)-1)/2

        freq_at50<-sample_size_pop1+sample_size_pop2
 
        for (i in 1:nrow(unfolded)) {
                for (j in 1:ncol(unfolded)) {
                        daf1=(i-1) ##take the values chosen from the rows and subtract 1 for the true daf
                        daf2=(j-1) ##take the values chosen from the columns and subtract 1 for the true daf
                        if ((daf1+daf2)>freq_at50) { ## when value is major
                                minor1=(2*sample_size_pop1)-daf1; ## converts to minor
                                minor2=(2*sample_size_pop2)-daf2; ## converts to minor
                                folded[(minor1+1),(minor2+1)]=folded[(minor1+1), (minor2+1)]+unfolded[(daf1+1), (daf2+1)]; ##add 1 back to the values to get to the true spot in the matrix and replace it by adding the values that were in the majors to the values in the minors
                                folded[(daf1+1),(daf2+1)]=NA ##replace all values where majors originally were with NA
                        }
                }
        }

        folded

}


doFST<-function(sfs) {

	# this function compute the FST from a 2D-SFS, using the Reynold's et al estimator

        sfs=sfs/sum(sfs, na.rm=T)
        nind1=(nrow(sfs)-1)/2
        nind2=(ncol(sfs)-1)/2

        nums=denoms=fsts=matrix(NA, ncol=ncol(sfs), nrow=nrow(sfs))
        for (i in 1:nrow(sfs)) {
                for (j in 1:ncol(sfs)) {
                        f1=(i-1)/(nind1*2); f2=(j-1)/(nind2*2)
                        tmp=reynolds(f1,f2,nind1*2,nind2*2)
                        nums[i,j]=tmp[1]
                        denoms[i,j]=tmp[2]
                        fsts[i,j]=tmp[3]
                }
        }

        fst=sum(sfs*nums,na.rm=T)/sum(sfs*denoms,na.rm=T)
        fst
}

reynolds<-function(pl1,pl2,n1,n2) {

	# this functions compute the FST estimator from Reynolds et al. from sample allele frequencies
	somma=sommaden=0;

	alfa1=1-((pl1^2)+((1-pl1)^2))
	alfa2=1-((pl2^2)+((1-pl2)^2))
	Al = (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) - (((n1+n2)*(n1*alfa1+n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))
	AlBl= (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) + (((4*n1*n2 - n1 - n2)*(n1*alfa1 + n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

	if (!is.na(Al) & !is.na(AlBl)) {
  		somma=somma+Al
  		sommaden=sommaden+AlBl
 	}

	if (somma==0 & sommaden==0) {
      		reyn=NA
 	} else {
        	reyn=somma/sommaden
        	if(reyn<0) reyn=0
 	}

 	c(somma, sommaden, reyn)

}


fromMStoSFSwith1site<-function(msfile, nr_repetitions, nr_chroms_pop1, nr_chroms_pop2 ) {

	# this functions read ms output file (if only 1 site is simulated) and returns a matrix containing the 2D-SFS
	# originally written by S.D. Gopal

        nr_chromosomes<-nr_chroms_pop1+nr_chroms_pop2

	# read ms file
        sequencedata<-readLines(msfile)##read the ms command output into R
        sequencedata<-suppressWarnings(as.numeric(sequencedata)) ##make output numeric
        sequencedata<-sequencedata[!is.na(sequencedata)] ##remove NAs which are words in this case

        sequencedata_matrix<-matrix(sequencedata, nrow= (nr_chromosomes), ncol= (nr_repetitions))

        chroms_pop1<-sequencedata_matrix[1:nr_chroms_pop1, (1:nr_repetitions)] ### subset pop1 data
        chroms_pop2<-sequencedata_matrix[(nr_chroms_pop1+1):(nr_chroms_pop1+nr_chroms_pop2), (1:nr_repetitions)] ### subset pop2 data

        sfs<-matrix(0, nrow= (nr_chroms_pop1+1), ncol= (nr_chroms_pop2+1)) ### create a matrix filled with zeroes

        for (j in 1:nr_repetitions) {

                ### sum the value across each site in pop1 and pop2 respectively to calculate the derived allele frequencies 

                daf_pop1<-sum(chroms_pop1[,j]) + 1 # +1 to convert into coordinates
                daf_pop2<-sum(chroms_pop2[,j]) + 1

                sfs[daf_pop1, daf_pop2]<-sfs[daf_pop1, daf_pop2]+1 ## place summations in the matrix 

        }

        sfs


}





dopbs<-function(Fgc,Fge,Fce) {
        Tgc= -log(1-Fgc)
        Tge= -log(1-Fge)
        Tce= -log(1-Fce)
        pbs= (Tgc + Tge - Tce)/2
        pbs
}

readMs<-function(filein, nsam, len=c()) {

 lfile<-readLines(con=filein)
 inds<-which(lfile=="//")+3
 inde<-which(lfile=="//")+3+nsam-1
 # cat("Sample:",length(inds),"\n")
 indpos<-inds-1

 res<-pos<-list()
 for (i in 1:length(inds)) {
  res[[i]]<-trim(lfile[inds[i]:inde[i]])
  ss<-strsplit(lfile[indpos[i]],split=" ")[[1]]
  pos[[i]]<-as.numeric(ss[2:length(ss)])
  if (length(len)>0) pos[[i]]<-as.integer(pos[[i]]*len)
 }

 readMs<-list(hap=res, pos=pos)

}

trim<-function (sequence, inner = FALSE) {
## Author: Giorgia Menozzi
    sequence <- sub("^ +", "", sequence)
    sequence <- sub(" +$", "", sequence)
    if (inner) {
        sequence <- gsub(" +", "", sequence)
    }
    sequence
}


chroms2fst<-function(haplos, verbose=FALSE) {

 # calcola numerosità campione
 npop=length(haplos)
 new_nsam=c(); for (p in 1:npop) new_nsam[p]=length(haplos[[p]])

 # output
 reyn=c()

 # calcola frequenze alleliche
 sfs=countFreq(haplos, haplos[[1]][1], fixed.na=FALSE, plot=FALSE)

 name=c()

 t=0
 for (p1 in 1:(npop-1)) {
  for (p2 in (p1+1):npop) {

   name=c(name, paste(p1,p2,sep="-",collaps=""))

   t=t+1
   somma=0; sommaden=0
   n1=new_nsam[p1]; n2=new_nsam[p2]

   for (i in 1:length(sfs[[1]])) {

    pl1=sfs[[p1]][i]
    pl2=sfs[[p2]][i]
    if (verbose) cat("\n\npl",i,pl1,pl2)

    alfa1=1-((pl1^2)+((1-pl1)^2))
    alfa2=1-((pl2^2)+((1-pl2)^2))

    if (verbose) cat("\nalfa", alfa1, alfa2)

    Al = (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) - (((n1+n2)*(n1*alfa1+n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

    AlBl= (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) + (((4*n1*n2 - n1 - n2)*(n1*alfa1 + n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

if (verbose) cat("\nA\tB\tA+B", Al, AlBl-Al, AlBl)

    if (!is.na(Al) & !is.na(AlBl)) {
     somma=somma+Al
     sommaden=sommaden+AlBl
    }

   }

# se 1 pop e' singleton e altra fissa e' 0
   if (somma==0 & sommaden==0) reyn[t]=NA else reyn[t]=somma/sommaden

  }
 }
#
#  cat(" NA:",length(which(is.na(reyn[g,]))))
#  cat(" ZERO:",length(which(reyn[g,]==0)))
 names(reyn)=name
 reyn

} # fine function




countFreq<-function(haplos, outgroup=c(), plot=F, fixed.na=FALSE) {

 if (!is.list(haplos)) haplos<-list(haplos)
 npop<-length(haplos)

 maf<-list()

 if (plot) { x11(); par(mfrow=c(npop,1)) }
 for (p in 1:npop) {

  nsam<-length(haplos[[p]])
  len<-nchar(haplos[[p]][1])
  maf[[p]]<-rep(NA, len)
  for (l in 1:len) {
   hs<-substring(haplos[[p]],l,l)
   hs<-hs[which(hs!="N")]
   if (length(outgroup)>0) {
    maf[[p]][l]<-length(which(hs!=substring(outgroup,l,l)))/length(which(hs!="N"))
   } else {
    maf[[p]][l]<-min(table(hs))/length(which(hs!="N"))
   }
  }
  if (fixed.na==TRUE) maf[[p]][which(maf[[p]]==1)]<-NA
  if (plot) if (length(outgroup)==0) hist(maf[[p]],breaks=10,xlim=c(0,0.5),sub=paste("nsam/2 is",nsam/2),main="MAF histogram") else hist(maf[[p]],breaks=10,xlim=c(0,1),sub=paste("nsam/2 is",nsam/2),main="DAF histogram")
 }

 countFreq<-maf

}



