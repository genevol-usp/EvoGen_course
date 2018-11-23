
# Rscript plotTreemix_multi.R Treemix_dir results.likes

args=commandArgs(T)

args=commandArgs(T)
tmix=args[1]
flikes=args[2]
rm(args)

source(paste(tmix,"/src/plotting_funcs.R", sep="",collapse=""))

likes=readLines(flikes)

tmp=as.character(unlist(strsplit(unlist(strsplit(likes[seq(1,length(likes),3)],paste("==> Results/hapmap.frq.strat.tree.",sep="",collapse=""))),split=".llik <==")))

inds=matrix(as.numeric(unlist(strsplit(tmp,split="\\."))),ncol=2,byrow=T)
likes=as.numeric(likes[seq(2,length(likes),3)])
rm(tmp)

ms=sort(unique(inds[,1]))

for (m in ms) {

	ii=which(inds[,1]==m)
	so=sort(likes[ii],dec=T)[1:10]
	im=which.max(likes[ii])
	irep=inds[ii[im],2]
	 uu=sort(unique(likes[ii]),dec=T)

	cat("\n",m,im,irep)

	pdf(file=paste("Results/treemix.m",m,".FULL.pdf", sep="", collapse=""))

	hist(so, xlab="", main="Top 10 likelihoods", sub=paste("min of 100 reps is",min(likes[ii]),";min on top10 is",min(so),"\nmax(this) is",max(likes[ii]),";2nd max unique (last plot) is",uu[2]))
	abline(v=so[1],lty=2)

	plot_tree(paste("Results/hapmap.frq.strat.tree.",m,".",irep, sep="", collapse=""))

	plot_resid(paste("Results/hapmap.frq.strat.tree.",m,".",irep, sep="", collapse=""), "Data/poporder.txt")	

	if (length(uu)>1) {
		im=which(likes[ii]==uu[2])[1]
		irep=inds[ii[im],2]
		plot_tree(paste("Results/hapmap.frq.strat.tree.",m,".",irep, sep="", collapse=""))
	}
	
	dev.off()

}

for (m in ms) {

        ii=which(inds[,1]==m)
        so=sort(likes[ii],dec=T)[1:10]
        im=which.max(likes[ii])
        irep=inds[ii[im],2]
         uu=sort(unique(likes[ii]),dec=T)

        cat("\n",m,im,irep)

        pdf(file=paste("Results/treemix.m",m,".pdf", sep="", collapse=""))

        plot_tree(paste("Results/hapmap.frq.strat.tree.",m,".",irep, sep="", collapse=""))

        dev.off()

}






