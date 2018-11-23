
args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

res=read.table(fin, sep="\t", head=T)

# odd-even colors for chroms
cols=rep("grey", nrow(res))
cols[which( (res$chrom %% 2) == 1)]="lightgrey"

jpeg(fout)
par(mfrow=c(3,1))

for (i in 8:10) {

	# cutoff outliers
	#ths=c()
	#ths=c(ths, as.numeric(quantile(res[,i], seq(0.99,0.99,0.01), na.rm=T)))
	#ths=c(ths, as.numeric(quantile(res[,i], seq(0.999,0.999,0.001), na.rm=T)))

	plot(x=res$cpos, y=res[,i], col=cols, frame=F, xlab="", xaxt="n", ylab="FST", main=colnames(res)[i], ylim=c(0,1), pch=16)
	#abline(h=ths, lty=2)

}

dev.off()



