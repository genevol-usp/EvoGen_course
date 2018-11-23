
args=commandArgs(T)
tmix=args[1]
fin=args[2]
rm(args)

source(paste(tmix,"/src/plotting_funcs.R", sep="",collapse=""))
pdf(paste(fin,".pdf",sep="",collapse=""))
par(mfrow=c(1,2))
plot_tree(fin)
plot_resid(fin, "../Data/poporder.txt")
dev.off()




