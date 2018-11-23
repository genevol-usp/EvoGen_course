
# plink files
fin="../Data/hapmap.frq.strat"
fbim="../Data/hapmap.bim"

source("../Scripts/functions.R")

# read allele frequencies and SNPs coordinates
res=read.table(fin, stringsAsFactors=F, header=T)
bim=read.table(fbim, stringsAsFactors=F, sep="\t")

snps=unique(res$SNP)

# record FST YRI-CEU YRI-CHB CEU-CHB
heads=c("chrom", "pos", "rs", "cpos", "MAF.YRI", "MAF.CEU", "MAF.CHB", "FST.YRI-CEU", "FST.YRI-CHB", "FST.CEU-CHB")
cat(heads, sep="\t", file="Results/hapmap.fst")
cat("\n", file="Results/hapmap.fst", append=T)

offset=0
current_chrom=1
last_cpos=0

for (s in 1:length(snps))
{

	if ((s %% 500)==0) cat(s,"/",length(snps),"\n")
	
        val=unlist(c(bim[which(bim[,2]==snps[s]),])[c(1,4,2)])

	if (!is.na(as.numeric(val[1]))) {

		if (as.numeric(val[1])>current_chrom) {
			current_chrom=as.numeric(val[1]);
			offset=last_cpos;
		}

		val=c(val, as.numeric(val[2])+offset)

	        ind=which(res$SNP==snps[s]);

		val=c(val, res[ind[which(res[ind,]$CLST=="YRI")],]$MAF, res[ind[which(res[ind,]$CLST=="CEU")],]$MAF, res[ind[which(res[ind,]$CLST=="CHB")],]$MAF)

		ii=which(res[ind,]$CLST=="YRI");  jj=which(res[ind,]$CLST=="CEU")
	        tmp=reynolds( res[ind[ii],]$MAF, res[ind[jj],]$MAF, res[ind[ii],]$NCHROBS, res[ind[jj],]$NCHROBS);
       		val=c(val, tmp[3]) # FST     

        	jj=which(res[ind,]$CLST=="CHB")
        	tmp=reynolds( res[ind[ii],]$MAF, res[ind[jj],]$MAF, res[ind[ii],]$NCHROBS, res[ind[jj],]$NCHROBS);
        	val=c(val, tmp[3]) # FST

        	ii=which(res[ind,]$CLST=="CEU");  jj=which(res[ind,]$CLST=="CHB")
        	tmp=reynolds( res[ind[ii],]$MAF, res[ind[jj],]$MAF, res[ind[ii],]$NCHROBS, res[ind[jj],]$NCHROBS);
        	val=c(val, tmp[3]) # FST

		cat(val, sep="\t", file="Results/hapmap.fst", append=T)
        	cat("\n", file="Results/hapmap.fst", append=T)

		last_cpos=as.numeric(val[4])
        	rm(val)

	}

}

