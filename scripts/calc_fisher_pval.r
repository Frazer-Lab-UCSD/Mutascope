options(warn=-1)

in_file = commandArgs()[4]
out_file = commandArgs()[5]

in.fish <- read.table(file=in_file)

pvals.less =  fisher.test(rbind(c(in.fish[1,3], in.fish[1,4]),c(in.fish[1,5], in.fish[1,6])), alternative = "less")$p.val

if(dim(in.fish)[1] > 1){
	for(i in 2:dim(in.fish)[1]){
		
		pvals.less = c(pvals.less, fisher.test(rbind(c(in.fish[i,3], in.fish[i,4]),c(in.fish[i,5], in.fish[i,6])), alternative = "less")$p.val)

	}
}

pvals.greater =  fisher.test(rbind(c(in.fish[1,3], in.fish[1,4]),c(in.fish[1,5], in.fish[1,6])), alternative = "greater")$p.val

if(dim(in.fish)[1] > 1){
	for(i in 2:dim(in.fish)[1]){
		
		pvals.greater = c(pvals.greater, fisher.test(rbind(c(in.fish[i,3], in.fish[i,4]),c(in.fish[i,5], in.fish[i,6])), alternative = "greater")$p.val)

	}
}

write.table(cbind(in.fish, pvals.less, pvals.greater), file=out_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


