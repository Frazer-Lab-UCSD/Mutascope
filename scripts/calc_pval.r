options(warn=-1)

in_file = commandArgs()[4]
out_file = commandArgs()[5]


in.pval <- read.table(file=in_file)

pvals = -1*pbinom(in.pval[,3], in.pval[,4], in.pval[,5], lower.tail=FALSE,log.p=TRUE)

write.table(cbind(in.pval, pvals), file=out_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



