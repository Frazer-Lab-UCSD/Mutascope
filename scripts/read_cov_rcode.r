
infile = commandArgs()[4]
ms_file = commandArgs()[5]
unif_file = commandArgs()[6]
sens_file = commandArgs()[7]
pdf_file = commandArgs()[8]

x <- read.table(file=infile)

ms = c("Mean", "SD")
ms = rbind(ms, c(mean(x[,5]), sd(x[,5])))

min = mean(x[,5])/2
max = mean(x[,5])*2

unif = c("Cov2X")
unif = rbind(unif, length(which(x[,5] >= min & x[,5] <= max))/dim(x)[1])


sens = rbind(c("< 50X", "< 100X", "< 500X", "< 1000X", ">= 1000X"), c(length(which(x[,5]<50)), length(which(x[,5]<100)), length(which(x[,5]<500)), length(which(x[,5]<1000)), length(which(x[,5]>=1000))), c(length(which(x[,5]<50))/dim(x)[1], length(which(x[,5]<100))/dim(x)[1], length(which(x[,5]<500))/dim(x)[1], length(which(x[,5]<1000))/dim(x)[1], length(which(x[,5]>=1000))/dim(x)[1]))

write.table(ms, file=ms_file, row.names=F, col.names=F, quote=F, sep="\t")
write.table(unif, file=unif_file, row.names=F, col.names=F, quote=F, sep="\t")
write.table(sens, file=sens_file, row.names=F, col.names=F, quote=F, sep="\t")

pdf(file=pdf_file)
plot(ecdf(x[,5]/mean(x[,5])), xlab="Normalized Read Counts", ylab="% of bases with <= X Normalized Read Count", main="", xlim=c(0,5))
dev.off()

