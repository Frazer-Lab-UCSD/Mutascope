options(warn=-1)

in_file = commandArgs()[4]
out_file = commandArgs()[5]

x = read.table(file=in_file)

x[,3] = x[,3]/100
x[,4] = x[,4]/100

ng = dim(x[x[,5] != 2,])[1]
ns = dim(x[x[,5] == 2,])[1]

pdf(file=out_file)

par(mar=c(5,5,4,1))

plot(x[(x[,6]=="SNP" & x[,5] != 2),3],x[(x[,6]=="SNP" & x[,5] != 2),4], xlim=c(0,1), ylim=c(0,1), col="blue", ylab="Tumor Allelic Fraction", xlab="Normal Allelic Fraction", main=c(paste(ng, "Germline Variants", sep=" "), paste(ns, "Somatic Variants", sep=" ")), pch=1, cex=2, cex.axis=1.5, cex.lab=1.5)

points(x[(x[,6]=="SNP" & x[,5] == 2),3],x[(x[,6]=="SNP" & x[,5] == 2),4], col="red", pch=1, cex=2)

points(x[(x[,6]!="SNP" & x[,5] != 2),3],x[(x[,6]!="SNP" & x[,5] != 2),4], pch=0, col="blue", cex=2)

points(x[(x[,6]!="SNP" & x[,5] == 2),3],x[(x[,6]!="SNP" & x[,5] == 2),4], pch=0, col="red", cex=2)

legend("bottomright", legend=c("SNP", "Indel", "Germline", "Somatic"), col=c("black", "black", "blue", "red"), pch=c(1, 0, 15, 15), pt.cex=c(1.8,1.8,1.8,1.8), cex=1.3)

dev.off()

