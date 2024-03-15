options(warn=-1)


in_file = commandArgs()[4]
out_file = commandArgs()[5]

in.rgb <- read.table(file=in_file)

options(warn=-1)

if(dim(matrix(as.double(strsplit(toString(in.rgb[1,3]),"_")[[1]]), nrow=2))[2]>=2){
		x = chisq.test(matrix(as.double(strsplit(toString(in.rgb[1,3]),"_")[[1]]), nrow=2))$p.value
		if(is.na(x)){
			rgb.pvals = 1
		}else{
			rgb.pvals = x			
		}
}else{
	rgb.pvals = c(1)
}

if(dim(in.rgb[1] >= 2)){
	for(i in 2:dim(in.rgb)[1]){
		if(dim(matrix(as.double(strsplit(toString(in.rgb[i,3]),"_")[[1]]), nrow=2))[2] == 1){
			rgb.pvals = rbind(rgb.pvals, 1)
		}else{
			x = chisq.test(matrix(as.double(strsplit(toString(in.rgb[i,3]),"_")[[1]]), nrow=2))$p.value
			if(is.na(x)){
				rgb.pvals = rbind(rgb.pvals, 1)
			}else{
				rgb.pvals = rbind(rgb.pvals, x)			
			}
		}
	}
}

write.table(cbind(in.rgb[,1:2], as.double(rgb.pvals)), file=out_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


