options(warn=-1)

cmd_args = commandArgs()[4]
out = commandArgs()[5]

is_installed = "mclust" %in% installed.packages()[,"Package"]

if(is_installed){

	library(mclust)

	tmp.abq = read.table(file=cmd_args)

	stmp.abq = sort(tmp.abq[,1])
	stmp.abq.mclust <- Mclust(stmp.abq, G=2)
	stmp.abq.mclust.inSecond = stmp.abq.mclust$z[,2] >= .99


	write(min(stmp.abq[stmp.abq.mclust.inSecond]), file=out)

}else{
	write(25, file=out)
}
