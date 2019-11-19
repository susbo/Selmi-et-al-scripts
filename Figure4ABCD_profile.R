args = commandArgs(TRUE)

setwd("/servers/frye-bioinf/smb208/Nsun6/Riboseq/profiles")

my_mean = function(list) {
	mean(list,trim=0.025)
}

file = args[1]
from = args[2]

filename = sub(".replicate.gz","",file)
replicates = 2
colors = c("#d7191cFF","#1d587c")
for (i in seq(1,8,1)) colors[i] = gsub("FF$","CC",colors[i])
ltys = c(1,1)

dat.plus.sense = read.table(gzfile(paste("matrix/",from,"/",file,".plus.sense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)
dat.minus.sense = read.table(gzfile(paste("matrix/",from,"/",file,".minus.sense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)
dat.plus.anti = read.table(gzfile(paste("matrix/",from,"/",file,".plus.antisense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)
dat.minus.anti = read.table(gzfile(paste("matrix/",from,"/",file,".minus.antisense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)

max_coord = (dim(dat.plus.sense)[2]-6)/replicates
for (repl in 1:replicates) {
	region = (max_coord*(repl-1)+1):(max_coord+max_coord*(repl-1))
	dat.minus.sense[,6+region] = rev(dat.minus.sense[,6+region])
	dat.minus.anti[,6+region] = rev(dat.minus.anti[,6+region])
}

dat = rbind(dat.plus.sense,dat.minus.sense)
dat.anti = rbind(dat.plus.anti,dat.minus.anti)

rows=dim(dat)[1]
outname = paste(filename," [n=",rows,"]",sep="")

dat[dat=="NaN"]=0

pdf(paste("R/",from,"2/Profile.",filename,".pdf",sep=""),width=4,height=3)
par( mfrow=c(1,1) )
par(oma = c(0,0,0,0) + 0.5)
par(mar = c(2.5,2.5,3,0), mgp=c(1.5,0.6,0))

y.all = list()
max = 0; min=100000;
for (repl in 1:1) {
	region = (max_coord*(repl-1)+1):(max_coord+max_coord*(repl-1))
	data_matrix = dat[,6+region]
	y = as.numeric(sapply(data_matrix,my_mean))
	x=1:length(y)-length(y)/2
	x=x*5

	y.all[[repl]] = y
	max = max(max,y)
	min = min(min,y)
}

ylim=c(min,max)

for (repl in 1:replicates) {
	lty = ltys[(repl-1)%%2+1]

	if (repl == 1) {
		plot(x,y.all[[repl]],type="l",bty="n",xlab="miCLIP site",ylab="Fraction with feature",lwd=1,col=colors[repl],main=outname,lty=lty,ylim=ylim,xpd=NA)
	} else {
		# Show only CDS, not exon annotations
#		lines(x,y.all[[repl]],type="l",lwd=1,col=colors[repl],xpd=NA,lty=lty,xpd=NA)
	}
}
segments(0,0,0,100000,lty=2)
legend("topright",col=colors,c("CDSs"),lwd=2,bty="n",cex=0.8)

dev.off()
