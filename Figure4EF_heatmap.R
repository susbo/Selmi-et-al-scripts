setwd("/servers/frye-bioinf/smb208/Nsun6/Riboseq/profiles")

set = "all"
file = paste("hg38.",set,sep="")
from = "all"

replicates = 7
colors = c("#ea484bFF","#d7191cFF","#a01315FF","#5b0b0cFF","#5aa8d8","#2b83baFF","#1d587c")
for (i in seq(1,8,1)) colors[i] = gsub("FF$","CC",colors[i])
ltys = c(1,1)

dat.plus.sense = read.table(gzfile(paste("matrix/",from,"/",file,".plus.sense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)
dat.minus.sense = read.table(gzfile(paste("matrix/",from,"/",file,".minus.sense.gz",sep="")),sep="\t",skip=1,stringsAsFactors = FALSE)

max_coord = (dim(dat.plus.sense)[2]-6)/replicates
for (repl in 1:replicates) {
	region = (max_coord*(repl-1)+1):(max_coord+max_coord*(repl-1))
	dat.minus.sense[,6+region] = rev(dat.minus.sense[,6+region])
}

dat = rbind(dat.plus.sense,dat.minus.sense)

mat = as.matrix(dat[,7:dim(dat)[2]])
mat = log(mat+1)
mat[mat>0] = 1
mat = t(mat)

image(mat,col  = gray(c(1,0)))

left = c(); right = c();
for (i in 1:replicates) { left = c(left,seq(max_coord*(i-1)+1,max_coord*(i-0.5))); right = c(right,seq(max_coord*(i-0.5)+1,max_coord*(i-0))); }
order = order(colSums(mat[left,])-colSums(mat[right,]))

threshold = sd(colSums(mat[left,])-colSums(mat[right,]))/sqrt(dim(mat)[2])*1.96
end = (colSums(mat[left,])-colSums(mat[right,])) > threshold
start = (colSums(mat[left,])-colSums(mat[right,])) < -threshold
mid = ((colSums(mat[left,])-colSums(mat[right,])) > -threshold & (colSums(mat[left,])-colSums(mat[right,])) < threshold)

write.table(cbind(dat[end,1:3],rep("1",sum(end)),c(rep("+",dim(dat.plus.sense)[1]),rep("-",dim(dat.minus.sense)[1]))[end]),"regions/hg38.end.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(cbind(dat[start,1:3],rep("1",sum(start)),c(rep("+",dim(dat.plus.sense)[1]),rep("-",dim(dat.minus.sense)[1]))[start]),"regions/hg38.start.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(cbind(dat[mid,1:3],rep("1",sum(mid)),c(rep("+",dim(dat.plus.sense)[1]),rep("-",dim(dat.minus.sense)[1]))[mid]),"regions/hg38.mid.bed",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

pdf(paste("Heatmap.",set,".direction.pdf",sep=""),width=4,height=2)
par( mfrow=c(1,1) )
par(oma = c(0,0,0,0) + 0.5)
par(mar = c(1,1,0,0)+0.5, mgp=c(0.5,0.6,0))
image(mat[,order],col = gray(c(1,0)),xaxt="n",yaxt="n",ylab=paste(set," miCLIP sites",sep=""))
segments(1:7/7,0,1:7/7,1)
text(1:7/7-1/14,-0.05,c("KO.D8.2","KO.H95.2","KO.H96.1","KO.H97.2","WT.C8.2","WT.H9.1","WT.H9.2"),xpd=NA,cex=0.6)
dev.off()

png(paste("Heatmap.",set,".direction.png",sep=""),width = 4, height = 2, units = "in", pointsize = 10, bg = "white",  res = 300)
par( mfrow=c(1,1) )
par(oma = c(0,0,0,0) + 0.5)
par(mar = c(1,1,0,0)+0.5, mgp=c(0.5,0.6,0))
image(mat[,order],col = gray(c(1,0)),xaxt="n",yaxt="n",ylab=paste(set," miCLIP sites",sep=""))
segments(1:7/7,0,1:7/7,1)
text(1:7/7-1/14,-0.05,c("KO.D8.2","KO.H95.2","KO.H96.1","KO.H97.2","WT.C8.2","WT.H9.1","WT.H9.2"),xpd=NA,cex=0.6)
dev.off()
