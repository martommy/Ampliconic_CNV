
require(reshape2)
require(ggplot2)

######haplotype networkanalysis####

#load cnv data
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#round mean copy number to nearest integer

copy.integer<-dat4
copy.integer[,c(2:10)]<-apply(copy.integer[,c(2:10)],2,round)

##testing effect of rounding mean copy number to integer on number of haplotypes
mcopy.number<-melt(dat4[,c(1:10,13)],id.vars = c("IID","major_haplo"))
colnames(mcopy.number)[c(3,4)]<-c("Gene","copy_number")

#create limits of window within which values will be randomly rounded up or down
#lower limit
mcopy.number$floor<-floor(mcopy.number$copy_number)+0.25
#upper limit
mcopy.number$ceiling<-ceiling(mcopy.number$copy_number)-0.25
#new column indicating which observations fall within this window
mcopy.number$hairy<-NA
mcopy.number$hairy[which(mcopy.number$copy_number<mcopy.number$ceiling & mcopy.number$copy_number>mcopy.number$floor)]<-1
mcopy.number$hairy[-which(mcopy.number$copy_number<mcopy.number$ceiling & mcopy.number$copy_number>mcopy.number$floor)]<-0

#function to randomly round up or down an observation x
rand.round<-function(x){
  flip<-rbinom(1,1,0.5)
  if(flip==1){y=ceiling(x)}
  if(flip==0){y=floor(x)}
  return(y)
}

#function to randomly round up or down an observation x
rand.round<-function(x){
  flip<-rbinom(1,1,0.5)
  if(flip==1){y=ceiling(x)}
  if(flip==0){y=floor(x)}
  return(y)
}

#apply rand.round function to data 100 times
rand.mat<-matrix(NA,900,100)
for(i in 1:100){
  rand.mat[which(mcopy.number$hairy==1),i]<-sapply(mcopy.number[which(mcopy.number$hairy==1),'copy_number'],rand.round)
  rand.mat[which(mcopy.number$hairy==0),i]<-round(mcopy.number[which(mcopy.number$hairy==0),'copy_number'])
}

#cleanup
rand.mat<-as.data.frame(rand.mat)
rand.mat$IID<-mcopy.number$IID
rand.mat$Gene<-mcopy.number$Gene
mrand.mat<-melt(rand.mat,id.vars=c("Gene","IID"))
dmrand.mat<-dcast(mrand.mat,IID+variable~Gene,value.var="value")
dmrand.mat<-dmrand.mat[order(dmrand.mat$variable),]

#write randomly rounded haplotypes to file
write.table(dmrand.mat,'random_round_haplotypes.txt',sep="\t",col.names=T,row.names=F,quote=F)

######calculating number of differences between pairs of haplotypes within and between haplogroups######

#create vector of haplogroups
haplos<-levels(dat4$major_haplo)

#randomly select from within haplogroups
wn.diffs<-matrix(NA,1000,9)
wn.haplo<-matrix(NA,1000,1)
for(i in 1:1000){
  test.haplo<-sample(haplos,1)
  wn.haplo[i,1]<-test.haplo
  test.dat<-copy.integer[which(copy.integer$major_haplo==test.haplo),]
  test.dat<-test.dat[sample(nrow(test.dat),2),]
  wn.diffs[i,]<-as.matrix(abs(test.dat[1,c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")]-test.dat[2,c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")]))
}

#randomly select from between haplogroups
bw.diffs<-matrix(NA,1000,9)
bw.haplo<-matrix(NA,1000,2)
for(i in 1:1000){
  test.haplos<-sample(haplos,2)
  bw.haplo[i,]<-test.haplos
  test.dat1<-copy.integer[which(copy.integer$major_haplo==test.haplos[1]),]
  test.dat2<-copy.integer[which(copy.integer$major_haplo==test.haplos[2]),]
  test.dat<-rbind(test.dat1[sample(nrow(test.dat1),1),],test.dat2[sample(nrow(test.dat2),1),])
  bw.diffs[i,]<-as.matrix(abs(test.dat[1,c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")]-test.dat[2,c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")]))
}

colnames(wn.diffs)<-colnames(bw.diffs)<-c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")
bw.diffs<-as.data.frame(bw.diffs)
wn.diffs<-as.data.frame(wn.diffs)

bw.diffs$comparison<-"Between Haplogroups"
wn.diffs$comparison<-"Within Haplogroups"

comb.diffs<-rbind(bw.diffs,wn.diffs)
mcomb.diffs<-melt(comb.diffs,id.vars=c("comparison"))
colnames(mcomb.diffs)<-c("comparison","Gene","Difference")
fig_8<-ggplot()+geom_boxplot(data=mcomb.diffs,aes(Gene,Difference,fill=comparison))+theme_bw()+labs(y="Copy number difference b/w two randomly picked haplotypes")

ggsave("../Figures/Fig_8.pdf",fig_8,height=5,width=7)
