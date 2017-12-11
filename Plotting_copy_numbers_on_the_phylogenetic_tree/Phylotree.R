
require(ggplot2)
require(ggtree)
require(plyr)

#read ampliconic gene numbers - cleaned
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#read y haplogrup tree
ytree<-read.tree('../Data_files/ytree.nwk')

#convert tree to dataframe for easier plotting and labeling
ydat<-fortify(ytree)

#split label into individual id and haplogroup info
ydat$IID<-sapply(ydat$label,function(x){unlist(strsplit(x,split="_"))[1]})
ydat$haplogroup<-sapply(ydat$label,function(x){unlist(strsplit(x,split="_"))[2]})
ydat$major_haplo<-sapply(ydat$haplogroup,function(x){unlist(strsplit(x,split=""))[1]})

#add copy number data to phylogeny table and round off to 2 decimal places
ydat<-join(ydat,dat4,by="IID")

#plot y phylogeny with haplogroups colored
groupInfo<-split(ydat$label,ydat$major_haplo)
ytree<-groupOTU(ytree,groupInfo)

p<-ggtree(ytree,aes(color=group))+geom_segment(data=ydat[which(is.na(ydat$haplogroup)=="FALSE"),],aes(x=x,xend=0.09,y=y,yend=y,color=major_haplo),linetype="dotted")

#add scale
p<-p+geom_treescale(y=-2,offset=-2)

ydat[,c('BPY','CDY','DAZ','HSFY','PRY','RBMY','TSPY','VCY','XKRY')]<-apply(ydat[,c('BPY','CDY','DAZ','HSFY','PRY','RBMY','TSPY','VCY','XKRY')],2,function(x){format(round(x,2),nsmall=2)})
ydat<-ydat[which(is.na(ydat$major_haplo)=="FALSE"),]

#add copy number information - Figure 2
fig_2<-p+geom_text(data=ydat,aes(x=0.096,label=CDY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.100,label=BPY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.104,label=DAZ,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.108,label=PRY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.112,label=RBMY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.116,label=HSFY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.120,label=XKRY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.124,label=VCY,y=y,color=major_haplo))+geom_text(data=ydat,aes(x=0.128,label=TSPY,y=y,color=major_haplo))+geom_text(aes(x=0.096,y=106,label="CDY"),color="black",angle=90,size=10)+geom_text(aes(x=0.100,y=106,label="BPY"),color="black",angle=90,size=10)+geom_text(aes(x=0.104,y=106,label="DAZ"),color="black",angle=90,size=10)+geom_text(aes(x=0.108,y=106,label="PRY"),color="black",angle=90,size=10)+geom_text(aes(x=0.112,y=106,label="RBMY"),color="black",angle=90,size=10)+geom_text(aes(x=0.116,y=106,label="HSFY"),color="black",angle=90,size=10)+geom_text(aes(x=0.120,y=106,label="XKRY"),color="black",angle=90,size=10)+geom_text(aes(x=0.124,y=106,label="VCY"),color="black",angle=90,size=10)+geom_text(aes(x=0.128,y=106,label="TSPY"),color="black",angle=90,size=10)

#save figure
ggsave("../Figures/Fig_2.pdf",fig_2,height=15,width=15)
