
require(ggplot2)
require(ggtree)
require(plyr)

#read ampliconic gene numbers - cleaned
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#read y haplogrup tree
ytree<-read.tree('../Data_files/Tree_segsites_IDs.nwk')

#convert tree to dataframe for easier plotting and labeling
ydat<-fortify(ytree)

#split label into individual id and haplogroup info
#ydat$IID<-sapply(ydat$label,function(x){unlist(strsplit(x,split="_"))[1]})
#ydat$haplogroup<-sapply(ydat$label,function(x){unlist(strsplit(x,split="_"))[2]})
#ydat$major_haplo<-sapply(ydat$haplogroup,function(x){unlist(strsplit(x,split=""))[1]})

#read haplogroup information
haplo<-read.table("../Data_files/haplogroup_info_11202017.txt",sep="\t",header=T)
ydat$IID<-ydat$label
ydat<-join(ydat,haplo,by="IID")

#add copy number data to phylogeny table and round off to 2 decimal places
ydat<-join(ydat,dat4,by="IID")

##root phylogenetic tree - for visualization and for carrying out phylogenetic ANOVA (EVE)


#place root at the base of E haplogroup i.e. African haplogroup
mrca.nodes<-ydat$node[which(ydat$major_haplo=="E")]
mrca.no<-getMRCA(ytree,mrca.nodes) # ydat are node labels in ydat
rooted.ytree<-reroot(ytree,node=mrca.no) # this is not very important

#data corresponding to rooted phlogeny
ydat<-fortify(rooted.ytree)
ydat$IID<-ydat$label
ydat<-join(ydat,haplo,by="IID")

#add copy number data to phylogeny table and round off to 2 decimal places
ydat<-join(ydat,dat4,by="IID")
#plot y phylogeny with haplogroups colored
groupInfo<-split(ydat$label,ydat$major_haplo)
rooted.ytree<-groupOTU(rooted.ytree,groupInfo)

p<-ggtree(rooted.ytree,aes(color=group))+
  geom_segment(data=ydat[which(is.na(ydat$Haplogroup)=="FALSE"),],aes(x=x,xend=0.50,y=y,yend=y,color=major_haplo),linetype="dotted")

#add scale
p<-p+geom_treescale(y=-2,offset=-2)

ydat[,c('BPY','CDY','DAZ','HSFY','PRY','RBMY','TSPY','VCY','XKRY')]<-apply(ydat[,c('BPY','CDY','DAZ','HSFY','PRY','RBMY','TSPY','VCY','XKRY')],2,function(x){format(round(x,2),nsmall=2)})
ydat<-ydat[which(is.na(ydat$major_haplo)=="FALSE"),]

#add copy number information - Figure 2
#create setp and base to set placement and width of number columns
base<-max(ydat$x+0.01)
step<-0.03
fig_2<-p+geom_text(data=ydat,aes(x=base,label=CDY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+step,label=BPY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+2*step,label=DAZ,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+3*step,label=PRY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+4*step,label=RBMY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+5*step,label=HSFY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+6*step,label=XKRY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+7*step,label=VCY,y=y,color=major_haplo))+
  geom_text(data=ydat,aes(x=base+8*step,label=TSPY,y=y,color=major_haplo))+
  geom_text(aes(x=base,y=106,label="CDY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+step,y=106,label="BPY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+2*step,y=106,label="DAZ"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+3*step,y=106,label="PRY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+4*step,y=106,label="RBMY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+5*step,y=106,label="HSFY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+6*step,y=106,label="XKRY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+7*step,y=106,label="VCY"),color="black",angle=90,size=10)+
  geom_text(aes(x=base+8*step,y=106,label="TSPY"),color="black",angle=90,size=10)

#save figure
ggsave("../Figures/Fig_2.pdf",fig_2,height=15,width=15)


#write rooted tree to file
write.tree(rooted.ytree,'../Data_files//y_timetree_haplo.nwk')



