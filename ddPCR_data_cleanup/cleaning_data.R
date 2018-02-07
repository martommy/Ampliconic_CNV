
#load REQUIRED packages
require(plyr)
require(reshape2)
require(ggplot2)
require(ape)
require(ggtree)
require(ggrepel)
require(MASS)
require(nlme)

dat<-read.table('../Data_files/Table_S1A_ddPCR_outlier_analysis.txt',header=T,sep="\t")

#remove RPP30
#dat<-dat[-which(dat$Gene=="RPP30"),]

###removing one replicate observation for each sample####

#count number of nonmissing observations for each person

dat$n<-apply(dat[,c(3:8)],1,function(x){length(which(is.na(x)=="FALSE"))})

#calculate cv for each row before outlier removal
dat$sd1<-apply(dat[,c(3:8)],1,sd,na.rm=T)
dat$mean1<-apply(dat[,c(3:8)],1,mean,na.rm=T)
dat$cv1<-dat$sd1/dat$mean1
dat$median1<-apply(dat[,c(3:8)],1,median,na.rm=T)

#plot cv for each gene before outlier removal - Figure S1A
fig_s1a<-ggplot()+geom_boxplot(data=dat,aes(Gene,cv1),color="blue",outlier.shape=NA)+geom_point(data=dat,aes(Gene,cv1),position="jitter",alpha=0.7,color="blue")+theme_bw()+labs(x="Gene",y="Coefficient of Variation")+geom_hline(yintercept=median(dat$cv1,na.rm=T),color="red",linetype="dashed",size=1)+ylim(c(0,0.57))

ggsave('../Figures/Fig_S1A.pdf',fig_s1a,height=7,width=7)

#determine which of the three replicates is the most distant from the other two. 
dat$outlier<-NA
for(i in 1:nrow(dat)){
  if(dat$n[i]>2){
    if(length(unique(as.character(dat[i,c('a','b','c','d','e','f')])))>1){
      g1<-dat[i,c('a','b','c','d','e','f')]-dat[i,'median1']
      distant.obs<-which(abs(g1)==max(abs(g1),na.rm=T))
      if(length(distant.obs)>1){dat$outlier[i]<-sample(distant.obs,1)}else{dat$outlier[i]<-distant.obs}
  }
  }
}

#create new dataframe and replace outlier values with NA

dat2<-dat
for(i in 1:nrow(dat2)){
  if(is.na(dat$outlier[i])=="FALSE"){
    outlier.index<-dat$outlier[i]
    dat2[i,c('a','b','c','d','e','f')][outlier.index]<-NA
  }
}

#calculate mean, sd and coefficient of variation per gene per individual after outlier removal
dat2$n<-apply(dat2[,c(3:5)],1,function(x){length(which(is.na(x)=="FALSE"))})
dat2$mean<-apply(dat2[,c('a','b','c','d','e','f')],1,mean,na.rm=T)
dat2$sd<-apply(dat2[,c('a','b','c','d','e','f')],1,sd,na.rm=T)
dat2$cv<-dat2$sd/dat2$mean

#plot coefficient of variation for each gene and red line with median cv across all genes after outlier removal - Figure S1B
fig_s1b<-ggplot()+geom_boxplot(data=dat2,aes(Gene,cv),color="blue",outlier.shape=NA)+geom_point(data=dat2,aes(Gene,cv),position="jitter",alpha=0.7,color="blue")+theme_bw()+labs(x="Gene",y="Coefficient of Variation")+geom_hline(yintercept=median(dat2$cv,na.rm=T),color="red",linetype="dashed",size=1)

ggsave('../Figures/Fig_S1B.pdf',fig_s1b,height=7,width=7)

#add haplogroup information to copy number data
#read table with haplogroup information for each ID
haplo<-read.table('../Data_files/haplogroup_info_11202017.txt',sep="\t",header=T)

#merge copy number with haplogroup info in one dataframe
dat3<-merge(dat2,haplo,by="IID",sort=F)

#pivot table so that ids are rows and columns are genes. average copy number across replicates
dat4<-dcast(dat3,IID~Gene,value.var="mean")
dat4$IID<-as.character(dat4$IID)
dat4<-join(dat4,haplo,by="IID")

#write cleaned data to file 
write.table(dat4,'../Data_files/ddpcr_outliers_removed.txt',sep="\t",col.names=T,row.names=F,quote=F)

