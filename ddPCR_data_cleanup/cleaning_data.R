
#load REQUIRED packages
require(plyr)
require(reshape2)
require(ggplot2)
require(ape)
require(ggtree)
require(ggrepel)
require(MASS)
require(nlme)

dat<-read.table('Table_S1A_ddPCR_outlier_analysis.txt',header=T,sep="\t")

#remove RPP30
#dat<-dat[-which(dat$Gene=="RPP30"),]

###removing one replicate observation for each sample####

#count number of nonmissing observations for each person

dat$n<-apply(dat[,c(3:5)],1,function(x){length(which(is.na(x)=="FALSE"))})

#calculate cv for each row before outlier removal
dat$sd1<-apply(dat[,c(3:5)],1,sd,na.rm=T)
dat$mean1<-apply(dat[,c(3:5)],1,mean,na.rm=T)
dat$cv1<-dat$sd1/dat$mean1
dat$median1<-apply(dat[,c(3:5)],1,median,na.rm=T)

#plot cv for each gene before outlier removal - Figure S1A
fig_s1a<-ggplot()+geom_boxplot(data=dat,aes(Gene,cv1),color="blue",outlier.shape=NA)+geom_point(data=dat,aes(Gene,cv1),position="jitter",alpha=0.7,color="blue")+theme_bw()+labs(x="Gene",y="Coefficient of Variation")+geom_hline(yintercept=median(dat$cv1,na.rm=T),color="red",linetype="dashed",size=1)+ylim(c(0,0.57))

ggsave('Fig_S1A.pdf',fig_s1a,height=7,width=7)

#determine which of the three replicates is the most distant from the other two. 
dat$outlier<-NA
for(i in 1:nrow(dat)){
  if(dat$n[i]>2){
    if(length(unique(as.character(dat[i,c('a','b','c')])))>1){
      g1<-dat[i,c('a','b','c')]-dat[i,'median1']
      distant.obs<-which(abs(g1)==max(abs(g1)))
      if(length(distant.obs)>1){dat$outlier[i]<-sample(distant.obs,1)}else{dat$outlier[i]<-distant.obs}
  }
  }
}

#create new dataframe and replace outlier values with NA

dat2<-dat[,c(1:6)]
for(i in 1:nrow(dat2)){
  if(is.na(dat$outlier[i])=="FALSE"){
    outlier.index<-dat$outlier[i]
    dat2[i,c('a','b','c')][outlier.index]<-NA
  }
}

#calculate mean, sd and coefficient of variation per gene per individual after outlier removal
dat2$n<-apply(dat2[,c(3:5)],1,function(x){length(which(is.na(x)=="FALSE"))})
dat2$mean<-apply(dat2[,c('a','b','c')],1,mean,na.rm=T)
dat2$sd<-apply(dat2[,c('a','b','c')],1,sd,na.rm=T)
dat2$cv<-dat2$sd/dat2$mean

#plot coefficient of variation for each gene and red line with median cv across all genes after outlier removal - Figure S1B
fig_s1b<-ggplot()+geom_boxplot(data=dat2,aes(Gene,cv),color="blue",outlier.shape=NA)+geom_point(data=dat2,aes(Gene,cv),position="jitter",alpha=0.7,color="blue")+theme_bw()+labs(x="Gene",y="Coefficient of Variation")+geom_hline(yintercept=median(dat2$cv,na.rm=T),color="red",linetype="dashed",size=1)

ggsave('Fig_S1B.pdf',fig_s1b,height=7,width=7)

#write cleaned data to file 
write.table('ddpcr_outliers_removed.txt',sep="\t",col.names=T,row.names=F,quote=F)
