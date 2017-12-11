#calculate median 
require(ggplot2)
require(plyr)

#read ampliconic gene numbers - cleaned
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#plot mean and variance of ampliconic genes
sum.dat4<-data.frame(Gene=colnames(dat4)[c(2:10)],Median=apply(dat4[,c(2:10)],2,median,na.rm=T),Variance=apply(dat4[,c(2:10)],2,var,na.rm=T))
fig_1<-ggplot(sum.dat4,aes(log(Median),log(Variance),color=Gene))+geom_point()+stat_smooth(method="lm",se=F,color="grey")+theme_bw()+geom_text_repel(aes(log(Median),log(Variance),label=Gene))
fig_1<-fig_1+theme(legend.position = "none")+labs(x="Log10(Median copy number)",y="Log10(Variance in copy number)")

#save
ggsave("../Figures/Fig_1.pdf",fig_1,height=7,width=7)


