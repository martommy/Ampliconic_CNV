require(ggplot2)
library(reshape2)


#plotting median copy number per gene vs variance
#read ampliconic gene numbers - cleaned
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#plot mean and variance of ampliconic genes
sum.dat4<-data.frame(Gene=colnames(dat4)[c(2:10)],Median=apply(dat4[,c(2:10)],2,median,na.rm=T),Variance=apply(dat4[,c(2:10)],2,var,na.rm=T))
fig_1<-ggplot(sum.dat4,aes(log(Median),log(Variance),color=Gene))+geom_point()+stat_smooth(method="lm",se=F,color="grey")+theme_bw()+geom_text_repel(aes(log(Median),log(Variance),label=Gene))
fig_1<-fig_1+theme(legend.position = "none")+labs(x="Log10(Median copy number)",y="Log10(Variance in copy number)")

#save
ggsave("../Figures/Fig_1.pdf",fig_1,height=7,width=7)



#plot copy number variation across genes and haplogroup

#melt dat4 to long-format
dat3<-melt(dat4,id.vars = c("IID","Clade","Haplogroup","geographic_loc","major_haplo"))
colnames(dat3)[c(6,7)]<-c("Gene","mean")

#plot copy number by gene and haplogroup
fig_3<-ggplot(dat3,aes(major_haplo,mean))+geom_boxplot(outlier.shape=NA,position=position_dodge(width=1))+geom_point(alpha=0.5,color="blue",position="jitter",size=0.7)+facet_wrap(~Gene,scale="free_y")+theme_bw()+theme()+labs(x="Y Haplogroup",y="Copy Number")

ggsave("../Figures/Fig_3.pdf",height=7,width=7)


#one-way Anova for every gene
#BPY
anova(lm(data=dat4,BPY~major_haplo))
#CDY
anova(lm(data=dat4,CDY~major_haplo))
#DAZ
anova(lm(data=dat4,DAZ~major_haplo))
#HSFY
anova(lm(data=dat4,HSFY~major_haplo))
#PRY
anova(lm(data=dat4,PRY~major_haplo))
#RBMY
anova(lm(data=dat4,RBMY~major_haplo))
#TSPY
anova(lm(data=dat4,TSPY~major_haplo))
#VCY
anova(lm(data=dat4,VCY~major_haplo))
#XKRY
anova(lm(data=dat4,XKRY~major_haplo))


###phylogenetic Anova###

##Eve requires three files - phylogenetic file, trait data file, number of individuals within each haplogroup

#load ytree
ytree<-read.tree('../Data_files/ytree.nwk')

#1. write phylogenetic file

##root phylogenetic tree
#C and E haplogroups are the oldest in our tree. get MRCA for any two descendants from them
getMRCA(ytree,c(95,74)) # ydat are node labels in ydat
# 173 is the parent node for these two lineages
rooted.ytree<-reroot(ytree,node=173) # this is not very important


#collapse all branches within each clade
col.ytree<-drop.tip(rooted.ytree,c(1:19,21:24,26:28,30:33,35:47,49:52,54:67,69:72,74:77,79:99))
#relabel the tips with their haplogroup
col.ytree$tip.label<-c("R","Q","L","T","O","J","I","G","C","E")

#make tree ultrametric - i.e. calibrate tree based on timing of root = 261.5 kya
mycalibration <- makeChronosCalib(col.ytree, node="root", age.max=261.5)
cal.col.ytree <- chronos(col.ytree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )

#write.tree to file
write.tree(cal.col.ytree,'y_timetree_haplo.nwk')

#add 10 as the first line - number of haplogroups -EVE requires this  
system('echo 10 | cat - y_timetree_haplo.nwk > y_timetree_haplo_eve.nwk')


#2. write trait/copy number data file

#create new dataframe to work with
copy.number<-dat4[,c(1:10,13)]

#remove row with missing data for CDY
copy.number2<-copy.number[which(is.na(copy.number$CDY)=="FALSE"),]

#order rows by haplogroups in the order of the phylogenetic tree file
copy.number2$order<-NA
copy.number2[which(copy.number2$major_haplo=="T"),'order']<-4
copy.number2[which(copy.number2$major_haplo=="L"),'order']<-3
copy.number2[which(copy.number2$major_haplo=="Q"),'order']<-2
copy.number2[which(copy.number2$major_haplo=="R"),'order']<-1
copy.number2[which(copy.number2$major_haplo=="E"),'order']<-10
copy.number2[which(copy.number2$major_haplo=="C"),'order']<-9
copy.number2[which(copy.number2$major_haplo=="G"),'order']<-8
copy.number2[which(copy.number2$major_haplo=="I"),'order']<-7
copy.number2[which(copy.number2$major_haplo=="J"),'order']<-6
copy.number2[which(copy.number2$major_haplo=="O"),'order']<-5
copy.number3<-copy.number2[order(copy.number2$order),]

#transpose so genes are rows and individuals are columns
tcopy.number3<-t(copy.number3[,-c(1,2,12)])

#write to file
write.table(tcopy.number3,'y.exprdat',row.names=T,col.names=F,quote=F,sep=" ")

#add 9 - number of genes to as first line - EVE requires this.
system('echo 9 | cat - y.exprdat > y_eve.exprdat')



