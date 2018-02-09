require(ggplot2)
require(reshape2)
require(ggrepel)
require(plyr)
require(ape)
require(ggtree)
#plotting median copy number per gene vs variance
#read ampliconic gene numbers - cleaned
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#plot mean and variance of ampliconic genes
sum.dat4<-data.frame(Gene=colnames(dat4)[c(2:10)],Median=apply(dat4[,c(2:10)],2,median,na.rm=T),Variance=apply(dat4[,c(2:10)],2,var,na.rm=T),sd=apply(dat4[,c(2:10)],2,sd,na.rm=T),min=apply(dat4[,c(2:10)],2,function(x){return(range(x,na.rm=T)[1])}),max=apply(dat4[,c(2:10)],2,function(x){range(x,na.rm=T)[2]}))
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

#load rooted ytree
rooted.ytree<-read.tree('../Data_files/y_rooted_tree.nwk')

#get haplogroup and clade info
haplo<-read.table("../Data_files/haplogroup_info_11202017.txt",sep="\t",header=T)
ydat<-fortify(rooted.ytree)
ydat$IID<-ydat$label
ydat<-join(ydat,haplo,by="IID")

#collapse all branches within each clade
#for each haplogroup, list the range of node numbers
#write function for this
desc.nodes<-function(x){
  #where x is a string - major haplogroup
  nodes.range<-ydat$node[grep(x,ydat$major_haplo)]
  nodes.no.x<-seq(min(nodes.range),max(nodes.range)-1)
  return(nodes.no.x)
}

col.ytree<-drop.tip(rooted.ytree,c(desc.nodes("C"),
                                   desc.nodes("G"),
                                   desc.nodes("E"),
                                   desc.nodes("R"),
                                   desc.nodes("Q"),
                                   desc.nodes("L"),
                                   desc.nodes("T"),
                                   desc.nodes("O"),
                                   desc.nodes("J"),
                                   desc.nodes("I")
))
#relabel the tips with their haplogroup
col.ytree$tip.label<-as.character(ydat$major_haplo[which(ydat$label%in%col.ytree$tip.label)])

#make tree ultrametric - i.e. calibrate tree based on timing of root = 71.7 kya (timing of split between E and the rest of the haplogrups - Karmin 2015)
mycalibration <- makeChronosCalib(col.ytree, node="root", age.max=71.7,age.min=71.7)
cal.col.ytree <- chronos(col.ytree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )

#replace tip labels with integers for EVE
cal.col.ytree.int<-cal.col.ytree
cal.col.ytree.int$tip.label<-as.character(seq(1,10,1))

#write.tree to file
write.tree(cal.col.ytree.int,'y_timetree_ultrametric_haplo.nwk')

#add 10 as the first line - number of haplogroups -EVE requires this  
system('echo 10 | cat - y_timetree_ultrametric_haplo.nwk > y_timetree_haplo_eve.nwk')


#2. write trait/copy number data file

#create new dataframe to work with
copy.number<-dat4[,c(1:10,13)]

#remove row with missing data for CDY
copy.number2<-copy.number[which(is.na(copy.number$CDY)=="FALSE"),]

#order rows by haplogroups in the order of the phylogenetic tree file
tip.labs<-data.frame(major_haplo=cal.col.ytree$tip.label,order=seq(1,length(cal.col.ytree$tip.label)))
copy.number2<-join(copy.number2,tip.labs,by="major_haplo")

copy.number3<-copy.number2[order(copy.number2$order),]


#transpose so genes are rows and individuals are columns
tcopy.number3<-t(copy.number3[,-c(1,11,12)])

#write to file
write.table(tcopy.number3,'y.exprdat',row.names=T,col.names=F,quote=F,sep=" ")

#add 9 - number of genes to as first line - EVE requires this.
system('echo 9 | cat - y.exprdat > y_eve.exprdat')


#3. write individual file - number of individuals per haplogroup in the SAME order as phylogenetic tree

#order levels of major_haplogroup to match phylogeny
copy.number3$major_haplo<-factor(copy.number3$major_haplo,levels=c("R","Q","L","T","O","J","I","G","C","E"))

#write frequency of observations per haplogroup in dataframe to file
write.table(t(table(copy.number3$major_haplo)),'y_eve.nindv',col.names=F,row.names=F,quote=F)

##run EVE in terminal
#assuming EVE binary and input data are in current directory 
#Read README file for more detail
#./EVEmodel -S -n 12 -t y_timetree_haplo_eve.nwk -i y_eve.nindv -d y_eve.exprdat -f _trialRun -v 10


