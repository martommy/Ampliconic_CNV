require(plyr)
require(ape)
require(nlme)
#####Phenotype-CNV correlations#######

#load phenotype data
pheno<-read.table('../Data_files/phenotype_data.txt',header=T,sep="\t")

#load ytree
ytree<-read.tree('../Data_files/ytree.nwk')

#load cnv data
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

ytree$tip.label<-sapply(ytree$tip.label,function(x){unlist(strsplit(x,split="_"))[1]})
pheno<-pheno[which(pheno$IID%in%ytree$tip.label),]
pheno<-join(pheno,dat4,by="IID")
row.names(pheno)<-as.character(pheno$IID)

#running one-way ANOVA (no phylogenetic correction between height and BPY). Run separately for each gene
height.bpy<-lm(data=pheno,Height_cm~BPY,na.action="na.exclude")
height.cdy<-lm(data=pheno,Height_cm~CDY,na.action="na.exclude")
height.daz<-lm(data=pheno,Height_cm~DAZ,na.action="na.exclude")
height.hsfy<-lm(data=pheno,Height_cm~HSFY,na.action="na.exclude")
height.pry<-lm(data=pheno,Height_cm~PRY,na.action="na.exclude")
height.rbmy<-lm(data=pheno,Height_cm~RBMY,na.action="na.exclude")
height.tspy<-lm(data=pheno,Height_cm~TSPY,na.action="na.exclude")
height.vcy<-lm(data=pheno,Height_cm~VCY,na.action="na.exclude")
height.xkry<-lm(data=pheno,Height_cm~XKRY,na.action="na.exclude")

#running phylogenetic linear model between height and BPY. Run separately for each gene
height.phylo<-gls(Height_cm~BPY,correlation=corBrownian(phy=ytree),data=pheno,method="ML",na.action="na.exclude") #run full model

#running phylogenetic linear model between fmf and BPY. Run separately for each gene
#some missing data present in pheno. remove and run

missing_ids<-pheno[which(is.na(pheno$FMF)=="TRUE"),'IID']
ydat<-fortify(ytree)
nodes2drop<-ydat$node[which(ydat$label%in%missing_ids)]
ytree.red<-drop.tip(ytree,tip=nodes2drop)

pheno.red<-pheno[-which(pheno$IID%in%missing_ids),]

fmf.phylo<-gls(FMF~BPY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model



