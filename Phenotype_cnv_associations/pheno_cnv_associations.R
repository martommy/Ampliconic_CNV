require(plyr)
require(ape)
require(nlme)
require(ggtree)

#####Phenotype-CNV correlations#######

#load phenotype data
pheno<-read.table('../Data_files/phenotype_data.txt',header=T,sep="\t")

#load rooted 
ytree<-read.tree('../Data_files/y_rooted_tree.nwk')

#make tree ultrametric - i.e. calibrate tree based on timing of root = 71.7 kya (timing of split between E and the rest of the haplogrups - Karmin 2015)
mycalibration <- makeChronosCalib(ytree, node="root", age.max=71.7,age.min=71.7)
ytree <- chronos(ytree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )


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


#running one-way ANOVA (ultrametric ytree used to specify correlation among errors - Brownian(ie. lambda =1))
height.full.bpy<-gls(Height_cm~BPY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.cdy<-gls(Height_cm~CDY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.daz<-gls(Height_cm~DAZ,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.hsfy<-gls(Height_cm~HSFY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.pry<-gls(Height_cm~PRY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.rbmy<-gls(Height_cm~RBMY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.tspy<-gls(Height_cm~TSPY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.vcy<-gls(Height_cm~VCY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)
height.full.xkry<-gls(Height_cm~XKRY,correlation=corBrownian(phy=ytree),method="ML",data=pheno)


#running linear model between fmf and cnv
fmf.bpy<-lm(data=pheno,FMF~BPY,na.action="na.exclude")
fmf.cdy<-lm(data=pheno,FMF~CDY,na.action="na.exclude")
fmf.daz<-lm(data=pheno,FMF~DAZ,na.action="na.exclude")
fmf.hsfy<-lm(data=pheno,FMF~HSFY,na.action="na.exclude")
fmf.pry<-lm(data=pheno,FMF~PRY,na.action="na.exclude")
fmf.rbmy<-lm(data=pheno,FMF~RBMY,na.action="na.exclude")
fmf.tspy<-lm(data=pheno,FMF~TSPY,na.action="na.exclude")
fmf.vcy<-lm(data=pheno,FMF~VCY,na.action="na.exclude")
fmf.xkry<-lm(data=pheno,FMF~XKRY,na.action="na.exclude")


#running phylogenetic linear model between fmf and BPY. Run separately for each gene
#some missing data present in pheno. remove and run

missing_ids<-pheno[which(is.na(pheno$FMF)=="TRUE"),'IID']
ydat<-fortify(ytree)
nodes2drop<-ydat$node[which(ydat$label%in%missing_ids)]
ytree.red<-drop.tip(ytree,tip=nodes2drop)

pheno.red<-pheno[-which(pheno$IID%in%missing_ids),]

fmf.full.bpy<-gls(FMF~BPY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.cdy<-gls(FMF~CDY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.daz<-gls(FMF~DAZ,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.hsfy<-gls(FMF~HSFY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.pry<-gls(FMF~PRY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.rbmy<-gls(FMF~RBMY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.tspy<-gls(FMF~TSPY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.vcy<-gls(FMF~VCY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model
fmf.full.xkry<-gls(FMF~XKRY,correlation=corBrownian(phy=ytree.red),data=pheno.red,method="ML",na.action="na.exclude") #run full model



