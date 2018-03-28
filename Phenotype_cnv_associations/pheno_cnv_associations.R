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


genes<-c("BPY","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")

#running one-way ANOVA (no phylogenetic correction) between height and FMF separately for each gene

#formula to loop over genes
simple.lm<-function(phenotype,gene){
  form1<-as.formula(paste(phenotype,'~',gene,sep=""))
  l1<-lm(data=pheno,form1,na.action="na.exclude")
  beta1<-summary(l1)$coefficients[2,1]
  tstat1<-summary(l1)$coefficients[2,3]
  pvalue1<-summary(l1)$coefficients[2,4]
  return(data.frame(beta=beta1,tstat=tstat1,pvalue=pvalue1))
}

height.lm<-t(sapply(genes,function(x){simple.lm("Height_cm",x)}))
fmf.lm<-t(sapply(genes,function(x){simple.lm("FMF",x)}))


#running PGLS between Height/FMF and gene copy number while simultaneously estimating Pagel's lambda

#write function to do this for height
gls.pagel.height<-function(gene){
  form1<-as.formula(paste('Height_cm~',gene,sep=""))
  gls1<-gls(form1,correlation=corPagel(phy=ytree,value=1,fixed=F),method="ML",data=pheno)
  lambda1<-coef(gls1$modelStruct)
  beta1<-summary(gls1)$tTable[2,1]
  tstat1<-summary(gls1)$tTable[2,3]
  pvalue1<-summary(gls1)$tTable[2,4]
  return(data.frame(beta=beta1,tstat=tstat1,pvalue=pvalue1,lambda=lambda1))
}

#write function to do this for FMF
#some missing data present in FMF. remove and create new data frame

missing_ids<-pheno[which(is.na(pheno$FMF)=="TRUE"),'IID']
ydat<-fortify(ytree)
nodes2drop<-ydat$node[which(ydat$label%in%missing_ids)]
ytree.red<-drop.tip(ytree,tip=nodes2drop)

pheno.red<-pheno[-which(pheno$IID%in%missing_ids),]

gls.pagel.fmf<-function(gene){
  form1<-as.formula(paste('FMF~',gene,sep=""))
  gls1<-gls(form1,correlation=corPagel(phy=ytree.red,value=1,fixed=F),method="ML",data=pheno.red)
  lambda1<-coef(gls1$modelStruct)
  beta1<-summary(gls1)$tTable[2,1]
  tstat1<-summary(gls1)$tTable[2,3]
  pvalue1<-summary(gls1)$tTable[2,4]
  return(data.frame(beta=beta1,tstat=tstat1,pvalue=pvalue1,lambda=lambda1))
}

#PGLS for height - lambda also included in output
height.pagel<-t(sapply(genes,gls.pagel.height))

#PGLS for FMF - lambda also included in output
fmf.pagel<-t(sapply(genes,gls.pagel.fmf))

write.table(height.lm,'Height_simple_lmresults.txt',sep="\t",col.names=T,row.names=T,quote=F)
write.table(fmf.lm,"FMF_simple_lmresults.txt",sep="\t",col.names=T,row.names=T,quote=F)
write.table(height.pagel,'Height_PGLS_results.txt',sep="\t",col.names=T,row.names = T,quote=F)
write.table(fmf.pagel,'FMF_PGLS_results.txt',sep="\t",col.names=T,row.names=T,quote=F)

