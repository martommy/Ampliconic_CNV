#run pca on SNP genotype data with plink
#system('plink --bfile ched_polysnps2 --pca 300 --out pca_11182017')
require(ggplot2)
require(plyr)

#load eigenvector and eigenvalue files for SNPs into R
eigval.snps<-read.table('../Data_files/pca_snps.eigenval',header=F,sep="\t")
eigvec.snps<-read.table('../Data_files/pca_snps.eigenvec',header=F,sep="\t")
colnames(eigvec.snps)<-c("FID","IID",paste("PC",seq(1,length(eigvec.snps)-2,1),sep=""))

#plot screeplot/proportion of variance explained by each PC - snps - part of fig_s2
eigval.snps$prop<-eigval.snps$V1/sum(eigval.snps$V1)
eigval.snps$PC<-seq(1,nrow(eigval.snps),1)
screesnps<-ggplot(eigval.snps,aes(PC,prop))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(y="Prop. of variance explained")+
  scale_x_continuous(breaks=seq(1,9,1),limits=c(1,9))+ # limit no. of PCs in plot to 9
  scale_y_continuous(limits=c(0,0.31)) # set upper limit to max eigenvalue from SNP data
ggsave("../Figures/Fig_S2A.pdf",screesnps,height=5,width=7)

#load haplogroup info
haplo<-read.table('../Data_files/haplogroup_info_11202017.txt',sep="\t",header=T)


#add haplogroup info to PCs for plotting
eigvec.snps<-join(eigvec.snps,haplo,by="IID")

#plot PC1 v PC2 for SNPs - part of fig_4 - combine in illustrator
snps.p1vp2<-ggplot(eigvec.snps,aes(PC1,PC2,color=major_haplo))+geom_point(alpha=0.7)+theme_bw()+labs(color="Haplogroup")
ggsave("../Figures/Fig_4A1.pdf",height=7,width=7)
#plot PC1 v PC3 for SNPs - part of fig_4 - combine in illustrator
snps.p1vp3<-ggplot(eigvec.snps,aes(PC1,PC3,color=major_haplo))+geom_point(alpha=0.7)+theme_bw()+labs(color="Haplogroup")
ggsave("../Figures/Fig_4A2.pdf",height=7,width=7)

#load cnv data
dat4<-read.table('../Data_files/ddpcr_outliers_removed.txt',header=T,sep="\t")

#run pca on ampliconic gene copy number - remove rows with missing values first
dat.pca<-dat4[which(is.na(dat4$CDY)=="FALSE"),]
pca<-prcomp(dat.pca[,c(2:10)],center=T,scale=T)

#load eigenvalues and eigenvectors for ampliconic genes
eigval.cnv<-data.frame(V1=pca$sdev)
eigvec.cnv<-data.frame(pca$x)
eigvec.cnv$IID<-as.character(dat.pca$IID)
eigvec.cnv$major_haplo<-as.character(dat.pca$major_haplo)

#plot screeplot/proportion of variance explained by each PC - cnv - part of fig_s2
eigval.cnv$prop<-eigval.cnv$V1/sum(eigval.cnv$V1)
eigval.cnv$PC<-seq(1,nrow(eigval.cnv),1)
screecnv<-ggplot(eigval.cnv,aes(PC,prop))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(y="Prop. of variance explained")+
  scale_x_continuous(breaks=seq(1,9,1),limits=c(1,9))+ # limit no. of PCs in plot to 10
  scale_y_continuous(limits=c(0,0.31)) # set upper limit to max eigenvalue of SNP data
ggsave("../Figures/Fig_S2B.pdf",screecnv,height=5,width=7)

#plot PC1 v PC2 for CNVs - part of fig_4 - combine in illustrator
cnv.p1vp2<-ggplot(eigvec.cnv,aes(PC1,PC2,color=major_haplo))+geom_point(alpha=0.7)+theme_bw()+labs(color="Haplogroup")
#plot PC1 v PC3 for CNVs - part of fig_4 - combine in illustrator
cnv.p1vp3<-ggplot(eigvec.cnv,aes(PC1,PC3,color=major_haplo))+geom_point(alpha=0.7)+theme_bw()+labs(color="Haplogroup")

ggsave("../Figures/Fig_4B1.pdf",cnv.p1vp2,height=7,width=7)
ggsave("../Figures/Fig_4B2.pdf",cnv.p1vp3,height=7,width=7)

