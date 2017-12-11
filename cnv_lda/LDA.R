
#load cnv data
dat4<-read.table("../Data_files/ddpcr_outliers_removed.txt",header=T,sep="\t")
# run LDA with major haplogroup as category and gene counts as predictors
ld2<-lda(data=dat4,major_haplo~BPY+CDY+DAZ+HSFY+PRY+RBMY+TSPY+VCY+XKRY,CV=T,method="mle",prior=rep(1,10)/10)

#write function to plot posterior probabilities
plt.ld<-function(x,df){
  ldpx<-x$posterior
  ldpx<-as.data.frame(ldpx)
  ldpx$IID<-df$IID
  ldpx$major_haplo<-df$major_haplo
  haplos=c("R","Q","L","T","O","J","I","G","C","E")
  ldpx$pmatch<-NA
  ldpx$pmismatch<-NA
  for(i in 1:10){
    ldpx$pmatch[which(ldpx$major_haplo==haplos[i])]<-ldpx[which(ldpx$major_haplo==haplos[i]),haplos[i]]
    ldpx$pmismatch[which(ldpx$major_haplo==haplos[i])]<-apply(ldpx[which(ldpx$major_haplo==haplos[i]),haplos[-which(haplos==haplos[i])]],1,sum)
  }
  mldpx<-melt(ldpx[,c('IID','major_haplo','pmatch','pmismatch')],id.vars=c("IID",'major_haplo'))
  p<-ggplot(mldpx,aes(as.character(IID),value,fill=variable))+geom_bar(stat="identity",width=1)+facet_wrap(~major_haplo,scales="free_x",nrow=2)+theme_bw()+theme(axis.text.x=element_blank())
  p<-p+scale_fill_manual(values=c("#377eb8","#ff7f00"),labels=c("match","mismatch"))
  p<-p+labs(x="Individuals",y="Posterior Probability",fill="Match/Mismatch")
  return(p)
}

#plot!
pld2<-plt.ld(ld2,copy.number3)

ggsave("../Figures/Fig_5.pdf",pld2,width=7,height=5)
