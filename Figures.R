# manhattan plot was draw using qqman packages
setwd("")
list<-list.files()
for(i in list){
   setwd(paste0(i))
   tmp<-list.files(pattern="GWAS.Results.csv")
   tmp<-read.csv(tmp,header=T)
   nobs=unique(tmp$nobs)
   tmp<-tmp[tmp$maf>0.05,1:5]
   names(tmp)<-c("SNP","CHR","BP","P","maf")
   png(paste0(i,".png"),width=12,height=6,type="cairo",res=600,units="in")
       par(mai=c(1,1.2,0.5,0.5))
       manhattan(tmp,suggestiveline = F, genomewideline = -log10(0.05/dim(tmp)[1]),ylim=c(0,max(-log10(tmp$P))+3), cex.lab=2,cex.ax
is=1.5)
       dev.off()
   print(paste0("Pop size for trait:", i," is: ",nobs))
 }


