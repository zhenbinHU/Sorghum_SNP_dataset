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


# LD heatmmap plot using LDheatmap R package 
# Using Y1 gene as example
# Color from blue to red mean R^2 from 0 to 1

library(LDheatmap)
setwd("~/Desktop/heatmap")
load("/Users/zhenbin/Desktop/heatmap/snp_dis_y1.1.rda")
load("/Users/zhenbin/Desktop/heatmap/snp_data_y1.rda")
png("heatmap_Y1.png",width=6,height=6,res=600,units="in")
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
MyHeatmap<-LDheatmap(snp_data,snp_dis,flip=TRUE,color=rgb.palette(20),title="",name = "myLDgrob",SNP.name="S2_57610965")
LD.grob1 <- editGrob(MyHeatmap$LDheatmapGrob, gPath("heatMap", "title"), gp = gpar(cex = 1.25, col = "blue"))
LD.grob2 <- editGrob(LD.grob1, gPath("Key","title"), gp = gpar(cex = 1.5, col = "black"))
LD.grob3 <- editGrob(LD.grob2, gPath("Key", "labels"),gp = gpar(cex = 1.3, col = "black"))
grid.newpage() 
grid.draw(LD.grob3)
dev.off()

