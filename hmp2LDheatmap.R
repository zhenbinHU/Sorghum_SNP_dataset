# This function was used to convert hapmap format to the format LDheatmap package need
# The input is a subsetted hapmap by the region we interested, such as around 50 kb of y1 gene
# out is a names of the file we want to save as for our output

hmp2LDheatmap<-function(hapmap,out=""){
    hapmap[12:dim(hapmap)[2]]<-apply(hapmap[12:dim(hapmap)[2]],2,as.character)
    snp_dis<-as.numeric(as.character(hapmap$pos))
    hapmap[hapmap=="A"]<-"A/A"
    hapmap[hapmap=="T"]<-"T/T"
    hapmap[hapmap=="C"]<-"C/C"
    hapmap[hapmap=="G"]<-"G/G"
    hapmap[hapmap=="M"]<-"A/C"
    hapmap[hapmap=="R"]<-"A/G"
    hapmap[hapmap=="W"]<-"A/T"
    hapmap[hapmap=="S"]<-"C/G"
    hapmap[hapmap=="Y"]<-"C/T"
    hapmap[hapmap=="K"]<-"G/T"
    library(genetics)
    snp_infor<-hapmap[,1:11]
    snp_data<-hapmap[,12:dim(hapmap)[2]]
    snp_data<-t(snp_data)
    snp_data<-as.data.frame(snp_data)
    names(snp_data)<-as.character(snp_infor$rs.)
    num<-ncol(snp_data)
    for(i in 1:num){
        snp_data[,i]<-as.genotype(snp_data[,i]) # convert the columns of the original data frame into  genotype objects
    }

save(snp_data,file=paste0("snp_data_",out,".rda"))
save(snp_dis,file=paste0("snp_dis_",out,".rda"))
}

# LD heatmmap plot using LDheatmap R package 
# Using Y1 gene as example
# Color from blue to red mean R^2 from 0 to 1
hmp2LDheatmap(y1_hapmap,out="y1") # y1_hapmap is the regional SNP data based on y1 gene.

library(LDheatmap)
load("snp_dis_y1.1.rda")
load("snp_data_y1.rda")
png("heatmap_Y1.png",width=6,height=6,res=600,units="in",type="cairo")
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
MyHeatmap<-LDheatmap(snp_data,snp_dis,flip=TRUE,color=rgb.palette(20),title="",name = "myLDgrob",SNP.name="S2_57610965")
LD.grob1 <- editGrob(MyHeatmap$LDheatmapGrob, gPath("heatMap", "title"), gp = gpar(cex = 1.25, col = "blue"))
LD.grob2 <- editGrob(LD.grob1, gPath("Key","title"), gp = gpar(cex = 1.5, col = "black"))
LD.grob3 <- editGrob(LD.grob2, gPath("Key", "labels"),gp = gpar(cex = 1.3, col = "black"))
grid.newpage() 
grid.draw(LD.grob3)
dev.off()
