# This function was used to convert hapmap format to the format LDheatmap package need
# The input is a subsetted hapmap which is the region we interested
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
