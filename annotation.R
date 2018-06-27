# This code was used to analyze the SNP distribution on genome using GFF3
snp<-read.delim("/path/to/GBSclear.hmp.txt",header=T)
gff<-read.table("/path/to/Sbicolor_313_v3.1.gene.gff3",skip=2)
snp_inform<-snp[,1:11]
gff$V8<-NA
chr<-as.character(unique(gff$V1[grep("Chr",gff$V1)]))
features<-unique(gff$V3)[unique(gff$V3)!="mRNA"]
features<-as.character(features)
for(f in features){
SNP_GENE<-NULL
snp_number<-NULL
for (i in 1:10){
	snp_in_gene<-NULL
	gff_tmp<-gff[gff$V1==chr[i] &gff$V3==f,]
	snp_tmp<-snp_inform[snp_inform$chrom==i,]
	for (j in 1:dim(gff_tmp)[1]){
		tmp<-which(snp_tmp$pos%in%gff_tmp$V4[j]:gff_tmp$V5[j])
		gff_tmp$V8[j]<-length(tmp)
		snp_in_gene<-c(snp_in_gene,tmp)
	}
	snp_in_gene<-unique(snp_in_gene)
	snp_number_in_gene<-data.frame(chr=rep(i,length(snp_in_gene)),snp_in_gene)
	snp_number<-rbind(snp_number,snp_number_in_gene)
	cat("There are",dim(snp_number_in_gene)[1], "in ",f," for SNP from chromosome ",i,"\n")
	SNP_GENE<-rbind(SNP_GENE,gff_tmp)
	write.table(gff_tmp,paste0("snp_",f,"_chr",i),quote=F,row.names=F,col.names=F,sep="\t")
	if(i==10) cat("There are ",round(dim(snp_number)[1]/dim(snp)[1]*100,2)," percent of SNP in ",f,"! \n")
}
names(SNP_GENE) <- c('Chr','Source','feature','start','end','score','strand','snp_number','attributes')
write.table(SNP_GENE,paste0("SNP_",f,".txt"),quote=F,row.names=F,sep="\t")
}
