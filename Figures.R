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

# code for Figure 1
gff<-read.table("/path/Sbicolor_313_v3.1.gene.gff3",skip=1)
gff<-gff[grep("Chr",gff$V1),]
features<-unique(gff$V3)
space<-gff$V5-gff$V4

genome_length<-sum(c(80884392,77742459,74386277,68658214,71854669,61277060,65505356,62686529,59416394,61233695))

gene_length<-sum(space[gff$V3=="gene"])/genome_length
snp_ann<-c(100-47.86,47.86,47.86-27.08-5.48-3.9,27.08,5.48,3.9)
mRNA_length<-sum(space[grep("gene",gff$V3)+1])/genome_length
CDS_length<-sum(space[gff$V3=="CDS"])/genome_length
UTR5<-sum(space[gff$V3=="five_prime_UTR"])/genome_length
UTR3<-sum(space[gff$V3=="five_prime_UTR"])/genome_length
intron_length<-gene_length-CDS_length-UTR5-UTR3
intron_length<-gene_length-CDS_length-UTR5-UTR3
dat<-c(1-gene_length,gene_length,intron_length,CDS_length,UTR5,UTR3)
dat<-data.frame(snp_ann,dat)
dat$dat<-dat$dat*100

setwd("/path/Figure_1")
snp_dens<-read.table("snp_density_200k.snpden",header=T)
sim_dig<-read.table("density_cut_site.txt",header=T)

snp_dens$End<-snp_dens$BIN_START+200000
snp_dens$BIN_START<-snp_dens$BIN_START+1

sim_dig$CHR<-gsub("CHR","",sim_dig$CHR)
sim_dig$CHR<-as.numeric(sim_dig$CHR)

chr1<-nrow(sim_dig[sim_dig$CHR==1,])/2
chr2<-nrow(sim_dig[sim_dig$CHR<2,])+nrow(sim_dig[sim_dig$CHR==2,])/2
chr3<-nrow(sim_dig[sim_dig$CHR<3,])+nrow(sim_dig[sim_dig$CHR==3,])/2
chr4<-nrow(sim_dig[sim_dig$CHR<4,])+nrow(sim_dig[sim_dig$CHR==4,])/2
chr5<-nrow(sim_dig[sim_dig$CHR<5,])+nrow(sim_dig[sim_dig$CHR==5,])/2
chr6<-nrow(sim_dig[sim_dig$CHR<6,])+nrow(sim_dig[sim_dig$CHR==6,])/2
chr7<-nrow(sim_dig[sim_dig$CHR<7,])+nrow(sim_dig[sim_dig$CHR==7,])/2
chr8<-nrow(sim_dig[sim_dig$CHR<8,])+nrow(sim_dig[sim_dig$CHR==8,])/2
chr9<-nrow(sim_dig[sim_dig$CHR<9,])+nrow(sim_dig[sim_dig$CHR==9,])/2
chr10<-nrow(sim_dig[sim_dig$CHR<10,])+nrow(sim_dig[sim_dig$CHR==10,])/2

CHR1<-nrow(snp_dens[snp_dens$CHROM==1,])/2
CHR2<-nrow(snp_dens[snp_dens$CHROM<2,])+nrow(snp_dens[snp_dens$CHROM==2,])/2
CHR3<-nrow(snp_dens[snp_dens$CHROM<3,])+nrow(snp_dens[snp_dens$CHROM==3,])/2
CHR4<-nrow(snp_dens[snp_dens$CHROM<4,])+nrow(snp_dens[snp_dens$CHROM==4,])/2
CHR5<-nrow(snp_dens[snp_dens$CHROM<5,])+nrow(snp_dens[snp_dens$CHROM==5,])/2
CHR6<-nrow(snp_dens[snp_dens$CHROM<6,])+nrow(snp_dens[snp_dens$CHROM==6,])/2
CHR7<-nrow(snp_dens[snp_dens$CHROM<7,])+nrow(snp_dens[snp_dens$CHROM==7,])/2
CHR8<-nrow(snp_dens[snp_dens$CHROM<8,])+nrow(snp_dens[snp_dens$CHROM==8,])/2
CHR9<-nrow(snp_dens[snp_dens$CHROM<9,])+nrow(snp_dens[snp_dens$CHROM==9,])/2
CHR10<-nrow(snp_dens[snp_dens$CHROM<10,])+nrow(snp_dens[snp_dens$CHROM==10,])/2

#cor_snp_sim
out<-"path"
sim_dig$id<-paste0(sim_dig$CHR,"_",sim_dig$Start)
snp_dens$id<-paste0(snp_dens$CHR,"_",snp_dens$BIN_START)
cor_snp_sim<-merge(sim_dig,snp_dens,by.x="id",by.y="id")

library(beeswarm)

sim_dig$id<-paste0(sim_dig$CHR,"_",sim_dig$Start)
snp_dens$id<-paste0(snp_dens$CHR,"_",snp_dens$BIN_START)
cor_snp_sim<-merge(sim_dig,snp_dens,by.x="id",by.y="id")

snp_number<-c(61320,56334,56927,48197,46334,41295,36086,35112,36460,41224)
chr_length<-c(80884392,77742459,74386277,68658214,71854669,61277060,65505356,62686529,59416394,61233695)

individual_bef<-read.delim("Bef_merged_missing.imiss",header=T)
individual_merged<-read.delim("merged_missing.imiss",header=T)

png("Figure_1.png",width=10,height=8,res=600,units="in",type="cairo")
layout(matrix(c(1,1,2,3,4,4), 3, 2, byrow = TRUE))
par(mar= c(5, 5, 1.5, 2))
# digest

plot(snp_dens$SNP_COUNT,col=ifelse(snp_dens$CHROM%%2==0,"black","gray"),type="h",xaxt = "n",xlab="Chromosome",ylab="SNP number",cex.lab=1.5,axes=F)
axis(side=2,at=seq(0,600,length=7),labels=seq(0,600,length=7),cex.lab=1.5,cex.axis = 1.3,las=2)
axis(side=1,at=c(CHR1,CHR2,CHR3,CHR4,CHR5,CHR6,CHR7,CHR8,CHR9,CHR10),labels=c("1","2","3","4","5","6","7","8","9","10"),cex.axis=2)

sample<-individual_bef$INDV
sample<-gsub("C57CTACXX.+","",sample)
sample<-gsub("[:]$","",sample)
individual_bef$INDV<-sample
sample<-as.data.frame(table(sample))
sample<-sample[sample$Freq!=1,]
missing_bef<-individual_bef[individual_bef$INDV%in%sample$sample,]
missing_merged<-individual_merged[individual_merged$INDV%in%sample$sample,]
missing_bef<-missing_bef[-c(grep("Blank",missing_bef$INDV)),]
missing_bef$type<-1
missing_merged$type<-2
missing<-rbind(missing_bef,missing_merged)

beeswarm(F_MISS ~ type, data = missing,cex = 0.25, pch = 16,pwcol = as.numeric(type), xlab = "n",ylab = "n",cex.lab=5,cex.axis=2,axes = FALSE, ann = FALSE)
boxplot(F_MISS ~ type, data = missing,cex.lab=5,cex.axis=1, add = T,names = c("",""),col="#0000ff22",axes = FALSE, ann = FALSE)
axis(1, at = 1:2, labels = c("Before", "After"), cex.axis = 1.5)
axis(2, at=seq(0,1,length=6),labels=seq(0,1,length=6)*100,cex.axis = 1.5,las=2)
title(ylab = "Rate of missing (%)", cex.lab =1.5,line = 3)

snp_ann<-c(100-47.86,47.86,47.86-27.08-5.48-3.9,27.08,5.48,3.9)
nam<-c("Intergenic","Intragenic","Intron","CDS","5'-UTR","3'-UTR")
barplot(t(dat),beside=TRUE,las=2,col=c("deepskyblue","blue"),border = F, ylab="Percentage (%)",cex.axis=1.5,cex.lab=1.5)
text(seq(2,17,length=6),rep(-1,times=6),labels=nam,srt=45,cex=1.5,adj = 1,xpd = TRUE)
legend("topright",legend=c("SNP distribution","Genomic features"),col=c("deepskyblue","blue"),pch=15,bty="n",cex=2)

plot(sim_dig$number,col=ifelse(sim_dig$CHR%%2==0,"black","gray"),type="h",xaxt = "n",xlab="Chromosome",ylab="Enzyme cut sites", cex.lab=1.5,axes=F)
axis(side=2,at=seq(0,600,length=7),labels=seq(0,600,length=7),cex.lab=1.5,cex.axis = 1.3,las=2)
axis(side=1,at=c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10),labels=c("1","2","3","4","5","6","7","8","9","10"),cex.axis=2)
dev.off()

# Figure S4, the correlation between SNP density and enzyme cut site number
cor(cor_snp_sim$number,cor_snp_sim$SNP_COUNT)
#[1] 0.8165231

cor.test(cor_snp_sim$number,cor_snp_sim$SNP_COUNT)

#	Pearson's product-moment correlation

#data:  cor_snp_sim$number and cor_snp_sim$SNP_COUNT
#t = 82.616, df = 3412, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8050276 0.8274057
#sample estimates:
#      cor 
#0.8165231 

#str(cor.test(cor_snp_sim$number,cor_snp_sim$SNP_COUNT))
#List of 9
# $ statistic  : Named num 82.6
#

library(ggplot2)
png(paste0(out,"cor_sim_snp.png"),width=7,height=6,units="in",res=600,type="cairo")
p<-ggplot(data=cor_snp_sim,aes(x=number,y=SNP_COUNT))+geom_point(alpha = 0.3)
p<-p+xlab("Enzyme cut sites number")+ylab("SNP number")
p<-p+geom_smooth(method = "lm",color="red")

p<-p+annotate("text", label = "r:0.82 \n p-value < 2.2e-16", x = 100, y = 600, size = 5)
p<-p+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ylim(0,620)
p
dev.off()

# LD decay Figure 2
# The LD was analyzed using PopLDdecay (https://github.com/BGI-shenzhen/PopLDdecay)
# The code used to plot LD decay 
# The real dir was replaced with path
setwd(" ")
ld_gn<-"/path/PopLDdecay -InVCF LD_snp.vcf -MaxDist 600 -MAF 0.1 -OutStat LD01_LD_snp.vcf"
system(ld_gn)

bp_file<-c("ld_ld.stat","chr_1.stat","chr_2.stat","chr_3.stat","chr_4.stat","chr_5.stat","chr_6.stat","chr_7.stat","chr_8.stat","chr_9.stat","chr_10.stat")
div.file<-c("LD01_LD_snp.vcf.stat","LD01_1.stat","LD01_2.stat","LD01_3.stat","LD01_4.stat","LD01_5.stat","LD01_6.stat","LD01_7.stat","LD01_8.stat","LD01_9.stat","LD01_10.stat")

png("LD_decay_figure2.png",width=10,height=7,res=600,unit="in",type="cairo")
par(mai=c(1,1,0.5,0.5))
plot(0,0,xlim=c(0,600000),ylim=c(0,0.6),axes=F, type="n",xlab="Distance KB",ylab=expression(r^2),cex.lab=1.5)
axis(1,at=seq(0,600000,length=7),labels=seq(0,600,length=7))
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6),labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6))
col<-c("red","chocolate","coral4","cornflowerblue","cornsilk4","cyan2","darkblue","darkgoldenrod2","blue","blueviolet","darkolivegreen4")
for(i in c(2:11,1)){
ld<-read.table(bp_file[i],header=F)
smoothingSpline = smooth.spline(ld$V1,ld$V2, spar=0.05)
lines(smoothingSpline,col=col[i],lwd=ifelse(i==1,2,1),lty=2)
}

for(i in c(2:11,1)){
ld<-read.table(paste0("/path/",div.file[i]),header=F)
smoothingSpline = smooth.spline(ld$V1,ld$V2, spar=0.05)
lines(smoothingSpline,col=col[i],lwd=ifelse(i==1,2,1))
}

legend("topright",legend=c("Whole genome","Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),pch=15,col=col,bty="n",cex=1.5)
legend(300000,0.6,legend=c("RILS","Landraces"),lty=c(2,1),bty="n",lwd=3,cex=1.5)
dev.off()


# LD heatmmap plot using LDheatmap R package 
# Using Y1 gene as example
# Color from blue to red mean R^2 from 0 to 1

library(LDheatmap)
load("snp_dis_y1.1.rda")
load("snp_data_y1.rda")
png("heatmap_Y1.png",width=6,height=6,res=600,units="in")
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
MyHeatmap<-LDheatmap(snp_data,snp_dis,flip=TRUE,color=rgb.palette(20),title="",name = "myLDgrob",SNP.name="S2_57610965")
LD.grob1 <- editGrob(MyHeatmap$LDheatmapGrob, gPath("heatMap", "title"), gp = gpar(cex = 1.25, col = "blue"))
LD.grob2 <- editGrob(LD.grob1, gPath("Key","title"), gp = gpar(cex = 1.5, col = "black"))
LD.grob3 <- editGrob(LD.grob2, gPath("Key", "labels"),gp = gpar(cex = 1.3, col = "black"))
grid.newpage() 
grid.draw(LD.grob3)
dev.off()

