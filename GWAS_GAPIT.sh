# GWAS code with GAPIT using FT as example
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler)
library("scatterplot3d") 
source("http://zzlab.net/GAPIT/emma.txt")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

FT<-read.table("FT.blup.txt",header=T)
FT$Taxa<-as.character(FT$Taxa)
FT$Taxa[FT$Taxa=="Macia"]<-"macia"
FT$Taxa[FT$Taxa=="SC35-14E"]<-"SC35.14E"
FT$Taxa[FT$Taxa=="Segaolane"]<-"segaolane"
FT$Taxa[FT$Taxa=="Ajabsido"]<-"ajabsido"

hmp<-read.delim("SB_nam.hmp.txt",header=F)
id<-as.character(unlist(hmp[1,]))
FT<-FT[which(FT$Taxa%in%id),]
myGAPIT <- GAPIT(Y=FT,G=hmp,PCA.total=5,SNP.MAF=0.05) 

!#/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=gwas.txt
#SBATCH --time=10-10:10:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=killable.q
#SBATCH --mem-per-cpu=20GB
#SBATCH --mail-user=zhenbin@ksu.edu
module load R
R CMD BATCH gwas.sh
