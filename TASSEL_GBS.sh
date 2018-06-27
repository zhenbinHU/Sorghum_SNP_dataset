#!/bin/bash
#SBATCH --job-name=sra2fq
#SBATCH --output=sra2fq.txt
#SBATCH --time=10-10:10:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=killable.q
#SBATCH --mem-per-cpu=20GB
#SBATCH --mail-user=zhenbin@ksu.edu
module load Java
# the real dir was replaced using path
tassel="/path/to/run_pipeline.pl"
reference="/path/to/final.Sbicolor_313_v3.1.fa"
key="/path/to/keyfile"
bwa="/homes/zhenbin/bwa-0.7.10/bwa"
mkdir /path/GBS_calling
mkdir /path/GBS_calling/alignment
mkdir /path/GBS_calling/database
mkdir /path/GBS_calling/hapmap
mkdir /path/GBS_calling/h5
cd /path/GBS_calling
$tassel -Xms150g -Xmx200G -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i fastq -db database/GBSclear1.db -k clear_key.txt -kmerLength 64 -minKmerL 20 -mnQS 20 -mxKmerNum 100000000 -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -TagExportToFastqPlugin -db database/GBSclear1.db  -o alignment/tagsForAlign_clear1.fa.gz -endPlugin -runfork1 
$bwa aln -t 64 $reference alignment/tagsForAlign_clear1.fa.gz > alignment/tagsForAlign_clear1.sai
$bwa samse $reference alignment/tagsForAlign_clear1.sai alignment/tagsForAlign_clear1.fa.gz > alignment/tagsForAlign_clear1.sam
$tassel -Xms150g -Xmx200G -fork1 -SAMToGBSdbPlugin -i alignment/tagsForAlign_clear1.sam -db database/GBSclear1.db -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -DiscoverySNPCallerPluginV2 -db database/GBSclear1.db -sC 1 -eC 10 -mnLCov 0.1 -mnMAF 0.000001 -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -SNPQualityProfilerPlugin -db database/GBSclear1.db -statFile "sorghum_gbsclear1_Stats.txt" -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -UpdateSNPPositionQualityPlugin -db database/GBSclear1.db -qsFile h5/myclear1QsFile -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -GetTagSequenceFromDBPlugin -db database/GBSclear1.db -o h5/allclear1TagFile.txt -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -GetTagTaxaDistFromDBPlugin -db database/GBSclear1.db -o h5/clear1TagTaxaDistOutput.txt -endPlugin -runfork1
$tassel -Xms150g -Xmx200G -fork1 -ProductionSNPCallerPluginV2 -i fastq -k clear_key.txt -e ApeKI -db database/GBSclear1.db -o h5/GBSclear1.h5 -endPlugin -runfork1 
$tassel -Xms150g -Xmx200G -fork1 -h5 h5/GBSclear1.h5 -export hapmap/GBSclear1 -exportType Hapmap -runfork1
