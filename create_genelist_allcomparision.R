#! /usr/bin/env Rscript
#When: 2018-07-09
#Who: Shalu Jhanwar
#What: Append geneNames to geneId of Differentially expressed genes

#################
#Inputs and set paths
#################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
inputDir=args[1] #Dir containing DEG files
gtf=args[2]

##################
#Load libraries
##################
library(plyr)
library(refGenome)
library(data.table)
options(warn=1)

#Read the gtf file of the respective genomes
ens <- ensemblGenome()
read.gtf(ens, gtf, useBasedir=FALSE)
gtfDf <- getGenePositions(ens)
newDf=data.frame(gene_name=gtfDf$gene_name, geneId=gtfDf$gene_id, chrName=gtfDf$seqid, startPos=gtfDf$start, endPos=gtfDf$end, strandInfo=gtfDf$strand)

#Export all genes bedFile
setwd(inputDir)
bedGene=newDf[,c(3:5,2,1,6)]
write.table(bedGene, file="AllEnsembleGenes.bed", quote=F, col.names=F, row.names=F, sep='\t')

fc=list()
df=list()

#Create a master df for all the DEG in 4 comparisons
sampleList=list.files(inputDir, recursive = TRUE, pattern="_FC.txt")
sampleList=grep("GLM", sampleList, value = TRUE)
for (i in 1:length(sampleList)){
	colString=gsub('(.*)_.*','\\1', basename(sampleList[i]))
	df[[colString]]=read.table(sampleList[i], header=T, sep='\t')
	setnames(df[[colString]], "regNames", "geneId")
}

#Map ensemble Id to gene list for all the comparisons
for (i in 1:length(df)){
       dd = merge(x=newDf, y=df[[i]], by="geneId")
	dd = dd[,c(7:ncol(dd),3:5,1,2,6)]
	ddBed = dd[,c((ncol(dd)-5):ncol(dd),1:(ncol(dd)-6))]
	write.table(dd, file=paste0(names(df)[i],"_withCoordGeneName.txt"), quote=F, col.names=T, row.names=F, sep='\t')
	write.table(ddBed, file=paste0(names(df)[i],"_withCoordGeneName.bed"), quote=F, col.names=T, row.names=F, sep='\t')
}

