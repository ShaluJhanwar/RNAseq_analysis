#! /usr/bin/env Rscript
#When: 2018-07-09
#Who: Shalu Jhanwar
#What: Determine differentially expressed genes using fisher excat, glm methods
        #Remove scRNA genes
        #cpm>=1 in all 3 rep of any of the samples

#################
#Inputs and set paths
#################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
inputDir=args[1] #/scicore/home/zellerro/jhanwar/ERC/RNA/Chicken only three stages
scRNAgeneList=args[2] #/scicore/home/zellerro/GROUP/Resources/galGal5/RNASeq/smallRNA_geneId_ensembleToUCSC_Galgal_5_0_91.gtf
scRNA=read.csv(scRNAgeneList, header=F)
setwd(inputDir)

##################
#Load libraries
##################
library("edgeR")
library("gplots")
library("methods")
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(VennDiagram)
library(ggplot2)
options(warn=1)

#####################
#Output Files and Dir
####################
dir.create(file.path(inputDir, "differentialExpression"), showWarnings = TRUE)
fisherExactPath=paste0(inputDir,"/differentialExpression")

#####################
#Define Functions
#####################
pvalFDRwriteFiles_exact<- function(filePrefix005, lrtObj, tmmObj, DEGobj){
        pVal_lrtObj=lrtObj$table[,"PValue"]
        bhs_lrtObj=p.adjust(pVal_lrtObj, method = "BH")
        out_bhs_lrtObj=cbind(DEGobj$counts, lrtObj$table, bhs_lrtObj)
        colnames(out_bhs_lrtObj)[ncol(out_bhs_lrtObj)] <- "FDR"

        fdr05 = out_bhs_lrtObj[abs(out_bhs_lrtObj$logFC) >= fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]
        up_fdr05 = out_bhs_lrtObj[out_bhs_lrtObj$logFC >= fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]
        down_fdr05 = out_bhs_lrtObj[out_bhs_lrtObj$logFC <= -fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]

        a = fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_fdr05 = dfExport_fdr05[,c(6:11,2:5,1)]
        write.table(dfExport_fdr05, file=paste0(filePrefix005,"_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_fdr05'))

        a = up_fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_up_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_up_fdr05 = dfExport_up_fdr05[,c(6:11,2:5,1)]
        write.table(dfExport_up_fdr05, file=paste0(filePrefix005,"_UP_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_up_fdr05'))

        a = down_fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_down_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_down_fdr05 = dfExport_down_fdr05[,c(6:11,2:5,1)]
        write.table(dfExport_down_fdr05, file=paste0(filePrefix005,"_DOWN_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_down_fdr05'))


	#Plot MA
        out = paste0(filePrefix005, "_EXCAT_chicken_MA_plot.pdf")
        pdf(out, width=6, height=6)
        par(mar=c(6,4,6,4))
        FC_cutoff=out_bhs_lrtObj[,"logFC"]
        plot(out_bhs_lrtObj[, "logCPM"], out_bhs_lrtObj[, "logFC"], pch = "*", col = ifelse((abs(out_bhs_lrtObj$logFC) >= fc_t & out_bhs_lrtObj[,"FDR"] <= fdr_t_5), 2, 1), xlab = "logCPM", ylab = "logFC")
        legend("topright", sprintf("p-value <%0.2f", fdr_t_5), text.col = 2, bty = "n")
        dev.off()

        #Volcano plot
        out = paste0(filePrefix005, "_EXCAT_chicken_volcano_plot.pdf")
        pdf(out, width=7, height=7)
        par(mar=c(6,6,6,6))
        plot(out_bhs_lrtObj$logFC, -log(out_bhs_lrtObj$PValue,10), pch = "*", cex = 1, xlab="log2 (fold-change)", ylab="-log10 (p-value)", xlim=c(-6,6), ylim=c(0,100), cex.axis=1.5, cex.lab=2, font.lab=2)
        points(down_fdr05$logFC, -log(down_fdr05$PValue,10), pch = "*", col="blue")
        points(up_fdr05$logFC, -log(up_fdr05$PValue,10), pch = "*", col="red")
        legend(0.6,90, "Down-regulated", text.font=2, text.col = "blue", bty = "n", cex = 1.2)
        legend(0.6,98, "Up-regulated", text.font=2, text.col = "red", bty = "n", cex = 1.2)
        dev.off()

        return(fdr05)
}

pvalFDRwriteFiles <- function(filePrefix005, lrtObj, tmmObj, DEGobj){
	#"GLM_H1920vsH22_fdr05", lrt_H1920vsH22, tmmCounts, genes
        pVal_lrtObj=lrtObj$table[,"PValue"]
        bhs_lrtObj=p.adjust(pVal_lrtObj, method = "BH")
        out_bhs_lrtObj=cbind(DEGobj$counts, lrtObj$table, bhs_lrtObj)
        colnames(out_bhs_lrtObj)[ncol(out_bhs_lrtObj)] <- "FDR"

        fdr05 = out_bhs_lrtObj[abs(out_bhs_lrtObj$logFC) >= fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]
        up_fdr05 = out_bhs_lrtObj[out_bhs_lrtObj$logFC >= fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]
        down_fdr05 = out_bhs_lrtObj[out_bhs_lrtObj$logFC <= -fc_t & out_bhs_lrtObj$FDR <= fdr_t_5,]

        a = fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_fdr05 = dfExport_fdr05[,c(6:14,2:5,1)]
        write.table(dfExport_fdr05, file=paste0(filePrefix005,"_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_fdr05'))

        a = up_fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_up_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_up_fdr05 = dfExport_up_fdr05[,c(6:14,2:5,1)]
        write.table(dfExport_up_fdr05, file=paste0(filePrefix005,"_UP_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_up_fdr05'))

        a = down_fdr05[,c("logFC","logCPM","PValue","FDR")]
        a = cbind(a, regNames = rownames(a))
        dfExport_down_fdr05 = merge(x=a, y=tmmObj, by.x="regNames", by.y="regNames", all.x=TRUE)
        dfExport_down_fdr05 = dfExport_down_fdr05[,c(6:14,2:5,1)]
        write.table(dfExport_down_fdr05, file=paste0(filePrefix005,"_DOWN_FC.txt"), sep="\t", quote=F, col.names=T, row.names=F)
        rm(list = c('a', 'dfExport_down_fdr05'))

	                #Plot MA
        out = paste0(filePrefix005, "_GLM_chicken_MA_plot.pdf")
        pdf(out, width=6, height=6)
        par(mar=c(6,4,6,4))
        FC_cutoff=out_bhs_lrtObj[,"logFC"]
        plot(out_bhs_lrtObj[, "logCPM"], out_bhs_lrtObj[, "logFC"], pch = "*", col = ifelse((abs(out_bhs_lrtObj$logFC) >= fc_t & out_bhs_lrtObj[,"FDR"] <= fdr_t_5), 2, 1), xlab = "logCPM", ylab = "logFC")
        legend("topright", sprintf("p-value <%0.2f", fdr_t_5), text.col = 2, bty = "n")
        dev.off()

        #Volcano plot
        out = paste0(filePrefix005, "_GLM_chicken_volcano_plot.pdf")
        pdf(out, width=7, height=7)
        par(mar=c(6,6,6,6))
        plot(out_bhs_lrtObj$logFC, -log(out_bhs_lrtObj$PValue,10), pch = "*", cex = 1, xlab="log2 (fold-change)", ylab="-log10 (p-value)", xlim=c(-6,6), ylim=c(0,100), cex.axis=1.5, cex.lab=2, font.lab=2)
        points(down_fdr05$logFC, -log(down_fdr05$PValue,10), pch = "*", col="blue")
        points(up_fdr05$logFC, -log(up_fdr05$PValue,10), pch = "*", col="red")
        legend(0.6,90, "Down-regulated", text.font=2, text.col = "blue", bty = "n", cex = 1.2)
        legend(0.6,98, "Up-regulated", text.font=2, text.col = "red", bty = "n", cex = 1.2)
        dev.off()

        return(fdr05)

}

####################
#set cutoffs
####################
#fc_t <- log2(2)
fc_t <- log2(1.5)
fdr_t_5 <- 0.05
noSampleCutoff=2
cpmCutoff=1

####################
#Prepare DF of count matrix
####################
colString=c()
masterDf=c()
sampleList=list.files(inputDir, recursive = TRUE, pattern=".rsamQuant.genes.results")
for (i in 1:length(sampleList)){
        colString=c(colString, gsub('(.*)_.*','\\1', basename(sampleList[i])))
        df=read.table(sampleList[i], header=T, sep='\t')
        masterDf=cbind(masterDf, headName=df$expected_count)
}

#substitute name of HH1920 as HH20
colString = gsub("19", "", colString)
masterDf=round(masterDf)
rownames(masterDf)=df$gene_id
colnames(masterDf)=colString

####################
#Design DGElist and contrast
####################
#For GLM
group=gsub('.{1}$', '', colString)
stages<-c(rep("Stages", ncol(masterDf)))
targets <- as.data.frame(cbind(group, stages))
Group<-factor(paste(targets$group,targets$stages,sep="_"))
samples=paste0("S",seq(1:length(colnames(masterDf))))
rownames(targets)<-samples
targets<-cbind(targets, Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
#Creat contrast matrix
my.contrasts <- makeContrasts(H1920vsH22=HH22_Stages-HH20_Stages, H22vsH24=HH24_Stages-HH22_Stages, H1920vsH24=HH24_Stages-HH20_Stages, levels=design)
genes <- DGEList(masterDf, group = group)
setwd(fisherExactPath)

####################
#Export cpm and TMM counts for ALL the genes (including scRNA)
####################
cpmCounts <- 1e+06 * (genes$counts/expandAsMatrix(genes$samples$lib.size, dim(genes)))
out=paste0(fisherExactPath,"/cpmCountsAllSamples.txt")
write.table(cpmCounts, file=out, sep="\t", quote=F, col.names=T, row.names=T)
allGenes <- calcNormFactors(genes, method="TMM")

tmmFile = paste0(fisherExactPath,"/tmmCounts_allGenes.txt")
tmmCounts = as.data.frame(cpm(allGenes, normalized.lib.sizes=TRUE))
tmmCounts = cbind(tmmCounts, regNames = rownames(tmmCounts))
write.table(tmmCounts, file=tmmFile, sep="\t", quote=F, col.names=T, row.names=F)

###################
#Get the clustering of the genes
genesLog <- log2(genes$counts+1)
genesLogSD=apply(genesLog, 1, sd)
percentile99 = quantile(genesLogSD, 0.99)
#Filter out the outliers
genesLog=genesLog[genesLogSD <= percentile99,]
genesLogSD=genesLogSD[genesLogSD <= percentile99]
cat(sprintf("Number of genes after 1 CPM and 99 percentile filtering : %d", length(genesLogSD)))
#Get the most variables genes
genesLogTop = genesLog[order(genesLogSD, decreasing=T),][1:1000,]
dis <- as.dist(1 - cor(genesLogTop, method="spearman"))
hc <- hclust(dis)
clusF = paste(paste0(fisherExactPath, "/chicken_hc_Clustering_plot"), "pdf", sep='.')
pdf(clusF, width=6, height=6)
par(mar=c(4,7,4,4))
plot(hc, cex.axis=2, cex.lab=2, font.lab=2, cex=1.6, xlab="Distance", ylab="Height")
dev.off()

#PCA
pca.norm_masterDf_LogTop = prcomp(t(genesLogTop), center = TRUE, scale. = TRUE)
pcaSamples=paste(paste0(fisherExactPath, "/chicken_pcaSamples"), "pdf", sep='.')
pdf(pcaSamples, width=10, height=8)
par(mar=c(7,8,7,6))
ggplot(data.frame(pca.norm_masterDf_LogTop$x, samples = factor(c("HH20","HH20","HH20","HH22","HH22","HH22","HH24","HH24","HH24"))), aes(x = PC1, y = PC2)) + geom_point(aes(shape = samples, color=samples), size = 7) + theme_bw() + theme(legend.text = element_text(size=18), axis.text=element_text(size=18, face="bold"), axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=18, face="bold"))
dev.off()
eigs=pca.norm_masterDf_LogTop$sdev^2
varPCA=eigs/sum(eigs)
pcaVar=paste(paste0(fisherExactPath, "/chicken_pca_variancePlot"), "pdf", sep='.')
pdf(pcaVar, width=8, height=8)
par(mar=c(7,7,7,5))
plot(varPCA, type='b', lwd=2.5,xlab="PC components", ylab="Proportion of Variance (PCA)", cex = 1.5, font.axis=1.5, font.lab=1.8, cex.axis=2, cex.lab=2.1)
dev.off()

#######################
#Filtering data
#scRNA 
#No. of samples: to get the consistency between the samples, get the genes with cpm>=1 in all replicates of any stage
#######################
stage1 = rowSums(cpmCounts[,c(1:3)] >= cpmCutoff)
stage2 = rowSums(cpmCounts[,c(4:6)] >= cpmCutoff)
stage3 = rowSums(cpmCounts[,c(7:9)] >= cpmCutoff)
genes = genes[(stage1 >= noSampleCutoff | stage2 >= noSampleCutoff | stage3 >= noSampleCutoff),]
cat(sprintf("Number of genes after >=1 CPM cutoff on n samples: %d", dim(genes)[1]))
scRNA=read.csv(scRNAgeneList, header=F)
genes = genes[!rownames(genes) %in% scRNA$V1, ]
cat(sprintf("Number of genes after scRNA removal: %d", dim(genes)[1]))

# Plot PCA on the filtered genes (TMM normalized values)
mdsF = paste0(fisherExactPath, "/chicken_mds_plot.pdf")
pdf(mdsF, width=6, height=6)
par(mar=c(4,7,4,4))
plotMDS(genes, gene.selection = "common", col = as.numeric(group), cex.axis=2, cex.lab=2, font.lab=2, cex=1)
dev.off()

######################
#Apply GLM framework
######################
genes <- calcNormFactors(genes, method="TMM")
tmmCounts = as.data.frame(cpm(genes, normalized.lib.sizes=TRUE))
tmmCounts = cbind(tmmCounts, regNames = rownames(tmmCounts))
genes <- estimateDisp(genes,design)

#Prepare DF for excat test (so that dispersion is calculated from the entire data only once)
genes_1920_22 <-genes[, c(1:6), keep.lib.sizes=TRUE]
genes_22_24 <-genes[, c(4:9), keep.lib.sizes=TRUE]
genes_1920_24 <-genes[, c(1:3,7:9), keep.lib.sizes=TRUE]

#GLM framework
fit <- glmFit(genes,design)
#Dispersion plot
dispF <- paste0(fisherExactPath, "/chicken_dispersion_plot.pdf")
pdf(dispF)
plotBCV(genes)
dev.off()

#If ‘NULL’ will be extracted from ‘y’, with order of precedence: genewise dispersion, trended dispersions, common dispersion.
lrt_H1920vsH22 = glmLRT(fit, contrast=my.contrasts[,"H1920vsH22"])
lrt_H22vsH24 = glmLRT(fit, contrast=my.contrasts[,"H22vsH24"])
lrt_H1920vsH24 = glmLRT(fit, contrast=my.contrasts[,"H1920vsH24"])
GLM_DEG_H1920_H22 = pvalFDRwriteFiles("GLM_H1920vsH22_fdr05", lrt_H1920vsH22, tmmCounts, genes)
GLM_DEG_H22_H24 = pvalFDRwriteFiles("GLM_H22vsH24_fdr05", lrt_H22vsH24, tmmCounts, genes)
GLM_DEG_H1920_H24 = pvalFDRwriteFiles("GLM_H1920vsH24_fdr05", lrt_H1920vsH24, tmmCounts, genes)

#Get all DEG list 
allDEGnames = unique(c(rownames(GLM_DEG_H1920_H22), rownames(GLM_DEG_H22_H24), rownames(GLM_DEG_H1920_H24)))

##################
#Generate plots: 
##################
tmmDEG = tmmCounts[rownames(tmmCounts)%in% allDEGnames, -c(ncol(tmmCounts))]
scaled_tmmDEG = t(apply(tmmDEG, 1, scale))
colnames(scaled_tmmDEG) <- colnames(tmmCounts)[-length(colnames(tmmCounts))]
dis <- as.dist(1 - cor(t(scaled_tmmDEG), method="spearman"))
hc <- hclust(dis)

#method 1
##heatmap.2(scaled_tmmDEG, margin=c(20,5), Colv=NULL, Rowv=as.dendrogram(hc), dendrogram=c("row"), col=colorpanel(20,"blue","white","red"), trace=c("none"), cexRow=0.00001, colsep=c(0, 1:dim(scaled_tmmDEG)[2]), rowsep=c(0,dim(scaled_tmmDEG)[1]), sepcolor="black")

#method2
type = group
pdf("GLM_chicken_heatmap.pdf", width=10, height=12)
par(mar=c(10,10,8,8))
ha = HeatmapAnnotation(df = data.frame(type = type), col = list(type = c("HH20" =  "brown4", "HH22" = "magenta3", "HH24" = "chartreuse4")),
        annotation_legend_param = list(type = list(title = "Stages", title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14))))
Heatmap(scaled_tmmDEG, name = "Normalized Expression", km = 4, col = colorRamp2(c(-2, 0, 2), c("blue1", "azure2", "red1")),
   top_annotation = ha, top_annotation_height = unit(8, "mm"),
   show_row_names = FALSE, show_column_names = FALSE,
   heatmap_legend_param = list(title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 14)))
dev.off()

H1920_H22 = rownames(GLM_DEG_H1920_H22)
H22_H24 = rownames(GLM_DEG_H22_H24)
H1920_H24 = rownames(GLM_DEG_H1920_H24)
pdf("GLM_chicken_vennDiagram.pdf", width=8, height=8)
par(mar=c(8,8,8,8))
v1 <-venn.diagram(list(H20_H22=H1920_H22, H22_H24=H22_H24, H20_H24=H1920_H24), filename=NULL, fill=c("blue1", "red1", "green1"), alpha = c(0.20, 0.20, 0.20), cex = 2, cat.fontface = 4,lty =0, cat.cex=1.7)
grid.newpage()
grid.draw(v1)
dev.off()

#other color palattes
#Heatmap(scaled_tmmDEG, name = "expression", km = 4, col = colorRamp2(c(-2, 0, 2), c("#E6E668", "#5D9E95", "#4F1D96")),
#Heatmap(scaled_tmmDEG, name = "expression", km = 4, col = colorRamp2(c(-2, 0, 2), c("#E6E601", "#5D9E95", "#582291")),

#####################
#Apply exact test framework
#####################
tmmCounts_1920_22 <- as.data.frame(cpm(genes_1920_22, normalized.lib.sizes=TRUE))
tmmCounts_1920_22 <- cbind(tmmCounts_1920_22, regNames = rownames(tmmCounts_1920_22))
fet_H1920vsH22 <- exactTest(genes_1920_22, pair = c("HH20", "HH22"))
FET_DEG_H1920vsH22 = pvalFDRwriteFiles_exact("FET_H1920vsH22_fdr05", fet_H1920vsH22, tmmCounts_1920_22, genes_1920_22)

tmmCounts_22_24 <- as.data.frame(cpm(genes_22_24, normalized.lib.sizes=TRUE))
tmmCounts_22_24 <- cbind(tmmCounts_22_24, regNames = rownames(tmmCounts_22_24))
fet_H22vsH24 <- exactTest(genes_22_24, pair = c("HH22", "HH24"))
FET_DEG_H22vsH24 = pvalFDRwriteFiles_exact("FET_H22vsH24_fdr05", fet_H22vsH24, tmmCounts_22_24, genes_22_24)

tmmCounts_1920_24 <- as.data.frame(cpm(genes_1920_24, normalized.lib.sizes=TRUE))
tmmCounts_1920_24 <- cbind(tmmCounts_1920_24, regNames = rownames(tmmCounts_1920_24))
fet_H1920vsH24 <- exactTest(genes_1920_24, pair = c("HH20", "HH24"))
FET_DEG_H1920vsH24 = pvalFDRwriteFiles_exact("FET_H1920vsH24_fdr05", fet_H1920vsH24, tmmCounts_1920_24, genes_1920_24)

#Get all DEG list
FET_allDEGnames = unique(c(rownames(FET_DEG_H1920vsH22), rownames(FET_DEG_H22vsH24), rownames(FET_DEG_H1920vsH24)))

##################
#Generate plots: 
##################
tmmDEG = tmmCounts[rownames(tmmCounts)%in% FET_allDEGnames, -c(ncol(tmmCounts))]
scaled_tmmDEG = t(apply(tmmDEG, 1, scale))
colnames(scaled_tmmDEG) <- colnames(tmmCounts)[-length(colnames(tmmCounts))]
dis <- as.dist(1 - cor(t(scaled_tmmDEG), method="spearman"))
hc <- hclust(dis)
#method 1
##heatmap.2(scaled_tmmDEG, margin=c(20,5), Colv=NULL, Rowv=as.dendrogram(hc), dendrogram=c("row"), col=colorpanel(20,"blue","white","red"), trace=c("none"), cexRow=0.00001, colsep=c(0, 1:dim(scaled_tmmDEG)[2]), rowsep=c(0,dim(scaled_tmmDEG)[1]), sepcolor="black")

#method2
type = group
pdf("FET_chicken_heatmap.pdf", width=10, height=12)
par(mar=c(10,10,8,8))
ha = HeatmapAnnotation(df = data.frame(type = type), col = list(type = c("HH20" =  "brown4", "HH22" = "magenta3", "HH24" = "chartreuse4")),
        annotation_legend_param = list(type = list(title = "Stages", title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14))))
Heatmap(scaled_tmmDEG, name = "Normalized Expression", km = 4, col = colorRamp2(c(-2, 0, 2), c("blue1", "azure2", "red1")),
   top_annotation = ha, top_annotation_height = unit(8, "mm"),
   show_row_names = FALSE, show_column_names = FALSE,
   heatmap_legend_param = list(title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 14)))
dev.off()

#Venn diagram
H1920_H22 = rownames(FET_DEG_H1920vsH22)
H22_H24 = rownames(FET_DEG_H22vsH24)
H1920_H24 = rownames(FET_DEG_H1920vsH24)
pdf("FET_chicken_vennDiagram.pdf", width=8, height=8)
par(mar=c(8,8,8,8))
v1 <-venn.diagram(list(H20_H22=H1920_H22, H22_H24=H22_H24, H20_H24=H1920_H24), filename=NULL, fill=c("blue1", "red1", "green1"), alpha = c(0.20, 0.20, 0.20), cex = 2, cat.fontface = 4,lty =0, cat.cex=1.7)
grid.newpage()
grid.draw(v1)
dev.off()






