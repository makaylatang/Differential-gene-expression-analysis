## transcript-levels DE Analysis
## merge 
## macrophage - tissues & monophage 
## kidney_macrophage, Macrophage(lung), Macrophage(spleen)

### only remove 0
### using sva
### p.value threshold 0.05

##----------------------------
## loand the library
library(ggplot2)
library(factoextra)
library(limma)
library(DESeq2)
library(heatmaply)
library(RColorBrewer)
library(edgeR)
library(dplyr)
library(sva)


##----------------------------
## load the pheno table
pheno_macro_t = read.table("pheno1109_mm_tissue.txt", header=TRUE) 
head(pheno_macro_t)

## load the transcript data
trans.expr = read.table("matrix_expr_transcript_bg_FPKM_1109.txt", 
                        header=T, na.strings = "NA")
head(trans.expr)
dim(trans.expr)    # 116596   41

## select only the sample that we want to analyze 
only_macro_t = (colnames(trans.expr) %in% pheno_macro_t$Sample_name) 
table(only_macro_t)

mRNA_macro_t = trans.expr[,only_macro_t]
dim(mRNA_macro_t)   # 116596    8
head(mRNA_macro_t)

## add transcript_ID 
transcript = as.data.frame(trans.expr$transcript_ID)
colnames(transcript) = "transcript_ID"
head(transcript)

## combine transcript_ID + data 
mRNA2_macro_t = cbind(transcript, mRNA_macro_t)
head(mRNA2_macro_t)
dim(mRNA2_macro_t)  # 116596    9

## transform rownames
row.names(mRNA2_macro_t) = mRNA2_macro_t[,1] 
head(mRNA2_macro_t)

##
mRNA2_macro_t = mRNA2_macro_t[,-1]          
str(mRNA2_macro_t)
head(mRNA2_macro_t)
View(mRNA2_macro_t)

##---------------------
## filter zero: remove genes non-expressed and low expressed
## total counts per gene
Totalcounts = rowSums(mRNA2_macro_t)
## genes with zero count?
table(Totalcounts == 0)
## filter genes with 0 counts
rm = rowMeans(mRNA2_macro_t) == 0
mRNA_macro_t.filter = mRNA2_macro_t[!rm,]
dim(mRNA_macro_t.filter) # 76791   8
head(mRNA_macro_t.filter)

##----------------------------
## Using SVA
expr.fpkm = round(mRNA_macro_t.filter) ## only to create the matrix
head(expr.fpkm)
dim(expr.fpkm)  #76791   8
##
pheno_macro_t$Tissue = as.factor(pheno_macro_t$Tissue)
##
dds = DESeqDataSetFromMatrix(countData = expr.fpkm,
                             colData = pheno_macro_t, design = ~ Tissue)
##
dds
colData(dds)$Tissue
#models
mod  <- model.matrix( ~ Tissue, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
head(mod)
head(mod0)
##calculating the variables
n.sv <- num.sv(mRNA_macro_t.filter, mod, method="leek")
n.sv #4

##
#Selecting the components to limma
svseq = svaseq(as.matrix(mRNA_macro_t.filter), mod, mod0, n.sv = 1)
svseq$sv

modSV = cbind(mod, svseq$sv)
mod0SV = cbind(mod0,svseq$sv)

## limma to fitted the values
## transform to log
log.fpkm = log(mRNA_macro_t.filter + 1, base = 2)

## factors
sva = as.numeric(svseq$sv)
Tissue = as.factor(pheno_macro_t$Tissue)
##---------------------
# limma
fit = lmFit(log.fpkm, mod0SV) ##without the tissue

## get residuals
data.residuals = residuals.MArrayLM(object = fit, y = log.fpkm)
head(data.residuals)

##---------------------
## box plot
pdf("boxplot_macrophage_t_sva.pdf", width = 12, height = 10) 
boxplot(data.residuals, cex = 0.8, cex.axis = 0.4,
        main="macrophage_tissues and monophage - log_fpkm value") +
        theme(axis.text.x = element_text(angle = 90))
dev.off()

##---------------------
## PCA 
pca.counts2 = prcomp(t(data.residuals), scale=F)
pdf(file = "pca_macrophage_t_sva1.pdf", width = 15, height = 10)   
par(mar = c(8,6,8,6))   #bottom, left, top, right
par (cex=2)
par(las=2)
fviz_pca_ind(pca.counts2,
             col.ind = "red",
             repel=FALSE)
dev.off() 

## PCA cluster by the tissue 
col=rainbow(8)

tiff(filename = "PCA_groups_residuals_macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10)

fviz_pca_ind(pca.counts2,
             pointsize = 2,
             geom.ind = "point",   
             col.ind = pheno_macro_t$Tissue,   
             palette = col,
             legend.title = "Groups",
             mean.point = FALSE)
dev.off()   

##---------------------
## DE nalysis 
mod1  <- model.matrix( ~ 0 + Tissue , colData(dds)) ## the new model. matrix in this case remove the intercept 

modSV1 = cbind(mod1, svseq$sv) ## including SVA vectors in the full model without the intercept 
head(modSV) ## with intercept
head(modSV1) ## without intercept

colnames(modSV1) = c("macrophage_kidney", "Macrophage_lung", "Macrophage_spleen",
                     "Monocyte", "sva1") ## replace for the correct order based on your model.matrix, this only to remove the "Sample" before the sample name id 

head (modSV1) ## check is the order is right 	

#--------------------
## fit the linear models
## including SVA vectors in the full model without the intercept
fit1 = lmFit(log.fpkm, modSV1)


contrast.matrix = cbind(
        "T1" = c(1,-1,0,0,rep(0,svseq$n.sv)),# deg for treatment 1
        "T2" = c(1,0,-1,0,rep(0,svseq$n.sv)), #deg for treatment 2
        "T3" = c(1,0,0,-1,rep(0,svseq$n.sv)), #deg for treatment 3
        "T4" = c(0,1,-1,0, rep(0,svseq$n.sv)),
        "T5" = c(0,1,0,-1, rep(0,svseq$n.sv)),
        "T6" = c(0,0,1,-1, rep(0,svseq$n.sv))
)

contrast.matrix ##save this order and send to me 
##
names= c("macrophage_kidney", "Macrophage_lung", "Macrophage_spleen",
         "Monocyte", "sva1") ##replace the names in the same other  , check the order

rownames(contrast.matrix) = names 
contrast.matrix

## analysis
## second fit

fitC = contrasts.fit(fit1,contrast.matrix)
fitC = eBayes(fitC)
results = summary(decideTests(fitC, method = "separate", adjust.method = "fdr", 
                              p.value = 0.05))

write.table(results, file = "result_DE_analysis_macrophage_t_sva1.txt", sep ="\t")


##---------------------
## MAplot

pdf(file = "plotMA_macrophage_t_sva1.pdf", width = 20, height = 15)
par(mar = c(8,6,8,6)) #bottom, left, top, right
par (cex=2)
par(las=2)
limma::plotMA(fitC, coef = "T1")
limma::plotMA(fitC, coef = "T2")
limma::plotMA(fitC, coef = "T3")
limma::plotMA(fitC, coef = "T4")
limma::plotMA(fitC, coef = "T5")
limma::plotMA(fitC, coef = "T6")

dev.off()

### results
#T1
resuts.T1 = topTable(fitC, adjust.method = "fdr",
                     coef= "T1",
                     number = nrow(log.fpkm))
write.table(resuts.T1, file = "Results_T1_macrophage_t_sva1.txt", sep ="\t")

#T2
resuts.T2 = topTable(fitC, adjust.method = "fdr",
                     coef= "T2",
                     number = nrow(log.fpkm))
write.table(resuts.T2, file = "Results_T2_macrophage_t_sva1.txt", sep ="\t")

#T3
resuts.T3 = topTable(fitC, adjust.method = "fdr",
                     coef= "T3",
                     number = nrow(log.fpkm))
write.table(resuts.T3, file = "Results_T3_macrophage_t_sva1.txt", sep ="\t")

#T4
resuts.T4 = topTable(fitC, adjust.method = "fdr",
                     coef= "T4",
                     number = nrow(log.fpkm))
write.table(resuts.T4, file = "Results_T4_macrophage_t_sva1.txt", sep ="\t")

#T5
resuts.T5 = topTable(fitC, adjust.method = "fdr",
                     coef= "T5",
                     number = nrow(log.fpkm))
write.table(resuts.T5, file = "Results_T5_macrophage_t_sva1.txt", sep ="\t")

#T6
resuts.T6 = topTable(fitC, adjust.method = "fdr",
                     coef= "T6",
                     number = nrow(log.fpkm))
write.table(resuts.T6, file = "Results_T6_macrophage_t_sva1.txt", sep ="\t")

#
##*****************************************************************************

#T1
resg1 = resuts.T1
LFC = abs(resg1$logFC)
plot.data1 = cbind(resg1,LFC)
plot.data1$thershold = plot.data1$adj.P.Val < 0.05
plot.data1$thershold = as.factor(plot.data1$thershold)

tiff(filename = "T1-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data1,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

#T2
resg2 = resuts.T2
LFC = abs(resg2$logFC)
plot.data2 = cbind(resg2,LFC)
plot.data2$thershold = plot.data2$adj.P.Val < 0.05
plot.data2$thershold = as.factor(plot.data2$thershold)

tiff(filename = "T2-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data2,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

#T3
resg3 = resuts.T3
LFC = abs(resg3$logFC)
plot.data3 = cbind(resg3,LFC)
plot.data3$thershold = plot.data3$adj.P.Val < 0.05
plot.data3$thershold = as.factor(plot.data3$thershold)

tiff(filename = "T3-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data3,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

#T4
resg4 = resuts.T4
LFC = abs(resg4$logFC)
plot.data4 = cbind(resg4,LFC)
plot.data4$thershold = plot.data4$adj.P.Val < 0.05
plot.data4$thershold = as.factor(plot.data4$thershold)

tiff(filename = "T4-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data4,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

#T5
resg5 = resuts.T5
LFC = abs(resg5$logFC)
plot.data5 = cbind(resg5,LFC)
plot.data5$thershold = plot.data5$adj.P.Val < 0.05
plot.data5$thershold = as.factor(plot.data5$thershold)

tiff(filename = "T5-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data5,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

#T6
resg6 = resuts.T6
LFC = abs(resg6$logFC)
plot.data6 = cbind(resg6,LFC)
plot.data6$thershold = plot.data6$adj.P.Val < 0.05
plot.data6$thershold = as.factor(plot.data6$thershold)

tiff(filename = "T6-macrophage_t_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data6,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
dev.off()

##******************************************************************************
#T1
library(plyr)
common1 = ddply(plot.data1,. (thershold), nrow) 
res1 = plot.data1[plot.data1$thershold == "TRUE",]
dim(res1) #30237
write.table(res1, file = "T1_results_sig_macrophage_t_sva1.txt", sep ="\t")
##
#T2
common2 = ddply(plot.data2,. (thershold), nrow) 
res2 = plot.data2[plot.data2$thershold == "TRUE",]
dim(res2) #30932 
write.table(res2, file = "T2_results_sig_macrophage_t_sva1.txt", sep ="\t")
##
#T3
common3 = ddply(plot.data3,. (thershold), nrow) 
res3 = plot.data3[plot.data3$thershold == "TRUE",]
dim(res3) #36275
write.table(res3, file = "T3_results_sig_macrophage_t_sva1.txt", sep ="\t")

#T4
common4 = ddply(plot.data4,. (thershold), nrow) 
res4 = plot.data4[plot.data4$thershold == "TRUE",]
dim(res4) #1301
write.table(res4, file = "T4_results_sig_macrophage_t_sva1.txt", sep ="\t")

#T5
common5 = ddply(plot.data5,. (thershold), nrow) 
res5 = plot.data5[plot.data5$thershold == "TRUE",]
dim(res5) #5270
write.table(res5, file = "T5_results_sig_macrophage_t_sva1.txt", sep ="\t")

#T6
common6 = ddply(plot.data6,. (thershold), nrow) 
res6 = plot.data6[plot.data6$thershold == "TRUE",]
dim(res6) #7172
write.table(res6, file = "T6_results_sig_macrophage_t_sva1.txt", sep ="\t")
##******************************************************************************
