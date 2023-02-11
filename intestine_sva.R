## transcript-levels DE Analysis
## merge 
## intestine
## cecum, ileum, jejunum 

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
pheno_intestine = read.table("pheno1109_intestine.txt", header=TRUE)
head(pheno_intestine)

## load the transcript data
trans.expr = read.table("matrix_expr_transcript_bg_FPKM_1109.txt", 
                        header=T, na.strings = "NA")
head(trans.expr)
dim(trans.expr)    # 116596   41

## select only the sample that we want to analyze 
only_intestine = (colnames(trans.expr) %in% pheno_intestine$Sample_name) 
table(only_intestine)

mRNA_intestine = trans.expr[,only_intestine]
dim(mRNA_intestine)   # 116596    6
head(mRNA_intestine)

## create transcript_ID 
transcript = as.data.frame(trans.expr$transcript_ID)
colnames(transcript) = "transcript_ID"
head(transcript)

## combine transcript_ID + data 
mRNA2_intestine = cbind(transcript, mRNA_intestine)
head(mRNA2_intestine)
dim(mRNA2_intestine)  # 116596  7

## transform rownames
row.names(mRNA2_intestine) = mRNA2_intestine[,1] 
head(mRNA2_intestine)

##
mRNA2_intestine = mRNA2_intestine[,-1]          
str(mRNA2_intestine)
head(mRNA2_intestine)
View(mRNA2_intestine)


##---------------------
## filter zero: remove genes non-expressed and low expressed
## total counts per gene
Totalcounts = rowSums(mRNA2_intestine)
## genes with zero count?
table(Totalcounts == 0)
## filter genes with 0 counts
rm = rowMeans(mRNA2_intestine) == 0
mRNA_intestine.filter = mRNA2_intestine[!rm,]
dim(mRNA_intestine.filter) #64922     6
head(mRNA_intestine.filter)
View(mRNA_intestine.filter)

##----------------------------
## Using SVA
expr.fpkm = round(mRNA_intestine.filter) ## only to create the matrix
expr.fpkm[1:5,1:5]
dim(expr.fpkm)  #64922     6
##
pheno_intestine$Tissue = as.factor(pheno_intestine$Tissue)
##
dds = DESeqDataSetFromMatrix(countData = expr.fpkm,
                             colData = pheno_intestine, design = ~ Tissue)
##
dds
colData(dds)$Tissue
#models
mod  <- model.matrix( ~ Tissue, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
head(mod)
head(mod0)
##calculating the variables
n.sv <- num.sv(mRNA_intestine.filter, mod, method="leek")
n.sv #3

##
#Selecting the components to limma
## only n.sv=1 works
svseq = svaseq(as.matrix(mRNA_intestine.filter), mod, mod0, n.sv = 1)
svseq$sv

modSV = cbind(mod, svseq$sv)
mod0SV = cbind(mod0,svseq$sv)
## limma to fitted the values
## transform to log
log.fpkm = log(mRNA_intestine.filter + 1, base = 2)

## factors
sva = as.numeric(svseq$sv)
Tissue=as.factor(pheno_intestine$Tissue)
##---------------------
# limma
fit = lmFit(log.fpkm, mod0SV) ##without the tissue

## get residuals
data.residuals = residuals.MArrayLM(object = fit, y = log.fpkm)
data.residuals[1:5,1:5]

## PCA 
pca.counts2 = prcomp(t(log.fpkm), scale=F)
pdf(file = "pca_insteine_sva1.pdf", width = 15, height = 10)   
par(mar = c(8,6,8,6))   #bottom, left, top, right
par (cex=2)
par(las=2)
fviz_pca_ind(pca.counts2,
             col.ind = "red",
             repel=FALSE)
dev.off() 

## PCA cluster by the tissue 
col = rainbow(8)
tiff(filename = "PCA_residuals_intestine_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 15)

fviz_pca_ind(pca.counts2,
             pointsize = 3.5,
             geom.ind = "point",   
             col.ind = pheno_intestine$Tissue,   
             palette = col,
             legend.title = "Groups",
             mean.point = FALSE)
dev.off()   

##---------------------
## DE nalysis 
#--
mod1  <- model.matrix( ~ 0 + Tissue , colData(dds)) ## the new model. matrix in this case remove the intercept 

modSV1 = cbind(mod1, svseq$sv) ## including SVA vectors in the full model without the intercept 
head(modSV) ## with intercept
head(modSV1) ## without intercept

colnames(modSV1) = c("Ileum", "Jejunum", "Proximal_Cecum ", "sva") 
head (modSV1) ## check is the order is right 	

#--------------------
## fit the linear models
## including SVA vectors in the full model without the intercept
fit1 = lmFit(log.fpkm, modSV1)

contrast.matrix_intestine = cbind(
        "T1" = c(1,-1,0,rep(0,svseq$n.sv)), #deg for treatment 1
        "T2" = c(1,0,-1,rep(0,svseq$n.sv)), #deg for treatment 2
        "T3" = c(0,1,-1,rep(0,svseq$n.sv))  #deg for treatment 3
)

contrast.matrix_intestine  
##
names = c("Ileum", "Jejunum", "Proximal_Cecum ", "sva") ##replace the names in the same other, check the order

rownames(contrast.matrix_intestine) = names 
contrast.matrix_intestine

##----
## analysis
## second fit

fitC = contrasts.fit(fit1,contrast.matrix_intestine)
fitC = eBayes(fitC)
results = summary(decideTests(fitC, method = "separate", adjust.method = "fdr", 
                              p.value = 0.05))

write.table(results, file = "result_DE_analysis_intestine_sva1.txt", sep ="\t")

##---------------------
## MAplot

pdf(file = "plotMA_intestine_sva1.pdf", width = 20, height = 15)
par(mar = c(8,6,8,6)) #bottom, left, top, right
par (cex=2)
par(las=2)
limma::plotMA(fitC, coef = "T1")
limma::plotMA(fitC, coef = "T2")
limma::plotMA(fitC, coef = "T3")

dev.off()

### results
#T1
resuts.T1 = topTable(fitC, adjust.method = "fdr",
                     coef= "T1",
                     number = nrow(log.fpkm))
write.table(resuts.T1, file = "Results_T1_intestine_sva1.txt", sep ="\t")

#T2
resuts.T2 = topTable(fitC, adjust.method = "fdr",
                     coef= "T2",
                     number = nrow(log.fpkm))
write.table(resuts.T2, file = "Results_T2_intestine_sva1.txt", sep ="\t")

#T3
resuts.T3 = topTable(fitC, adjust.method = "fdr",
                     coef= "T3",
                     number = nrow(log.fpkm))
write.table(resuts.T3, file = "Results_T3_intestine_sva1.txt", sep ="\t")

#
##*****************************************************************************

#T1
resg1 = resuts.T1
LFC = abs(resg1$logFC)
plot.data1 = cbind(resg1,LFC)
plot.data1$thershold = plot.data1$adj.P.Val < 0.05
plot.data1$thershold = as.factor(plot.data1$thershold)

tiff(filename = "T1-intestine_sva1.tiff", res=600,
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

tiff(filename = "T2-intestine_sva1.tiff", res=600,
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

tiff(filename = "T3-intestine_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 10,
)
ggplot(data=plot.data3,
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
dim(res1) #0 
write.table(res1, file = "T1_results_sig_intestine_sva1.txt", sep ="\t")
##
#T2
common2 = ddply(plot.data2,. (thershold), nrow) 
res2 = plot.data2[plot.data2$thershold == "TRUE",]
dim(res2) #0 
write.table(res2, file = "T2_results_sig_intestine_sva1.txt", sep ="\t")
##
#T3
common3 = ddply(plot.data3,. (thershold), nrow) 
res3 = plot.data3[plot.data3$thershold == "TRUE",]
dim(res3) #0
write.table(res3, file = "T3_results_sig_intestine_sva1.txt", sep ="\t")

##******************************************************************************
##Diagrama de venn

