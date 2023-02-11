## transcript-levels DE Analysis
## merge 
## meat
## dark meat, white meat 

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

##----------------------------
## load the pheno table
pheno_meat = read.table("pheno1109_meat.txt", header=TRUE) 
head(pheno_meat)

## load the transcript data
trans.expr = read.table("matrix_expr_transcript_bg_FPKM_1109.txt", 
                        header=T, na.strings = "NA")
head(trans.expr)
dim(trans.expr)    # 116596   41

## select only the sample that we want to analyze 
only_meat = (colnames(trans.expr) %in% pheno_meat$Sample_name) 
table(only_meat)

mRNA_meat = trans.expr[,only_meat]
dim(mRNA_meat)   # 116596    4
head(mRNA_meat)

## add transcript_ID 
transcript = as.data.frame(trans.expr$transcript_ID)
colnames(transcript) = "transcript_ID"
head(transcript)

## combine transcript_ID + data 
mRNA2_meat = cbind(transcript, mRNA_meat)
head(mRNA2_meat)
dim(mRNA2_meat)  # 116596    5

## transform rownames
row.names(mRNA2_meat) = mRNA2_meat[,1] 
head(mRNA2_meat)

##
mRNA2_meat = mRNA2_meat[,-1]          
str(mRNA2_meat)
head(mRNA2_meat)
View(mRNA2_meat)

##---------------------
## filter zero: remove genes non-expressed and low expressed
## total counts per gene
Totalcounts = rowSums(mRNA2_meat)
## genes with zero count?
table(Totalcounts == 0)
## filter genes with 0 counts
rm = rowMeans(mRNA2_meat) == 0
mRNA_meat.filter = mRNA2_meat[!rm,]
dim(mRNA_meat.filter) #57045     4
head(mRNA_meat.filter)

##----------------------------
## Using SVA
expr.fpkm = round(mRNA_meat.filter) ## only to create the matrix
head(expr.fpkm)
dim(expr.fpkm)  #57045   4
##
pheno_meat$Tissue = as.factor(pheno_meat$Tissue)
##
dds = DESeqDataSetFromMatrix(countData = expr.fpkm,
                             colData = pheno_meat, design = ~ Tissue)
##
dds
colData(dds)$Tissue
#models
mod  <- model.matrix( ~ Tissue, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
head(mod)
head(mod0)
##calculating the variables
n.sv <- num.sv(mRNA_meat.filter, mod, method="leek")
n.sv #2

##
#Selecting the components to limma
svseq = svaseq(as.matrix(mRNA_meat.filter), mod, mod0, n.sv = 1)
svseq$sv

modSV = cbind(mod, svseq$sv)
mod0SV = cbind(mod0,svseq$sv)

## limma to fitted the values
## transform to log
log.fpkm_meat = log(mRNA_meat.filter + 1, base = 2)

## factors
sva = as.numeric(svseq$sv)
Tissue = as.factor(pheno_meat$Tissue)
##---------------------
# limma
fit = lmFit(log.fpkm_meat, mod0SV) ##without the tissue

## get residuals
data.residuals = residuals.MArrayLM(object = fit, y = log.fpkm_meat)
head(data.residuals)

##---------------------
## box plot
pdf("boxplot_meat_sva.pdf", width = 12, height = 10) 
boxplot(data.residuals, cex = 0.8, cex.axis = 0.8,
        main="meat - log_fpkm value") +
        theme(axis.text.x = element_text(angle = 90))
dev.off()

##---------------------
## PCA 
pca.counts2 = prcomp(t(data.residuals), scale=F)
pdf(file = "pca_meat_sva.pdf", width = 15, height = 10)   
par(mar = c(8,6,8,6))   #bottom, left, top, right
par (cex=2)
par(las=2)
fviz_pca_ind(pca.counts2,
             col.ind = "red",
             repel=FALSE)
dev.off() 

## PCA cluster by the tissue 
col=rainbow(8)

tiff(filename = "PCA_groups_residuals_meat_sva1.tiff", res=600,
     width = 4500, height = 3500, pointsize = 15)

fviz_pca_ind(pca.counts2,
             pointsize = 3.5,
             geom.ind = "point",   
             col.ind = pheno_meat$Tissue,   
             palette = col,
             legend.title = "Groups",
             mean.point = FALSE)
dev.off()   

##---------------------
## DE nalysis 
# mod  <- model.matrix( ~ Tissue, colData(dds))
# modSV = cbind(mod, svseq$sv)


## fit the linear models
fit1 = lmFit(log.fpkm_meat, modSV)

contrast.matrix = cbind("T" = c(1,-1, rep(0,svseq$n.sv)))
contrast.matrix

names = c("Dark_meat", "White_meat", "sva")
rownames(contrast.matrix) = names 
contrast.matrix

## analysis
## second fit

fitC = contrasts.fit(fit1,contrast.matrix)
fitC = eBayes(fitC)
results = summary(decideTests(fitC, method = "separate", adjust.method = "fdr", 
                              p.value = 0.05))

write.table(results, file = "result_DE_analysis_meat_sva.txt", sep ="\t")

##---------------------
## MAplot

pdf(file = "plotMA_meat_sva.pdf", width = 20, height = 15)
par(mar = c(8,6,8,6)) #bottom, left, top, right
par (cex=2)
par(las=2)
limma::plotMA(fitC, coef = "T")
dev.off()

### results
#T
resuts.T = topTable(fitC, adjust.method = "fdr",
                    coef= "T",
                    number = nrow(log.fpkm_meat))
write.table(resuts.T, file = "Results_T_meat_sva.txt", sep ="\t")

##*****************************************************************************

## plots
#T
resg = resuts.T
LFC = abs(resg$logFC)
plot.data = cbind(resg,LFC)
plot.data$thershold = plot.data$adj.P.Val < 0.05
plot.data$thershold = as.factor(plot.data$thershold)


tiff(filename = "T_meat_sva3.tiff", res=600,
     width = 4500, height = 3500, pointsize = 8)
ggplot(data=plot.data,
       aes(x=logFC, y=-log10(adj.P.Val))) +
        geom_point(aes(colour = thershold), size= 1.5) +
        scale_color_brewer(palette = "Set1",
                           name = "FDR < 0.05 ",
                           labels = c("Non-significant","Significant" ))
       
dev.off()

##******************************************************************************
##*
# T1
## adj.P.Val < 0.05
library(plyr)
common = ddply(plot.data,. (thershold), nrow) 
res = plot.data[plot.data$thershold == "TRUE",]
dim(res) # 3000 8
write.table(res, file = "T_results_sig_meat_sva.txt", sep ="\t")

##******************************************************************************
## prepare for running STRING
## i don't need for now
## skep below code! 
library(dplyr)
DE_T1 <- read.csv("Results_T1_2_p0.1.csv")
View(DE_T1)
DE_T1_sig <- DE_T1 %>% filter(adj.P.Val < 0.01)
DE_T1_sig_p0.1_F2 <- DE_T1 %>% filter(adj.P.Val < 0.1 & abs(logFC) > 2)

write.table(DE_T1_sig, file = "DE_T1_sig.txt")
write.table(DE_T1_sig_p0.1_F2, file = "DE_T1_sig_p0.1_F2.txt")

## T2
DE_T2 <- read.csv("Results_T2_2_p0.1.csv")
View(DE_T2)
DE_T2_sig_p0.1_F2 <- DE_T2 %>% filter(adj.P.Val < 0.1 & abs(logFC) > 2)
write.table(DE_T2_sig_p0.1_F2, file = "DE_T2_sig_p0.1_F2.txt")

r1 <- read.table("~/Desktop/Hawkins_merge/meat_nosva/T_meat_results_sig.txt")
r2 <- read.table("~/Desktop/Hawkins_merge/meat_nosva/Results_meat.txt")

nrow(r1)
nrow(r2)
View(r2)
View(r1)







nrow(r2 %>% filter(adj.P.Val < 0.05))

r3 <- read.table("~/Desktop/Hawkins_merge/T_results_sig_meat_sva.txt")
nrow(r3)
