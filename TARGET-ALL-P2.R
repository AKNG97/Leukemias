library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)

qry.rna <- GDCquery(project = "TARGET-ALL-P2",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",)

GDCdownload(qry.rna)

dat <- qry.rna[[1]][[1]]
table(as.factor(dat$sample_type))

rnas <- GDCprepare(qry.rna, summarizedExperiment = TRUE)
rnas.raw.TAP2 <- saveRDS(rnas, "rnas_raw_TAP2.RDS")

data <- assay(rnas)
rownames(data) <- rowData(rnas)$external_gene_name
rownames(rnas) <- rowData(rnas)$external_gene_name
head(rownames(data))

dim(data)
dataFilt <- TCGAanalyze_Filtering(tabDF = data,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(data)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))    
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ]
rnas <- rnas[!duplicated(rownames(rnas)),]

annot<-read.delim(file="mart_export.txt", sep="\t")
names(annot)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot$Length <- abs(annot$End - annot$Start)
inter <- intersect(rownames(rnas), annot$Gene.name)
rnas1 <- rnas[rownames(rnas) %in% inter,]
annot1 <- annot[annot$Gene.name %in% inter,]
annot1 <- annot1[!duplicated(annot1$Gene.name),]

bool <- (rnas1$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia") & (rnas1$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")
rnas1$grupo[bool] <- "BALL-PB"
bool <- (rnas1$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia") & (rnas1$sample_type == "Primary Blood Derived Cancer - Bone Marrow")
rnas1$grupo[bool] <- "BALL-BM"
bool <- (rnas1$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia") & (rnas1$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow")
rnas1$grupo[bool] <- "BALL-rBM"
bool <- (rnas1$primary_diagnosis == "T lymphoblastic leukemia/lymphoma") & (rnas1$sample_type == "Primary Blood Derived Cancer - Peripheral Blood")
rnas1$grupo[bool] <- "TALL-PB"
bool <- (rnas1$primary_diagnosis == "T lymphoblastic leukemia/lymphoma") & (rnas1$sample_type == "Primary Blood Derived Cancer - Bone Marrow")
rnas1$grupo[bool] <- "TALL-BM"
rnas1$grupo <- as.factor(rnas1$grupo)
rnas1 <- rnas1[,!is.na(rnas1$grupo)]


#Generar QCReport_Before

factorsALL <- rnas.before$grupo
factorsALL <- as.data.frame(factorsALL)
factorsALL <- cbind(factorsALL, rnas.before$sample)
colnames(factorsALL) <- c("grupo", "samples")

#Normalización

ln.data <- withinLaneNormalization(assay(rnas1), annot1$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot1$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(rnas1$grupo))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
assay(rnas1) <- exprs(mydata2corr1)

BALL_PB <- rnas1[, rnas1$grupo == "BALL-PB"]
BALL_BM <- rnas1[, rnas1$grupo == "BALL-BM"]
BALL_rBM <- rnas1[, rnas1$grupo == "BALL-rBM"]
TALL_PB <- rnas1[, rnas1$grupo == "TALL-PB"]
TALL_BM <- rnas1[, rnas1$grupo == "TALL-BM"]

saveRDS(BALL_PB, file="BALL_PB.RDS")
saveRDS(BALL_BM, file="BALL_BM.RDS")
saveRDS(BALL_rBM, file="BALL_rBM.RDS")
saveRDS(TALL_PB, file="tALL_PB.RDS")
saveRDS(TALL_BM, file="TALL_BM.RDS")

mydata2 = readData(assay(rnas1), factors = as.data.frame(rnas1$grupo))

myPCA = dat(mydata2, type = "PCA")
explo.plot(myPCA) 
