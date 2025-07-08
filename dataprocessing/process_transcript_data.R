##############################
### LOAD PACKAGES
##############################

library(biomaRt)
library(NOISeq)
library(edgeR)
library(cqn)
library(scales)
library(limma)
library(maSigPro)
library(sva)

##############################
### DATA: count tables and DESIGN
##############################

## Count data per samples
M.intersect.36 <- read.csv("../datasets/STATegra/Script_STATegra_RNAseq/RNA_seq/STATegra.RNAseq.allSamples.counts.csv",row.names=1)
#head(M.intersect.36[,1:4],3)
#ExpBatch_4_Ctr_0H ExpBatch_5_Ctr_0H ExpBatch_6_Ctr_0H ExpBatch_1_Ctr_2H
#ENSMUSG00000000001              3336              3584              6999              3565
#ENSMUSG00000000028              5801              7374             10889              6588
#ENSMUSG00000000031                33                51                50                17

## Design, organized as the 
DESIGN.36 <- read.csv("../datasets/STATegra/Script_STATegra_RNAseq/RNA_seq/DESIGN.36.csv",row.names = 1)
#head(DESIGN.36)
#EXP_BATCH LIBRARY_PREP HiSeq.RUN QubitCon nMQubit Kappa_Quant TIME IKAROS ERC.spikes
#ExpBatch_4_Ctr_0H         4      BATCH_6         3     9.14   49.14       70.49    0      0      0.971
#ExpBatch_5_Ctr_0H         5      BATCH_1         1    17.10   91.94      203.88    0      0      0.972
#ExpBatch_6_Ctr_0H         6      BATCH_4         1    12.70   68.28      121.57    0      0      0.971

## Annotation downloaded from BiomaRt during 2015: Gene name + Length + GC content
gene.INFO<-read.csv("../datasets/STATegra/Script_STATegra_RNAseq/RNA_seq/ENS_MOUSE_LENGHT_GC.txt",sep="\t",row.names=1)
#head(gene.INFO)
#GC LENGTH
#ENSMUSG00000095309 42.75    923
#ENSMUSG00000000126 52.53   6630
#ENSMUSG00000086196 42.02   3496

## Checking all sample same organization
rownames(DESIGN.36)==colnames(M.intersect.36)

##############################
### STEP 1, FILTER: 
##############################

# Decision criteria: "at least 10 reads in at least 20 samples"
# Result differ by <10 genes with cpm>1 in at least 3 sampels (minimal group size)

M.intersect.n10<-apply(M.intersect.36,1,function(x){sum(x>10)})
M.intersect.36F<-M.intersect.36[M.intersect.n10>=20,]

M.go<-M.intersect.36F
M.go<-M.go[rownames(M.go) %in% rownames(gene.INFO),] 
gene.INFO.GO<-gene.INFO[rownames(gene.INFO) %in% rownames(M.go),]

##############################
### STEP 2, NORMALIZATION: "at least 10 reads in at least 20 samples"
##############################

#### ORDER THE VARIABLES & DESIGN
M.go.ORD<-M.go[order(rownames(M.go)),]
gene.INFO.GO.ORD<-gene.INFO.GO[order(rownames(gene.INFO.GO)),]
sizeFactors.GO<-apply(M.go.ORD,2,sum)

DESIGN.ORD<-DESIGN.36[colnames(M.go.ORD),]
rownames(DESIGN.ORD)==colnames(M.go.ORD)

Time <- DESIGN.ORD$TIME
IKAROS<- DESIGN.ORD$IKAROS
BATCH<-DESIGN.ORD$EXP_BATCH
LIB<-DESIGN.ORD$LIBRARY_PREP

## We consider 12 groups: 6 times x 2 conditions
TIME_IK<-paste("T",Time,"_","Exp",IKAROS,sep="")

#### NORMALIZATION: CQN

cqn.subset <- cqn(M.go.ORD, lengths = gene.INFO.GO.ORD$LENGTH,
                  x = gene.INFO.GO.ORD$GC, sizeFactors = sizeFactors.GO,
                  verbose = TRUE)
M.go.ORD.cqn<- cqn.subset$y + cqn.subset$offset
head(M.go.ORD.cqn)

#### BATCH CORRECTION: COMBAT

design.cb<-model.matrix(~ TIME_IK )
quant.batch.adj.CQN.LIB<-ComBat(M.go.ORD.cqn, batch=LIB, mod=design.cb) 
n.sv.CQN = num.sv(M.go.ORD.cqn,design.cb,method="leek")
n.sv.CORRECTED.CQN = num.sv(quant.batch.adj.CQN.LIB,design.cb,method="leek")

#### Correction to set minimun value to 1 (a linear correction)
DE.set2<-quant.batch.adj.CQN.LIB+abs(min(quant.batch.adj.CQN.LIB)) + 1

#### Saving the file
write.csv(DE.set2,"../datasets/STATegra/Script_STATegra_RNAseq/RNA_seq/STATegra.RNAseq.CQN.Combat.Annotated.positive_2014_09.csv")