#######################################################################
#######         PRE-PROCESSING Proteomics STATegra DATA          #######       
#######################################################################

## By Andreas Schmidt A.Schmidt@med.uni-muenchen.de, David Gomez-Cabrero and Ana Conesa; 

# Raw data are available through ProteomeXchage PXD003263

####  FROM RAW DATA TO PROTEIN QUANTIFICATION   ####
#####################################################

##Proteomics data were acquired in three technical repeats of a 1D reversed phase separation of the 
# complete cellular proteome on an nano-chromatography system coupled to  an LTQ Orbitrap hybrid mass spectrometer.
# Initial protein mapping to a fasta protein sequence database generated from the RNA-seq data of the STATegra project 
# using the MaxQuant database search software (vs 1.5.0.0; https://maxquant.org/). ####

#MaxQuant database search settings

FDR rates
PSM FDR   0.01 
Protein FDR	0.02
Site FDR 	0.05

Result filters
Min. peptide length	7
Min. peptide score	15
Min. # of razor peptides	1

# Variable modifications
# Acetylation of protein N-terminus
# Oxidation of methionine

# Database search was performed against a concatenated forward/reversed library to determine the score distribution 
# of reversed peptide sequences (true negatives) serving as a model for FDR filtering. 
# Further parameters can be found in the corresponding parameter.txt and summary.txt file included in the repository.

# Protein quantitation is part of the search algorithm which extracts the chromatographic peak area of each peptide signal. 
# Peptide peaks of identified sequences are summed up to obtain the protein intensity values. 
# Intensity values were automatically converted to iBAQ values by dividing by the number of expected tryptic peptides for the 
# corresponding protein amino acid sequence. 
# Protein results are written into the protein.groups.txt file including quantitation results after FDR filtering. 


#??The results are written into the Proteomics_01_uniprot_canonical_normalized.txt


#### THIS PORTION OF THE SCRIPT IS RUN IN R USING A .R SCRIPT ####
##################################################################

folder_files <- "/path/to/data"
library(NOISeq)

#### STEP 1: LOAD THE DATA
#### 

proteomics= read.delim(paste(folder_files, "Proteomics_01_uniprot_canonical_normalized.txt",sep=""), as.is = TRUE, header = TRUE, sep = "\t"); dim(proteomics) # 2191 dim(proteomics)
colnames(proteomics)
data = proteomics[1:2527,51:86] ; dim (data)# 2527 36  # extract true values
names = paste("prot", formatC(1:nrow(data), digits = 3, flag = 0), sep = "_")
rownames(data) = names
all.names <- proteomics[1:2527,1:7]
rownames(all.names) = names
colnames(data) = lapply(colnames(data), function (x) strsplit(x, ".", fixed = TRUE)[[1]][3])
colnames(data) = lapply(colnames(data), function (x) substr(x, 4, nchar(x)))
head(data)

#### STEP 2:  DESIGN
####

data = data[,c(2*(1:18)-1, 2*(1:18))] 
colnames(data) = sapply(colnames(data), function (x) sub( "batch", "",x))
design = data.frame("treat" = sapply(colnames(data), function (x) substr(x, 1, 3)),
                     "time" = as.numeric(sapply(sapply(colnames(data), function (x) strsplit(x, "_")[[1]][2]),
                                                function (x) strsplit(x, "h")[[1]][1])),
                     "Brep" = sapply(colnames(data), function (x) substr(x, nchar(x), nchar(x))))
design = data.frame(design, "cond" = apply(design[,1:2], 1, paste, collapse = ""))
design


#### STEP 3: LOG-TRANSFORMATION
####

## Log transformation
data.norm = log2(data)
head(data.norm)


#### STEP 4: MISSING VALUES: DISCARD AND IDENTIFICATION
####

## 4.1 Replace 0 with NA
data.norm[data == 0] <- NA  

## 4.2 Identification of proteins with no missing values
noNA = na.omit(data.norm)
nrow(noNA)  # 612 proteins without any missing

## 4.3 Characterisation of missing values
# Computing a matrix of how many NAs we have per condition at each gene

NAstat = matrix(NA, ncol = 12, nrow = nrow(data.norm))
ncond = unique(design[,"cond"])
for (i in 1:nrow(NAstat)) {
  for (j in 1:length(ncond)) {
    NAstat[i,j] = sum(is.na(data.norm[i,grep(ncond[j], design[,"cond"])]))
  }
}
colnames(NAstat) = ncond
rownames(NAstat) = rownames(data.norm)
head(NAstat)

## 4.4 Filtering based on  criteria: 3 misings per group in at least 11 group
discard = which(apply(NAstat, 1, function (x) length(which(x>2))>10))  # 3 misings in at least 11 samples
data.norm2 = data.norm[-discard,] ;dim (data.norm2 ) #  2396 proteins remaining with at least 1 value in 10 conditions
NAstat2 = NAstat[-discard,]


#### STEP 5: IMPUTATION
####

## 5.1: Imputation of all NA in a given condition:
	#When IK is NA for all samples, the value is replaced with 0 + a "random value below discovery" 
	#When Control is NA for all samples, the value is replaced with 0 + a "random value below discovery" 

std.distr = apply(data.norm2,1, function (x) {tapply(x, design[,"cond"], sd, na.rm = T)} )# sd within conditions
allNACtr = allNATr  = NULL
for (i in 1 : nrow (data.norm2)) {
  if (all(NAstat2[i, 1:6] == 3)) { # if all controls are NA for all 3 replicates
    allNACtr <- c(allNACtr, i) # keep gene number
    distr = rnorm (180,0,median(std.distr, na.rm = T))  # 180 values around 0 with std our std
    data.norm2 [i,1:18] <- distr[which(distr>0)][1:18] # assign positive values
  }
  if (all(NAstat2[i, 7:12] == 3)) { # the same is done for the Treatment
    allNATr <- c(allNATr, i)
    distr= rnorm (180,0,median(std.distr, na.rm = T)) 
    data.norm2 [i,19:36] <- distr[which(distr>0)][1:18] 
  }
}

length(allNACtr) # 21 prots where all Control values were NA and now have values around 0
length(allNATr) # 36 prots where all Control values were NA and nowhave values around 0

## 5.2: Recomputing NA:

NAstat3 = matrix(NA, ncol = 12, nrow = nrow(data.norm2))
ncond = unique(design[,"cond"])

for (i in 1:nrow(NAstat3)) {
  for (j in 1:length(ncond)) {
    NAstat3[i,j] = sum(is.na(data.norm2[i,grep(ncond[j], design[,"cond"])]))
  }
}
colnames(NAstat3) = ncond
rownames(NAstat3) = rownames(data.norm2)
head (NAstat3)

## Identification of proteins with no missing values
noNA = na.omit(data.norm2)
nrow(noNA)  # still 612 proteins without any missing



## 5.3: Imputing for those conditions with only 1 NA: by the mean of the condition.

data.norm3 = data.norm2

for (i in 1:nrow (data.norm3)) {
  if (any(NAstat3[i,] == 1) ) {  
    conds = names(which(NAstat3[i,] == 1))
    for (j in 1:length(conds)){
      jcond = conds[j]
      data.cond = rownames(design[design[,4] == jcond,])
      mydata = data.norm3[i,data.cond]
      time = design[design[,4] == jcond,2][1]
      matching.cond = setdiff( unique(design[design[,2] == time,4]), jcond)
      matching.cond = rownames(design[design[,4] == matching.cond,])
      matching.data = data.norm3[i,matching.cond]
      mysd = sd(mydata, na.rm = T)
      input = runif(1,min(mydata, na.rm = T),max(mydata, na.rm = t))
      value = mean(as.numeric(mydata/matching.data), na.rm = T)
      value = value*matching.data[is.na(mydata)]
      input2 = rnorm (1,value,0.01) 
      remplace = which(is.na(mydata))
      if (all(is.na(matching.data))) {
        data.norm3[i,data.cond[remplace]] = input # if all NAs in matching data, then take an average value
      } else if (all(!is.na(matching.data))) { 
        data.norm3[i,data.cond[remplace]] = input2 # if no NAs in matching data, then make use of ratio info
      }  else { 
        data.norm3[i,data.cond[remplace]] = mean (c(input, input2), na.rm = T)
      } # if 2 NAs in matching data, then take average
    }
  }
}

noNA3 = na.omit(data.norm3); nrow(noNA3)  # 864  # we increased around 200 proteins now having complete values
boxplot (as.matrix(noNA3)~col(noNA3))
data.norm3 <- cbind(all.names[rownames(data.norm3),], data.norm3) # add all protein names

#### STEP 6: TMM normalisation of imputed data

boxplot(as.matrix(noNA3)~col(noNA3)) # distributions are different, normalisation across samples is required.
noNA3_tmmNorm = tmm(noNA3)
boxplot(as.matrix(noNA3_tmmNorm)~col(noNA3_tmmNorm)) # looks much better!
noNA3_tmmNorm2 <- cbind(all.names[rownames(noNA3_tmmNorm),], noNA3_tmmNorm) # add all protein names
####

#### STEP 7: OUTCOME

write.table(noNA3_tmmNorm2, "STATegra_Proteomics_imputed_no_missings_tmm.txt", sep= "\t", row.names = F, col.names = T, quote = F)

