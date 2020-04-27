#install.packages("gplots")
#install.packages("pheatmap")
#install.packages("quantable")

# normalization
library("quantable")
# heatmap
library("pheatmap")
# for color bluered() in heatmap
library("gplots")

pdata = read.delim("protein_level_data_mapDIA.txt", header=T, as.is=T, row.names = 1)
sdata = read.delim("sample_meta_data.txt", header=T) #, as.is=T)  ## We need to keep categorial vars as "factors"


## Protein level data (remove the last two columns -- protein & peptide info)
# log2, row wise mean centering
pdata[pdata == 0] = NA
nc = ncol(pdata)
pdata = pdata[,-c(nc-1,nc)]
pdata = log2(pdata + 2000)
mp = apply(pdata, 1, function(x) mean(x, na.rm=TRUE))
pdata = sweep(pdata, 1, mp)

#Split pdata into training set and test set. 
#Samples for test set have been randomly selected, but were not allowed to overlap with samples needed for analysis B (47 samples with complete meta data) in order to not decrease the number of samples for this set further.

pdata_testsamples <- read.csv("Test_Samples.csv", header = T)
m = match(pdata_testsamples$Test_Samples, colnames(pdata))
pdata_testset <- pdata[ ,m]
pdata_trainingset <- pdata[ ,-m]

#write.csv(pdata_trainingset, "Protein_data_mapDIA_norm_Trainingsset.csv")
#write.csv(pdata_testset, "Protein_data_mapDIA_norm_Testset.csv")

pdata_a1 <- pdata_trainingset #set for analysis A
pdata_a2 <- pdata_trainingset #set for analysis B


## Remove samples with incomplete phenotype data
## Error for IgG: So we leave IgG out: Coefficients: (1 not defined because of singularities)
sdata$File.Name = paste("X", sdata$File.Name, sep="")
m = match(colnames(pdata_trainingset), sdata$File.Name)
sdata_trainingset = sdata[m,]

attach(sdata_trainingset)
id1 = !is.na(Diagnosed_Type)
id2 = !is.na(Diagnosed_Type) & !is.na(IgM) & !is.na(Gender) & !is.na(Race) & !is.na(Days_since_onset_corrected) & !is.na(Age)
detach(sdata_trainingset)

sdata_a1 = sdata_trainingset[id1, ]
sdata_a2 = sdata_trainingset[id2, ]


#Analysis_A (no confounding factors)

attach(sdata_a1)

m = match(sdata_a1$File.Name, colnames(pdata_a1)) #get rid off samples with incomplete data
pdata_a1 = pdata_a1[ ,m]

np = nrow(pdata_a1)

pval.zika = rep(NA, np)
coef.zika = rep(NA, np)

resdata = pdata_a1

for(i in 1:np) {
  
  xx = as.numeric(pdata_a1[i,])
  tmp.fit = lm(xx ~ Diagnosed_Type) 
  tmp.coef = coef(tmp.fit)
  tmp.res = xx  * tmp.coef[2]
  resdata[i,] = tmp.res
  
  ## record p-values
  pval.zika[i] = summary(tmp.fit)$coef[2,4]
  
  
  coef.zika[i] = summary(tmp.fit)$coef[2,1]
}

pmat = data.frame(Protein=rownames(pdata_a1), 
                  pval_Zika=pval.zika, 
                  coef_Zika=coef.zika, 
                  stringsAsFactors = FALSE)

pmat[,2] = round(pmat[,2], 5)
pmat[,3] = round(pmat[,3], 2)

write.table(pmat, "pvalue_matrix_AnalysisA_Diagnosis_Training.txt", sep="\t", quote=F, row.names=F)

## Impute missing data for visualization and clustering
library(impute)
pdata2_a1 = impute.knn(as.matrix(resdata))$data   
### pdata2 object holds the "residual" data after removing confounder effects
write.csv(pdata2_a1, "Residual_data_AnalysisA.csv")

#### Significant proteins only
psig = pmat$pval_Zika <= 0.01
pdata3_a1 = pdata2_a1[psig, ]

detach(sdata_a1)



#Analysis B (remove confounding factors)
attach(sdata_a2)

m = match(sdata_a2$File.Name, colnames(pdata_a2)) #get rid off samples with incomplete data
pdata_a2 = pdata_a2[ ,m]

np = nrow(pdata_a2)

pval.zika = pval.IgM = pval.Gender = pval.Race = pval.Days_since_onset_corrected = pval.Age = rep(NA, np)
coef.zika = coef.IgM = coef.Gender = coef.Race = coef.Days_since_onset_corrected = coef.Age = rep(NA, np)

mm = median(Days_since_onset_corrected, na.rm = TRUE)
nn = median(Age, na.rm = TRUE)

Days_since_onset_corrected = Days_since_onset_corrected - mm
Age = Age - nn

resdata = pdata_a2

for(i in 1:np) {
  
  xx = as.numeric(pdata_a2[i,])
  tmp.fit = lm(xx ~ Diagnosed_Type + IgM + Gender + Race + Days_since_onset_corrected + Age) 
  tmp.coef = coef(tmp.fit)
  tmp.res = xx - IgM * tmp.coef[3]
  tmp.res = tmp.res - ifelse(Gender=="Male", 1, 0) * tmp.coef[4]
  tmp.res = tmp.res - ifelse(Race=="Indo", 1, 0) * tmp.coef[5]
  tmp.res = tmp.res - Days_since_onset_corrected * tmp.coef[6]
  tmp.res = tmp.res - Age * tmp.coef[7]
  resdata[i,] = tmp.res
  
  ## record p-values
  pval.zika[i] = summary(tmp.fit)$coef[2,4]
  pval.IgM[i] = summary(tmp.fit)$coef[3,4]
  pval.Gender[i] = summary(tmp.fit)$coef[4,4]
  pval.Race[i] = summary(tmp.fit)$coef[5,4]
  pval.Days_since_onset_corrected[i] = summary(tmp.fit)$coef[6,4]
  pval.Age[i] = summary(tmp.fit)$coef[7,4]
  
  coef.zika[i] = summary(tmp.fit)$coef[2,1]
  coef.IgM[i] = summary(tmp.fit)$coef[3,1]
  coef.Gender[i] = summary(tmp.fit)$coef[4,1]
  coef.Race[i] = summary(tmp.fit)$coef[5,1]
  coef.Days_since_onset_corrected[i] = summary(tmp.fit)$coef[6,1]
  coef.Age[i] = summary(tmp.fit)$coef[7,1]
}

pmat = data.frame(Protein=rownames(pdata_a2), 
                  pval_Zika=pval.zika, pval_IgM=pval.IgM, 
                  pval_GenderMale=pval.Gender, pval_RaceIndo=pval.Race, pval_Days_since_onset_corrected=pval.Days_since_onset_corrected, pval_Age=pval.Age,
                  coef_Zika=coef.zika, coef_IgM=coef.IgM, 
                  coef_GenderMale=coef.Gender, coef_RaceIndo=coef.Race, coef_Days_since_onset_corrected=coef.Days_since_onset_corrected, coef_Age = coef.Age,
                  stringsAsFactors = FALSE)

pmat[,2:7] = round(pmat[,2:7], 5)
pmat[,8:13] = round(pmat[,8:13], 2)

write.table(pmat, "pvalue_matrix_AnalysisB_Diagnosis_Training.txt", sep="\t", quote=F, row.names=F)

## Impute missing data for visualization and clustering
library(impute)
pdata2_a2 = impute.knn(as.matrix(resdata))$data   
### pdata2 object holds the "residual" data after removing confounder effects
write.csv(pdata2_a2, "Residual_data_AnalysisB.csv")

#### Significant proteins only
psig = pmat$pval_Zika <= 0.01
pdata3_a2 = pdata2_a2[psig, ]

detach(sdata_a2)
