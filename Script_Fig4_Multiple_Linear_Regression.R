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
sdata = read.delim("meta_data_with_PC.txt", header=T) #, as.is=T)  ## We need to keep categorial vars as "factors"

## Remove samples with incomplete phenotype data
## Error for IgG: So we leave IgG out: Coefficients: (1 not defined because of singularities)
sdata = sdata[,c(1:33)]
attach(sdata)

id = !is.na(Diagnosed_Type) & !is.na(IgM) & !is.na(Gender) & !is.na(Race) & !is.na(Days_since_onset_corrected) & !is.na(Age)
detach(sdata)

sdata = sdata[id, ]

## Sort by virus and individual
ord = order(sdata$Diagnosed_Type, sdata$PatientID)
sdata = sdata[ord, ]

attach(sdata)

## Protein level data (remove the last two columns -- protein & peptide info)
# log2, row wise mean centering
pdata[pdata == 0] = NA
nc = ncol(pdata)
pdata = pdata[,-c(nc-1,nc)]
pdata = log2(pdata + 2000)
mp = apply(pdata, 1, function(x) mean(x, na.rm=TRUE))
pdata = sweep(pdata, 1, mp)

pdata_all_patients <- pdata
m2 = match(sdata_all_patients$File.Name, colnames(pdata_all_patients)) 
pdata_all_patients = pdata_all_patients[ ,m2]

m = match(sdata$File.Name, colnames(pdata)) #get rid off samples with incomplete data
pdata = pdata[ ,m]

np = nrow(pdata)

pval.zika = pval.IgM = pval.Gender = pval.Race = pval.Days_since_onset_corrected = pval.Age = rep(NA, np)
coef.zika = coef.IgM = coef.Gender = coef.Race = coef.Days_since_onset_corrected = coef.Age = rep(NA, np)

mm = median(Days_since_onset_corrected, na.rm = TRUE)
nn = median(Age, na.rm = TRUE)

Days_since_onset_corrected = Days_since_onset_corrected - mm
Age = Age - nn

resdata = pdata

for(i in 1:np) {
  
  xx = as.numeric(pdata[i,])
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

pmat = data.frame(Protein=rownames(pdata), 
                  pval_Zika=pval.zika, pval_IgM=pval.IgM, 
                  pval_GenderMale=pval.Gender, pval_RaceIndo=pval.Race, pval_Days_since_onset_corrected=pval.Days_since_onset_corrected, pval_Age=pval.Age,
                  coef_Zika=coef.zika, coef_IgM=coef.IgM, 
                  coef_GenderMale=coef.Gender, coef_RaceIndo=coef.Race, coef_Days_since_onset_corrected=coef.Days_since_onset_corrected, coef_Age = coef.Age,
                  stringsAsFactors = FALSE)

pmat[,2:7] = round(pmat[,2:7], 5)
pmat[,8:13] = round(pmat[,8:13], 2)

write.table(pmat, "pvalue_matrix.txt", sep="\t", quote=F, row.names=F)

## Impute missing data for visualization and clustering
library(impute)
pdata2 = impute.knn(as.matrix(resdata))$data   
### pdata2 object holds the "residual" data after removing confounder effects
#### Significant proteins only
psig = pmat$pval_Zika <= 0.01
pdata3 = pdata2[psig, ]

## Fig 4 with pheatmap

##residual data, 26 significant proteins, patients with complete metadata######
#export pdata3 to add gene names
write.table(pdata3, "pdata3.txt", sep="\t", quote=F)
#read in table with gene names: pdata3__gene_names.xlsx
#set first column as rownames
pdata3_gene_names <- data.frame(X4_pdata3_190509_gene_names, row.names = 1)


# data
data_sig <- pdata3_gene_names
colnames(data_sig) <- gsub(".*patient", "patient", colnames(data_sig))

#annotations
annotations <- sdata

annotations$File.Name <- gsub(".*patient", "patient", annotations$File.Name)
annotations$File.Name <- gsub("patient", "p", annotations$File.Name)
annotations$File.Name <- gsub("Dengue", "D", annotations$File.Name)
annotations$File.Name <- gsub("Zika", "Z", annotations$File.Name)

data_sig_mat <- as.matrix(data_sig)

colnames(data_sig_mat) <- gsub(".*patient", "patient", colnames(data_sig_mat))
colnames(data_sig_mat) <- gsub("patient", "p", colnames(data_sig_mat))
colnames(data_sig_mat) <- gsub("Dengue", "D", colnames(data_sig_mat))
colnames(data_sig_mat) <- gsub("Zika", "Z", colnames(data_sig_mat))

pheatmap(data_sig_mat,
         main = "", fontsize = 7,
         color = bluered(20),
         show_rownames=T,
         scale = "row",
         breaks=seq(-2,2,by=0.2),
         na_col="grey",
         cluster_cols=F,
         #cluster_rows=F,
         annotation_legend = F,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,
         cutree_cols = 1,
         annotation_col = data.frame(row.names = annotations[,2], 
                                     Diagnosis = annotations[,8]),
         annotation_names_col = F,
         annotation_colors = list(Diagnosis = c(`Dengue` = '#fec44f', `Zika` = '#31a354')))


####same for p</=0.05 (Supplement Figure)
##residual data, 71 significant proteins, patients with complete metadata######
#export pdata4 to add gene names
write.table(pdata4, "4_pdata4_190610.txt", sep="\t", quote=F)
#read in table with gene names: 4_pdata4_190610_gene_names.xlsx
#set first column as rownames
pdata4_gene_names <- data.frame(X4_pdata4_190610_gene_names, row.names = 1)

data_sig0.05 <- pdata4_gene_names
data_sig_mat0.05 <- as.matrix(data_sig0.05)

pheatmap(data_sig_mat0.05,
         main = "", fontsize = 6,
         color = bluered(20),
         show_rownames=T,
         scale = "row",
         breaks=seq(-2,2,by=0.2),
         na_col="grey",
         cluster_cols=F,
         #cluster_rows=F,
         annotation_legend = F,
         border_color = NA,
         treeheight_row = 0, treeheight_col = 0,
         cutree_cols = 1, 
         annotation_col = data.frame(row.names = annotations[,2], 
                                     Diagnosis = annotations[,8]),
         annotation_names_col = F,
         annotation_colors = list(Diagnosis = c(`Dengue` = '#fec44f', `Zika` = '#31a354')))
