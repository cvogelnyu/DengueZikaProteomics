pdata = read.delim("protein_level_data_mapDIA.txt", header=T, as.is=T, row.names = 1)

## Protein level data (remove the last two columns)
# log2, row wise mean centering
pdata[pdata == 0] = NA
nc = ncol(pdata)
pdata = pdata[,-c(nc-1,nc)]
pdata = log2(pdata)
mm = apply(pdata, 1, function(x) mean(x,na.rm=TRUE))
pdata = sweep(pdata, 1, mm)

## Impute missing data for clustering
library(impute)
pdata.i = impute.knn(as.matrix(pdata))$data

#################################
### Protein-level PCA

## use all proteins
tmp.pca = prcomp(t(pdata.i))
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)

ccc = rep('#31a354', ncol(pdata.i))
ccc[1:68] = '#fec44f'

vv
#PC percentages: 
# [1] 21.40  8.03  6.24  4.10  3.09  2.84  2.27  2.05  1.84  1.70  1.63
# [12]  1.61  1.55  1.50  1.41  1.36  1.25  1.24  1.14  1.12  1.05  1.02

pdf("PCAplot_protein_level.pdf", height=8.5, width=8, useDingbats = FALSE)
LM = 20

par(mfrow=c(2,2))

plot(tmp.pca$x[,1], tmp.pca$x[,2], col=ccc, pch=19,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1 (21%)", ylab="PC2 (8%)", 
     main="Protein abundance data", cex=0.5)
legend(-LM,LM, c("Dengue","Zika"), pch=19, col=c('#fec44f', '#31a354'), cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,1], tmp.pca$x[,3], col=ccc, pch=19,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1 (21%)", ylab="PC3 (6%)", 
     main="Protein abundance data", cex=0.5)
legend(-LM,LM, c("Dengue","Zika"), pch=19, col=c('#fec44f', '#31a354'), cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

LM = 10
plot(tmp.pca$x[,2], tmp.pca$x[,3], col=ccc, pch=19,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC2 (8%)", ylab="PC3 (6%)", 
     main="Protein abundance data", cex=0.5)
legend(-LM,LM, c("Dengue","Zika"), pch=19, col=c('#fec44f', '#31a354'), cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

LM = 10
plot(tmp.pca$x[,3], tmp.pca$x[,4], col=ccc, pch=19,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC3 (6%)", ylab="PC4 (4%)", 
     main="Protein abundance data", cex=0.5)
legend(-LM,LM, c("Dengue","Zika"), pch=19, col=c('#fec44f', '#31a354'), cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

dev.off()


#####################################
### Analysis of PCs: what drives the most of the variations? (clinical parameters)
scores = tmp.pca$x[,1:5]  ### look at first 5 PC scores

npc = 5
tmp = data.frame(Sample=rownames(scores), scores[,1:npc])
write.table(tmp, "PCscores_Protein_Level.txt", sep="\t", quote=F, row.names=F)

phenodata = read.delim("sample_meta_data.txt", header=T, as.is=T)
phenodata$File.Name = paste("X", phenodata$File.Name, sep="")

vars = colnames(phenodata)[c(8,10,12:22,33)]
nvar = length(vars)
PCsig = matrix(NA, npc, nvar)
colnames(PCsig) = vars
rownames(PCsig) = paste("PC", 1:npc, sep="")

qcID = grep("QC", phenodata$File.Name)
phenodata = phenodata[-qcID, ]

badsample = grep("patient_48", phenodata$File.Name)
phenodata = phenodata[-badsample, ]

m = match(phenodata$File.Name, rownames(scores))
scores = scores[m, ]

phenodata$PC1 = scores[, 1]
phenodata$PC2 = scores[, 2]
phenodata$PC3 = scores[, 3]
phenodata$PC4 = scores[, 4]
phenodata$PC5 = scores[, 5]

write.table(phenodata, "meta_data_with_PC_protein_level.txt", sep="\t", quote=F, row.names=F, na="")

#### Diagnosis
varType = rep("Discrete", nvar)
varType[vars %in% c("Days_since_onset","Age", "Days_since_onset_corrected")] = "Continuous"

for(j in 1:nvar) {
  m = match(vars[j], colnames(phenodata))
  uid = !is.na(phenodata[,m])
  
  ### Discrete variable
  if(varType[j] == "Discrete") {
    for(i in 1:npc) {
      cid = grep(paste("PC", i, sep=""), colnames(phenodata))
      tmp.fit1 = lm(phenodata[uid,cid] ~ 1)
      tmp.fit2 = lm(phenodata[uid,cid] ~ as.factor(phenodata[uid,m]))
      PCsig[i,j] = anova(tmp.fit1, tmp.fit2)$'Pr(>F)'[2]
      PCsig[i,j] = round(PCsig[i,j], 4)
    }
  }
  
  ### Continuous variable
  if(varType[j] == "Continuous") {
    for(i in 1:npc) {
      cid = grep(paste("PC", i, sep=""), colnames(phenodata))
      tmp.fit = lm(phenodata[uid,cid] ~ phenodata[uid,m])
      PCsig[i,j] = summary(tmp.fit)$coefficient[2,4]
      PCsig[i,j] = round(PCsig[i,j], 4)
    }
  }
  
}

res = data.frame(PCname = rownames(PCsig), PCsig, stringsAsFactors = FALSE)
write.table(res, "PrincipalComponent_Protein_Level_ClinicalVar.txt", sep="\t", quote=F, row.names=F)


