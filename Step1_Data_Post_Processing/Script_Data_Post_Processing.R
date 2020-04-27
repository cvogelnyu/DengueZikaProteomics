# Input file for script is Spectronaut BGS_Report output file which cannot be provided on github due to file size (4.4GB).
# We are therefore providing the reduced output file Spectronaut_Output_Zika_BGS_Report_sample instead, which consists of random sample of 1000 rows of the orginal file.
# From line 75 on, the complete file merged_table_input.txt is used instead of the sampled version.

## This part is to clean up unnecessary data (reduce file size)
d <- read.csv("Spectronaut_Output_Zika_BGS_Report_reduced_sample.csv")

## Make unique identifier to get a single value for a fragment in a peptide
gg = c(1,2,7,12,30,42,43,44,46,47,48,49,50,51,54,56,58)
d = d[, gg]

d = d[d$F.ExcludedFromQuantification == "False", ]

## Order the data in decreasing order of fragment ion peak areas
ord = order(d$F.NormalizedPeakArea, decreasing=TRUE)
d = d[ord, ]

## Make sure that there are no replicate fragment ions 
## (by taking maximum intensity of the ions of the same type)
## type meaning: same peptide sequence, precursor ion charge, fragment ion mz and charge,
## fragment ion loss type (e.g. H2O)

d$F.FrgMz = round(d$F.FrgMz, 1)
id = paste(d$R.FileName, d$PG.ProteinGroups, d$PEP.GroupingKey, d$FG.Charge, 
           d$F.FrgIon, d$F.Charge, d$F.FrgLossType)
uid = unique(id)
m = match(uid, id)
d = d[m, ]

id = paste(d$PG.ProteinGroups, d$PEP.GroupingKey, d$FG.Charge, 
           d$F.FrgIon, d$F.Charge, d$F.FrgLossType)
uid = unique(id)
nid = length(uid)

samples = unique(d$R.FileName)
nsamples = length(samples)

PAdat = matrix(NA, nid, nsamples) 

for(j in 1:nsamples) {
  print(j)
  wid = d$R.FileName == samples[j]
  tmpd = d[wid, ]
  tmpid = id[wid]
  m = match(tmpid, uid)
  PAdat[m,j] = tmpd$F.NormalizedPeakArea
}

colnames(PAdat) = samples

Protein = Peptide = PepPrecursorCharge = FragIon = FragCharge = FragMZ = FragLoss = rep(NA, nid)
for(i in 1:nid) {
  tmp = strsplit(uid[i], " ")[[1]]
  Protein[i] = tmp[1]
  Peptide[i] = tmp[2]
  PepPrecursorCharge[i] = tmp[3]
  FragIon[i] = tmp[4]
  FragCharge[i] = tmp[5]
  FragLoss[i] = tmp[6]
  if(i %% 10000 == 0) print(i)
}

dat = data.frame(Protein=Protein, Peptide=Peptide, PrecursorCharge=PepPrecursorCharge, 
                 FragmentIon=FragIon, FragmentCharge=FragCharge, 
                 FragmentLoss=FragLoss, round(PAdat,2), stringsAsFactors = FALSE)

ord = order(dat$Protein, dat$Peptide, dat$PrecursorCharge, dat$FragmentIon, dat$FragmentCharge)
dat = dat[ord, ]

write.table(dat, "merged_table_input_sample.txt", sep="\t", quote=F, row.names=F, na="")


# Quality Control, Filtering

d = read.delim("merged_table_input.txt", header=T, as.is=T)

## Divide the data into fragment ion meta data and quantitative data (tmpdata = QC + samples)
metadata = d[,1:6]
tmpdata = d[,-c(1:6)]

sdata = read.delim("ConditionSetup_meta.txt", header=T, as.is=T)
sdata$File.Name = paste("X", sdata$File.Name, sep="")

m = match(sdata$File.Name, colnames(tmpdata))
tmpdata = tmpdata[,m]

## Remove artefacts from the abundance data
tmpdata[tmpdata <= 50] = NA

write.csv(tmpdata, "data_boxplot1.csv")

pdf("boxplot.pdf", height=10, width=25)
par(mar=c(20,5,2,3))
boxplot(log2(tmpdata), las=2, ylab="Log2 Peak Area", cex=.3)
dev.off()

### At this point, we decide not to use patient 48
rid = c(14,31)
tmpdata = tmpdata[,-rid]
sdata = sdata[-rid, ]

#########################
## Remove missing data and try PCA first
gg = grep("QC", colnames(tmpdata))
wid1 = apply(tmpdata[,gg], 1, function(x) mean(!is.na(x))) >= 0.6
wid2 = apply(tmpdata[,-gg], 1, function(x) mean(!is.na(x))) >= 0.6

tmpdata2 = tmpdata[wid2, ] 
metadata2 = metadata[wid2, ] 

length(unique(metadata$Protein))
#expected number of proteins: 517
length(unique(metadata2$Protein))
#expected number of proteins: 303

### Impute missing data first 
nr = nrow(tmpdata2)
x = as.matrix(log2(tmpdata2))

#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)
x.i = impute.knn(x)$data  ### This object is solely for PCA.

ppp = c(rep(19,78), rep(15,64))

### PCA: let's see what the overall trend looks like first

tmp.pca = prcomp(t(x.i))
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)

sdata$color = 1
sdata$color[grep("Dengue",sdata$Condition)] = 2
sdata$color[grep("Zika",sdata$Condition)] = 4

LM = 150
pdf("PCAplot_before_batch_correction.pdf", height=5.5, width=5, useDingbats = FALSE)
plot(tmp.pca$x[,1], tmp.pca$x[,2], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1", ylab="PC2", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
abline(v=0, lty=2)
abline(h=0, lty=2)
dev.off()

###################################################

### New column (but same type) was used after sample 78
### Use all samples to fit lowess curve (for now fit one curve for all, without separately fitting for two columns)

pt1 = 1:78 #batch1
pt2 = 79:142 #batch2

pdf("fragment_ion_plots_before_batch_correction.pdf", height=12, width=8, useDingbats = FALSE)
par(mfrow=c(4,1))
for(i in 1:nrow(x)) {
  tmp = as.numeric(x[i,])
  tmp1 = tmp[pt1]
  tmp2 = tmp[pt2]
  
  plot(tmp, col=sdata$color, pch=19, cex=.5, main=paste(metadata2$Protein[i], metadata2$Peptide[i]))
  uid1 = !is.na(tmp1)
  uid2 = !is.na(tmp2)
  lines(lowess(pt1[uid1], tmp1[uid1], f=0.25))
  lines(lowess(pt2[uid2], tmp2[uid2], f=0.25))
  abline(v=78.5)
  if(i %% 1000 == 0) print(i)
}
dev.off()


#####################################
### PC2 shows serious batch effect -- need to correct batch effect first
### For each fragment, fit lowess and subtract in each batch
### Then set the average to a common overall average peak area value

sdata$batch = 2
sdata$batch[1:78] = 1   ### old batch

sdata$Group = "QC"
sdata$Group[grep("Dengue",sdata$Condition)] = "Dengue"
sdata$Group[grep("Zika",sdata$Condition)] = "Zika"
sdata$Group = as.factor(sdata$Group)

sdata$color = 1
sdata$color[grep("Dengue",sdata$Condition)] = 2
sdata$color[grep("Zika",sdata$Condition)] = 4

xx = x   ### Create a new object
nr = nrow(x)
for(i in 1:nr) {
  zz = as.numeric(x[i,])
  
  zz1 = zz[pt1]
  zz2 = zz[pt2]
  uid1 = !is.na(zz1)
  uid2 = !is.na(zz2)
  
  ### Batch 1 and 2 lowess fit
  fit1 = lowess(pt1[uid1], zz1[uid1], f=0.25)
  fit2 = lowess(pt2[uid2], zz2[uid2], f=0.25)
  
  ### Equalize the mean
  med1 = median(fit1$y)
  med2 = median(fit2$y)
  mm = (med1 + med2) / 2
  
  zz1[uid1] = zz1[uid1] - fit1$y + mm
  zz2[uid2] = zz2[uid2] - fit2$y + mm
  zz[pt1] = zz1
  zz[pt2] = zz2
  xx[i,] = as.numeric(zz)
  
  if(i %% 1000 == 0) print(i)
}


### Remove outliers: 3 standard deviation
nr = nrow(xx)

for(i in 1:nr) {
  zz = as.numeric(xx[i,])
  mm = mean(zz, na.rm=TRUE)
  ss = sd(zz, na.rm=TRUE)
  wid = abs(zz - mm) > (3*ss)
  wid[is.na(wid)] = TRUE
  zz[wid] = NA
  xx[i,] = zz
}

pdf("fragment_ion_plots_after_batch_correction.pdf", height=12, width=8, useDingbats = FALSE)

par(mfrow=c(4,1))
for(i in 1:nrow(xx)) {
  tmp = as.numeric(xx[i,])
  tmp1 = tmp[pt1]
  tmp2 = tmp[pt2]
  
  plot(tmp, col=sdata$color, pch=19, cex=.5, main=paste(metadata2$Protein[i], metadata2$Peptide[i]))
  uid1 = !is.na(tmp1)
  uid2 = !is.na(tmp2)
  lines(lowess(pt1[uid1], tmp1[uid1], f=0.25))
  lines(lowess(pt2[uid2], tmp2[uid2], f=0.25))
  abline(v=78.5)
  if(i %% 1000 == 0) print(i)
}

dev.off()


#### Another round of PCA
library(impute)
xx.i = impute.knn(xx)$data

tmp.pca = prcomp(t(xx.i))
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)

LM = 130
pdf("PCAplot_after_batch_correction.pdf", height=5.5, width=5, useDingbats = FALSE)
plot(tmp.pca$x[,1], tmp.pca$x[,2], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1", ylab="PC2", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)
dev.off()

tmpdata = 2^xx ## data without imputation!

write.csv(tmpdata, "data_boxplot2.csv")

####Second boxplot
pdf("boxplot_after_batch_correction.pdf", height=10, width=25)
par(mar=c(20,5,2,3))
boxplot(log2(tmpdata), las=2, ylab="Log2 Peak Area", cex=.3)
dev.off()

### CoV: computed in log2 scale
qcID = grep("QC", colnames(tmpdata))
tmpdata2 = tmpdata[,qcID]
sd1 = apply(tmpdata2, 1, function(x) mad(x, na.rm=TRUE))   
mean1 = apply(tmpdata2, 1, function(x) median(x, na.rm=TRUE))  
CoV = sd1 / mean1 * 100

pdf("CoV_plot.pdf", height=5.5, width=5)
hist(CoV[CoV < 200], breaks=100, xlim=c(0,100), xlab="Coefficient of variation (%)", main="Variability of QCs")
dev.off()

write.csv(CoV[CoV < 200], "cov_data.csv")

#########################################
## Another round of PCA

tmp.pca = prcomp(t(xx.i))
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)

### PCA summary 
vv
#expected results for vv (first 6 values): 19.64  6.39  4.49  3.45  2.88  2.77


pdf("PCAplot_summary.pdf", height=11, width=10.5, useDingbats = FALSE)
LM = 70

par(mfrow=c(3,3))

### First row
plot(tmp.pca$x[,1], tmp.pca$x[,1], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1", ylab="PC1", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,1], tmp.pca$x[,2], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1", ylab="PC2", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,1], tmp.pca$x[,3], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC1", ylab="PC3", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

### Second row
plot(tmp.pca$x[,2], tmp.pca$x[,1], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC2", ylab="PC1", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,2], tmp.pca$x[,2], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC2", ylab="PC2", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,2], tmp.pca$x[,3], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC2", ylab="PC3", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

### Third row
plot(tmp.pca$x[,3], tmp.pca$x[,1], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC3", ylab="PC1", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,3], tmp.pca$x[,2], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC3", ylab="PC2", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

plot(tmp.pca$x[,3], tmp.pca$x[,3], col=sdata$color, pch=ppp,
     xlim=c(-LM,LM), ylim=c(-LM,LM), xlab="PC3", ylab="PC3", 
     main="Fragment ion peak area data", cex=0.5)
legend(-LM,LM, c("QC-1","Dengue-1","Zika-1", "QC-2","Dengue-2","Zika-2"), pch=c(19,19,19,15,15,15), col=c(1,2,4), ncol=2, cex=0.8)
#text(tmp.pca$x[,1], tmp.pca$x[,2], colnames(xx), cex=.4)
abline(v=0, lty=2)
abline(h=0, lty=2)

dev.off()


#######################################
### Export the data

### First sort the columns from Dengue to Zika
### and order each person's replicate side by side

ord = order(sdata[["Diagnosed_Type"]], sdata[["PatientID"]], sdata[["Sample_TimePoint"]]) 
sdata = sdata[ord, ]
m = match(sdata$File.Name, colnames(tmpdata))

tmpdata = tmpdata[,m]

metadata.new = metadata2 
qcID = grep("QC", colnames(tmpdata))
tmpdata.new = round(tmpdata[, -qcID], 4) 
tmpdataQC.new = round(tmpdata[, qcID], 4) 

tmpdata.new = data.frame(metadata.new, tmpdata.new, stringsAsFactors = FALSE)
tmpdataQC.new = data.frame(metadata.new, tmpdataQC.new, stringsAsFactors = FALSE)

dir.create("normalized_QCed")
write.table(metadata.new, "normalized_QCed/meta_data.txt", sep="\t", quote=F, row.names=F, na="")
write.table(tmpdata.new, "normalized_QCed/quant_data.txt", sep="\t", quote=F, row.names=F, na="")
write.table(tmpdataQC.new, "normalized_QCed/quantQC_data.txt", sep="\t", quote=F, row.names=F, na="")

### How many proteins
prot = metadata.new$Protein
pep = metadata.new$Peptide
frag = metadata.new$FragmentIon

length(unique(prot))
length(unique(paste(prot, pep)))
length(unique(paste(prot, pep, frag)))

#Expected results: 303 proteins, 3329 peptides, 12887 fragments



