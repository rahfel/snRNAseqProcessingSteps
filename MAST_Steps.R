
library(colorout)
library(DropletUtils)
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(sctransform)
library(data.table)
library(conos)
library(pagoda2)
library(MAST)
library(edgeR)
library(graphics)
library(gplots)

tmp = list.files(pattern="*OE.NonNormalised.filteredcounts.dgCMatrix.rds")
Allfiles = lapply(tmp,readRDS)
filenames=list.files(pattern="*OE.NonNormalised.filteredcounts.dgCMatrix.rds", full.names=TRUE)
names(Allfiles) = gsub(".*/", "", filenames)
names(Allfiles) = gsub(".OE.NonNormalised.filteredcounts.dgCMatrix.rds", "", names(Allfiles))
############################################################################################
############################################################################################
######################## METADATA ##########################################
######### CONTROL AND CASES ####################################################
MetaData.C_PD=read.csv("MetaData.C_PD.csv",sep=",",row.names=1)
rownames(MetaData.C_PD)=gsub("-1",".1",rownames(MetaData.C_PD))
MetaData.C_PD$Cell_IDs=rownames(MetaData.C_PD)
Table=Reduce(merge, lapply(Allfiles, function(x) data.frame(x, rn = row.names(x))))
rownames(Table)=Table$rn
Table$rn=NULL
Table.MX.cpm=cpm(Table)
################################Extract cell-types############################################################
A=read.csv("C_PD_CellTypes.csv",sep=" ",row.names=1)
A2=rownames(A)[which(A$x=="Microglia")]
A2=gsub("-1",".1",A2)
head(A2)
Table.MX.cpm3=Table.MX.cpm[,intersect(A2,colnames(Table.MX.cpm))]

############################################################################################
######################## DataMatrix ##########################################
######### DataMatrix CONTROL AND CASES ####################################################
Ex.MetaData.C_PD=MetaData.C_PD[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender,

fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)

# Make a single cell assay object
sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData.C_PD,fdata)
# count number of detected genes and scale.
Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"PD")
colData(sca)$condition=cond

colData(sca)$Gender=factor(colData(sca)$Gender)
colData(sca)$AgeAtDeath=as.numeric(colData(sca)$AgeAtDeath)
colData(sca)$PMI=as.numeric(colData(sca)$PMI)
colData(sca)$IndID=factor(colData(sca)$IndID)

head(colData(sca))

zlmCond=zlm(~condition+DetRate+Gender+AgeAtDeath+PMI+IndID,sca)

#Run likelihood ratio test for the condition coefficient.

summaryCond = summary(zlmCond, doLRT='conditionControl')

FcThreshold = log2(1.5)

summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionControl' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='conditionControl' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
rownames(fcHurdleSig)=fcHurdleSig$primerid
length(rownames(fcHurdleSig))
fcHurdleSig$primerid=NULL
write.csv(fcHurdleSig,"PD_C.Sig.Genes.CDGAPI.Microglia.csv")
