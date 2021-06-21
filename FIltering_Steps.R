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

source("UseFulRfunction.R")
DataMatrix = Load.Data.Matrices(list(C36_Control='/rds/general/user/rf1116/ephemeral/PD_Samples/C36_raw_feature_bc_matrix'))

DataMatrix.B=DataMatrix[[1]]
set.seed(1234)
E.out = emptyDrops(DataMatrix.B,lower=300,retain=300)
E.keep = E.out$FDR <= 0.001
E.keep.na = E.keep
E.keep.na[is.na(E.keep)]= FALSE
summary(E.keep.na)


ControlSample1.data= CreateSeuratObject(counts =as.matrix(DataMatrix.B[,E.keep.na]), min.cells = 5,project = "pd")

ControlSample1.data
ControlSample1.data[["MT.PER"]] = PercentageFeatureSet(ControlSample1.data, pattern = "^MT-")
ControlSample1.data.A = subset(ControlSample1.data, subset = MT.PER < 5)
ControlSample1.data.A

# check how many genes have at least one transcript in each cell
ControlSample1.data.Matrix=GetAssayData(object = ControlSample1.data.A, assay = "RNA", slot = "data")

GenesDetected = apply(ControlSample1.data.Matrix, 1, function(x) sum(x>0))
table(GenesDetected>=5) # a lot of genes are not detected in 5 or more cells
keep = GenesDetected>=5
ControlSample1.data.Matrix.B= ControlSample1.data.Matrix[keep,]
dim(ControlSample1.data.Matrix.B)

ControlSample1.data.B= CreateSeuratObject(counts =as.matrix(ControlSample1.data.Matrix.B))
ControlSample1.data.B[["MT.PER"]] = PercentageFeatureSet(ControlSample1.data.B, pattern = "^MT-")

#ControlSample1.data.B.RawCountsMatrix=GetAssayData(object = ControlSample1.data.B, assay = "RNA", slot = "data")

ControlSample1.data.B= NormalizeData(ControlSample1.data.B, normalization.method = "LogNormalize", scale.factor = 10000)

## Find variable features,runPCA etc
ControlSample1.data.B = FindVariableFeatures(ControlSample1.data.B, selection.method = "mean.var.plot", nfeatures = 2000)
All.genes = rownames(ControlSample1.data.B)
ControlSample1.data.B = ScaleData(ControlSample1.data.B, features = All.genes)
ControlSample1.data.B = RunPCA(ControlSample1.data.B,npcs = 30, verbose = FALSE)
ControlSample1.data.B = FindNeighbors(ControlSample1.data.B, dims = 1:10)
ControlSample1.data.B= FindClusters(ControlSample1.data.B, resolution = 2)
ControlSample1.data.B = RunUMAP(ControlSample1.data.B, dims = 1:10)
pdf("C36.OE.UMAP2.pdf")
DimPlot(ControlSample1.data.B, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()
ControlSample1.data.FinalMatrix=GetAssayData(object = ControlSample1.data.B, assay = "RNA", slot = "data")


Con.res.list = paramSweep_v3(ControlSample1.data.B, PCs = 1:10)
Con.res.list.stats = summarizeSweep(Con.res.list, GT = FALSE)
FindPK = find.pK(Con.res.list.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations =ControlSample1.data.B@meta.data$RNA_snn_res.0.5
Prop.Homotypic = modelHomotypic(annotations)
nExp_poi = round(0.07*length(ControlSample1.data.B@active.ident))
nExp_poi.adj = round(nExp_poi*(1-Prop.Homotypic))

## Run DoubletFinder with varying classification stringencies -

ControlSample1.data.B = doubletFinder_v3(ControlSample1.data.B, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE,sct = FALSE)
ControlSample1.data.B = doubletFinder_v3(ControlSample1.data.B, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25_0.09_",nExp_poi,sep=""), sct = FALSE)

ControlSample1.data.B@meta.data[,"DF_Class"] = get('ControlSample1.data.B')@meta.data[[paste("DF.classifications_0.25_0.09_",nExp_poi,sep="")]]
ControlSample1.data.B@meta.data$DF_Class[which(ControlSample1.data.B@meta.data$DF_Class == "Doublet" & get('ControlSample1.data.B')@meta.data[[paste("DF.classifications_0.25_0.09_",nExp_poi,sep="")]] == "Singlet")] = "Doublet_LOW"
ControlSample1.data.B@meta.data$DF_Class[which(ControlSample1.data.B@meta.data$DF_Class == "Doublet")] = "Doublet_HIGH"
A.Doublets=as.character(colnames(ControlSample1.data.B)[which(ControlSample1.data.B@meta.data$DF_Class == "Doublet_HIGH")])

ControlSample1.data.C=ControlSample1.data.B[, !(colnames(ControlSample1.data.B) %in% A.Doublets), drop = FALSE]


ControlSample1.data.A.markers= FindAllMarkers(ControlSample1.data.C, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ControlSample1.data.A.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
ControlSample1.data.A.markers.LIST=unstack(ControlSample1.data.A.markers, ControlSample1.data.A.markers$gene ~ ControlSample1.data.A.markers$cluster)
saveRDS(ControlSample1.data.A.markers.LIST,"C36.OE.markers.rds")

DER.UMI=read.csv"DER-21_Single_cell_markergenes_UMI.csv",sep=",",check.names=FALSE)
DER.UMI.LIST=unstack(DER.UMI, (DER.UMI$GeneName) ~ DER.UMI$CellType)

# Create an object for wang annotation
ControlSample1.data.A.wang=ControlSample1.data.C
Filtered.Genes=unique(unlist(ControlSample1.data.A.markers.LIST))
ControlSample1.FisherTest.res.wang=lapply(ControlSample1.data.A.markers.LIST,FisherTest.wang)
TMP.wang = matrix(unlist(ControlSample1.FisherTest.res.wang), ncol = 25, byrow = TRUE)
colnames(TMP.wang)=names(DER.UMI.LIST)
write.csv(TMP.wang,"C36.OE.TMP.wang.csv")
TMP.LABELS.wang=CellTypeTest.wang(TMP.wang)
names(TMP.LABELS.wang)=names(ControlSample1.FisherTest.res.wang)
wang.cluster.combined.ids = TMP.LABELS.wang
names(wang.cluster.combined.ids) = levels(ControlSample1.data.A.wang)
ControlSample1.data.A.wang = RenameIdents(ControlSample1.data.A.wang,wang.cluster.combined.ids)
levels(ControlSample1.data.A.wang)
pdf("C36.OE.LABELLED.pdf")
DimPlot(ControlSample1.data.A.wang, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()
table(Idents(ControlSample1.data.A.wang))
levels(ControlSample1.data.A.wang)

ControlSample1.data.A.wang[["ClusterIdent"]] = Idents(object = ControlSample1.data.A.wang)
ClusterIdent=as.character(ControlSample1.data.A.wang@meta.data$ClusterIdent)
names(ClusterIdent)=rownames(ControlSample1.data.A.wang@meta.data)

saveRDS(as.factor(ClusterIdent),"C36.OE.CellID_Type.rds")

#Extract counts only for filtered matrix reads
ControlSample1.data.FinalMatrix=GetAssayData(object = ControlSample1.data.A.wang, assay = "RNA", slot = "data")
ControlSample1.data.FinalMatrix2=ControlSample1.data.B.RawCountsMatrix[,colnames(ControlSample1.data.FinalMatrix)]
dim(ControlSample1.data.FinalMatrix2)
dim(ControlSample1.data.FinalMatrix)
saveRDS(ControlSample1.data.FinalMatrix2,"C36.OE.NonNormalised.filteredcounts.dgCMatrix.rds")
head(colnames(ControlSample1.data.FinalMatrix2))
