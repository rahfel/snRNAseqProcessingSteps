

library(WebGestaltR)

tmp = list.files(pattern="*.PD_C.txt")
Allfiles = lapply(tmp,read.delim,sep="\t")
filenames=list.files(pattern="*.PD_C.txt", full.names=TRUE)
filenames
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)

Inpath="C:/Users/rahel/Documents/PD_Webgestalt/PD_C/"
DataBase="C:/Users/rahel/Documents/PD_Webgestalt/REACTOME_KEGG_PID_PATHWAYS.gmt"
REF="C:/Users/rahel/Documents/PD_Webgestalt/PD_C/PD_C_REF_GENES.txt"

for (i in 1:length(Allfiles)){
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabaseFile=DataBase,enrichDatabaseType="genesymbol",interestGeneFile=paste0(Inpath,names(Allfiles)[i],".txt"),
 interestGeneType="genesymbol", referenceGeneFile=REF,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory=Inpath, projectName=paste0("REACTOME_KEGG_PID_PATHWAYS.",names(Allfiles)[i]))}
