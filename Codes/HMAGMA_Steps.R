
PATH=$PATH:/Users/rf1116/Downloads/magma_v1
Gene_annot_out=/Users/rf1116/Documents/HMAGMA/ADULT_BRAIN.magma.genes.raw
DEGenes_File=/Users/rf1116/Documents/HMAGMA/PD_C.AllRegulated_ENTID.txt
magma --gene-results $Gene_annot_out --set-annot $DEGenes_File 'col=2,1' --out PD2019NALLS_HMAGMA_ADULT_BRAIN_PD_C


C=read.table("PD2019NALLS_HMAGMA_ADULT_BRAIN_PD_C.gsa.out",sep="",header=T)
C$FDR=p.adjust(C$P, method = "BH", n = length(C$P))

rownames(C)=C$VARIABLE
write.table(C,"/Users/rf1116/Documents/PD10X/MetaDatasets/umap_DE/PD_January/PD_FEB/MAGMA_Results/PD2019NALLS_HMAGMA_ADULT_BRAIN_PD_C.gsa.out.csv")
