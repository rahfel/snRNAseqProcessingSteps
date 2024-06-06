#
#set names euqal to the values
sn <- function(x) { names(x) <- x; return(x); }
# filter cells based on the gene/molecule dependency
FilterCells <- function(countMatrix,min.cell.size=500, max.cell.size=5e4,p.level=min(1e-3,1/ncol(countMatrix)),alpha=0.1,do.par=T) {
  if(do.par) { par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);}
  hist(log10(colSums(countMatrix)),col='wheat',xlab='log10[ molecules ]',main='') 
  # some of the cells are very large .. those can skew the analysis of more subtle populations (too much bias) .. letting them in here though
  abline(v=log10(c(min.cell.size,max.cell.size)),lty=2,col=2)
  # look at the number of genes vs. molecule size depenency
  df <- data.frame(molecules=colSums(countMatrix),genes=colSums(countMatrix>0)); 
  df <- df[df$molecules>=min.cell.size,];
  df <- log10(df);
  df <- df[order(df$molecules,decreasing=F),]
  plot(df,col=adjustcolor(1,alpha=alpha),cex=0.5,ylab='log10[ gene counts]',xlab='log10[ molecule counts]')
  abline(v=log10(c(min.cell.size,max.cell.size)),lty=2,col=2)
  #abline(lm(genes ~ molecules, data=df),col=4)
  require(MASS)  
  m <- rlm(genes~molecules,data=df)
  suppressWarnings(pb <- data.frame(predict(m,interval='prediction',level = 1-p.level,type="response")))
  polygon(c(df$molecules,rev(df$molecules)),c(pb$lwr,rev(pb$upr)),col=adjustcolor(2,alpha=0.1),border = NA)
  outliers <- rownames(df)[df$genes > pb$upr | df$genes < pb$lwr];
  points(df[outliers,],col=2,cex=0.6)
  # set of filtered cells to move forward with  
  valid.cells <- colSums(countMatrix)>min.cell.size & colSums(countMatrix)<max.cell.size & !(colnames(countMatrix) %in% outliers)
  countMatrix[,valid.cells,drop=F]
}
# load 10x matrices from a named list of result folders
Load.Data.Matrices <- function(matrixPaths) {
  require(parallel)
  require(Matrix)
  mclapply(sn(names(matrixPaths)),function(nam) {
    matrixPath <- matrixPaths[nam];
    # read all count files (*_unique.counts) under a given path
    #cat("loading data from ",matrixPath, " ");
    x <- as(readMM(paste(matrixPath,'matrix.mtx.gz',sep='/')),'dgCMatrix'); # convert to the required sparse matrix representation
    cat(".")
    gs <- read.delim(paste(matrixPath,'features.tsv.gz',sep='/'),header=F)
    rownames(x) <- gs[,2]
    cat(".")
    gs <- read.delim(paste(matrixPath,'barcodes.tsv.gz',sep='/'),header=F)
    colnames(x) <- gs[,1]
    cat(".")
    colnames(x) <- paste(nam,colnames(x),sep='_');
    x
  },mc.cores=30)
}

FisherTest.wang=function(x)
{
  TMP=matrix(ncol=25,nrow=1)
  Overlap.Genes=intersect(Filtered.Genes,x)
  for(j in 1:25)
  {
    TMP.Genes=intersect(Filtered.Genes,DER.UMI.LIST[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,Overlap.Genes))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,Overlap.Genes))
    TMP.MAT[2,1]=length(setdiff(Overlap.Genes,TMP.Genes))
    TMP.MAT[2,2]=length(Filtered.Genes)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

CellTypeTest.wang=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05)]})
  CellType.Names=c(rep("Unclassified",length(ControlSample1.FisherTest.res.wang)))
  CellType.Names[which(CellType=="Astro")]="Astro"
  CellType.Names[which(CellType=="Endo")]="Endo"
  CellType.Names[which(CellType=="Ex1")]="Excitatory"
  CellType.Names[which(CellType=="Ex2")]="Excitatory"
  CellType.Names[which(CellType=="Ex3e")]="Excitatory"
  CellType.Names[which(CellType=="Ex4")]="Excitatory"
  CellType.Names[which(CellType=="Ex5b")]="Excitatory"
  CellType.Names[which(CellType=="Ex6a")]="Excitatory"
  CellType.Names[which(CellType=="Ex6b")]="Excitatory"
  CellType.Names[which(CellType=="Ex8")]="Excitatory"
  CellType.Names[which(CellType=="Ex9")]="Excitatory"
  CellType.Names[which(CellType=="In1a")]="Inhibitory"
  CellType.Names[which(CellType=="In1b")]="Inhibitory"
  CellType.Names[which(CellType=="In1c")]="Inhibitory"
  CellType.Names[which(CellType=="In3")]="Inhibitory"
  CellType.Names[which(CellType=="In4a")]="Inhibitory"
  CellType.Names[which(CellType=="In4b")]="Inhibitory"
  CellType.Names[which(CellType=="In6a")]="Inhibitory"
  CellType.Names[which(CellType=="In6b")]="Inhibitory"
  CellType.Names[which(CellType=="In7")]="Inhibitory"
  CellType.Names[which(CellType=="In8")]="Inhibitory"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligo")]="Oligo"
  CellType.Names[which(CellType=="OPC")]="OPC"
  CellType.Names[which(CellType=="Per")]="Per"
  CellType.Names
}
 # # # # 

fisher_cell_type.U=function(x)
                {
                        mat.tmp=matrix(ncol=25,nrow=1)
                        genes.tmp=intersect(BG.genes,x)
                        for(j in 1:25)
                        {
                                Oligo.tmp=intersect(BG.genes,DER.UMI.LIST[[j]])
                                mat=matrix(ncol=2,nrow=2)
                                mat[1,1]=length(intersect(Oligo.tmp,genes.tmp))
                                mat[1,2]=length(setdiff(Oligo.tmp,genes.tmp))
                                mat[2,1]=length(setdiff(genes.tmp,Oligo.tmp))
                                mat[2,2]=length(BG.genes)-mat[1,1]-mat[1,2]-mat[2,1]
                                mat.tmp[1,j]=fisher.test(mat,alternative="greater")$p.value
                        }
                        mat.tmp
                }

most.sig.celltype.U=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
  CellType.Names=c(rep("Unclassified",length(FinalBarcodes)))
  CellType.Names[which(CellType=="Astro")]="Astro"
  CellType.Names[which(CellType=="Endo")]="Endo"
  CellType.Names[which(CellType=="Ex1")]="Excitatory"
  CellType.Names[which(CellType=="Ex2")]="Excitatory"
  CellType.Names[which(CellType=="Ex3e")]="Excitatory"
  CellType.Names[which(CellType=="Ex4")]="Excitatory"
  CellType.Names[which(CellType=="Ex5b")]="Excitatory"
  CellType.Names[which(CellType=="Ex6a")]="Excitatory"
  CellType.Names[which(CellType=="Ex6b")]="Excitatory"
  CellType.Names[which(CellType=="Ex8")]="Excitatory"
  CellType.Names[which(CellType=="Ex9")]="Excitatory"
  CellType.Names[which(CellType=="In1a")]="Inhibitory"
  CellType.Names[which(CellType=="In1b")]="Inhibitory"
  CellType.Names[which(CellType=="In1c")]="Inhibitory"
  CellType.Names[which(CellType=="In3")]="Inhibitory"
  CellType.Names[which(CellType=="In4a")]="Inhibitory"
  CellType.Names[which(CellType=="In4b")]="Inhibitory"
  CellType.Names[which(CellType=="In6a")]="Inhibitory"
  CellType.Names[which(CellType=="In6b")]="Inhibitory"
  CellType.Names[which(CellType=="In7")]="Inhibitory"
  CellType.Names[which(CellType=="In8")]="Inhibitory"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligo")]="Oligo"
  CellType.Names[which(CellType=="OPC")]="OPC"
  CellType.Names[which(CellType=="Per")]="Per"
  CellType.Names
}

