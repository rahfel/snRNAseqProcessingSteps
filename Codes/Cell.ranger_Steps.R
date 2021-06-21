
# STEP 1: MAKE GTF FILE

mkdir refdata-cellranger-GRCh38-1.2.0_premrna
cd refdata-cellranger-GRCh38-1.2.0_premrna

awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' \
       /rds/general/user/rf1116/ephemeral/Parkinsons/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf > GRCh38-1.2.0.premrna.gtf
       # /rds/general/user/ckhozoie/ephemeral/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf > GRCh38-1.2.0.premrna.gtf


export PATH=/rds/general/user/rf1116/ephemeral/Parkinsons/cellranger-3.0.2:$PATH

cellranger mkref --genome=GRCh38-1.2.0_premrna \
                   --fasta=/rds/general/user/rf1116/ephemeral/Parkinsons/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa \
                   --genes=GRCh38-1.2.0.premrna.gtf \
                   --nthreads=32


# STEP2: RUN Cellranger

export PATH=/rds/general/user/rf1116/home/SCRATCH/cellranger-3.0.2:$PATH


#ID sample ID
#NAME Folder name
FQ1=/rds/general/user/rf1116/ephemeral/MS/MS/BATCH2/H2JHHBBXY/IGFQ000667_matthews_5-11-2018_SC/$ID/
cellranger count  --id=$NAME --transcriptome=/rds/general/ephemeral/user/rf1116/ephemeral/Parkinsons/GRCh38-1.2.0_premrnaA --fastqs=$FQ1 --sample=$NAME --expect-cells=10000 --jobmode=pbspro --maxjobs=5 --localcores=10
