# Zebrafish retina GRNs

The codes for constructing GRNs of MG cell groups between LD and NMDA conditions.

### STEP0: Example datasets

We use LD and NMDA scRNAseq/scATACseq datasets as example datasets.

All example files can be downloaded in the following link: [Datasets](https://drive.google.com/drive/folders/1yYuWGWyFog8xhMxbpK26uhdEOh620sz3?usp=sharing)

All functions are in the function folder.

### STEP1: Inferring activators and repressors by expression and motif activity
``` r
### load fish motif name and gene name ###
load("Fish_All_motifs_table_short")
### load chromVAR object generated from chromVAR ###################
load(file='Injury_Motif_seurat_July23')
### 
LD_M_seurat = Injury_Motif_seurat[,which(Injury_Motif_seurat$orig.ident == "LD")]
NMDA_M_seurat = Injury_Motif_seurat[,which(Injury_Motif_seurat$orig.ident == "NMDA")]

### load RNA object ######################
load("Merge_seurat_precursor_final")
LD_seurat = Merge_seurat_precursor_final[,which(Merge_seurat_precursor_final$Condition == "LD")]
NMDA_seurat = Merge_seurat_precursor_final[,which(Merge_seurat_precursor_final$Condition == "NMDA")]

### 
LD_RNA_mat = LD_seurat[['RNA']]@data
NMDA_RNA_mat = NMDA_seurat[['RNA']]@data

######
LD_Motif_mat = LD_M_seurat[['RNA']]@data
NMDA_Motif_mat = NMDA_M_seurat[['RNA']]@data

######
LD_Corr_RNA_Motif_Res = Get_corr_gene_motif(LD_Motif_mat,LD_RNA_mat,Fish_All_motifs_table_short)
NMDA_Corr_RNA_Motif_Res = Get_corr_gene_motif(NMDA_Motif_mat,NMDA_RNA_mat,Fish_All_motifs_table_short)

LD_Corr_RNA_Motif_Res$TAG = 'LD'
save(LD_Corr_RNA_Motif_Res,file='LD_Corr_RNA_Motif_Res_July23')
NMDA_Corr_RNA_Motif_Res$TAG = 'NMDA'
save(NMDA_Corr_RNA_Motif_Res,file='NMDA_Corr_RNA_Motif_Res_July23')

colnames(LD_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(LD_Corr_RNA_Motif_Res)[c(3)],sep='_') 
colnames(NMDA_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(NMDA_Corr_RNA_Motif_Res)[c(3)],sep='_') 
LD_NMDA_TF_activity <- rbind(LD_Corr_RNA_Motif_Res,NMDA_Corr_RNA_Motif_Res)
library(writexl)
write_xlsx(LD_NMDA_TF_activity, path = "LD_NMDA_TFs_gene_motif_corr.xlsx")
```


### STEP2: Identifying Cis-regulatory elements
``` r

```


### STEP3: Predicting TF Binding Sites
``` r

```


### STEP4: TF-target correlation 
``` r

```


### STEP5: Construction of TF-peak-target links
``` r

```


### STEP6: Identification of enriched gene regulatory sub-networks
``` r

```


### STEP7: Identication of key activator TFs
``` r

```





