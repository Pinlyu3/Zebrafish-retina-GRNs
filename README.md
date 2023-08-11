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
``` python
import os
import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

import pandas as pd
import arboreto
import dask
from dask.distributed import Client
from distributed import LocalCluster, Client

ex_matrix = pd.read_csv("injury.LD.mat.txt", sep='\t')
col_sums = ex_matrix.sum(axis=0)
ex_matrix = ex_matrix.loc[:, col_sums != 0]
tf_names = load_tf_names("GRNs_TF.txt")

local_cluster = LocalCluster(n_workers=30,threads_per_worker=1)
custom_client = Client(local_cluster)

network = grnboost2(expression_data=ex_matrix,tf_names=tf_names,client_or_address=custom_client)
network.head()
network.to_csv('injury_LD_network_new.tsv', sep='\t', header=False, index=False)

```

``` r
### calculate the gene-gene correlation ###

sparse.cor3 <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    cSums <- colSums(x)
    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp
    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}


library(reshape2)
LD_mat = read.table("injury.LD.mat.txt",sep='\t',header=T)
LD_mat = as.matrix(LD_mat)
k = which(colSums(LD_mat) == 0)
LD_mat_cl = LD_mat[,-k]
LD_Corr_res = sparse.cor3(LD_mat_cl)
LD_Corr_res = melt(LD_Corr_res)
save(LD_Corr_res,file="LD_Corr_res")

LD_network <- read.table("injury_LD_network_new.tsv",sep="\t",header=F)
colnames(LD_network) <- c("TF","Gene","Score")
LD_network$index = paste(LD_network$TF,LD_network$Gene,sep="-->")
k1 = which(LD_Corr_res$Var1 %in% LD_network$TF == T & LD_Corr_res$Var2 %in% LD_network$Gene == T)
LD_Corr_res_cl = LD_Corr_res[k1,]
LD_Corr_res_index = paste(LD_Corr_res_cl$Var1,LD_Corr_res_cl$Var2,sep="-->")
m = match(LD_network$index,LD_Corr_res_index)
LD_network$Corr = LD_Corr_res_cl$value[m]

save(LD_network,file="LD_network")

#### save(NDMA_network,file="NMDA_network") in the google folder ####
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





