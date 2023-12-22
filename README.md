# Zebrafish retina GRNs

The test codes for constructing GRNs of MG cell groups between LD and NMDA conditions.

### Datasets ###
The fragment files of snATACseq can be downloaded in the following link:[Datasets](https://drive.google.com/drive/folders/1N646KmKq8NAB-_nISVeIqu30YcOannvy?usp=sharing)

The annotated snRNAseq single cell annotation and UMAPs in Figure5A can be downloaded in the following link:[Datasets](https://drive.google.com/drive/folders/1N646KmKq8NAB-_nISVeIqu30YcOannvy?usp=sharing).The cell type annotations are in the "celltype" column.
Please load the object with "readRDS" command.

### System requirements ###
1.The code can be run on any operating system that has R installed.
2.ArchR (v1.0.2) and Seurat (v4.0.6) package are required to test the code.
3.No non-standard hardware is required.

### Installation guide ###
1.No need to install in an R environment.

### Expected run time ###
1.The expected construction time for GRNs is approximately 2.5h.

### STEP0: Example datasets

We use LD and NMDA scRNAseq/scATACseq datasets as example datasets.

All example files can be downloaded in the following link: [Datasets](https://drive.google.com/drive/folders/1yYuWGWyFog8xhMxbpK26uhdEOh620sz3?usp=sharing)

### STEP1: Inferring activators and repressors by expression and motif activity
``` r
### load fish motif name and gene name ###
load("Fish_All_motifs_table_short")
### load chromVAR object generated from chromVAR ###################
load(file='Injury_Motif_seurat_July23')

LD_M_seurat = Injury_Motif_seurat[,which(Injury_Motif_seurat$orig.ident == "LD")]
NMDA_M_seurat = Injury_Motif_seurat[,which(Injury_Motif_seurat$orig.ident == "NMDA")]

### load RNA object ######################
load("Merge_seurat_precursor_final")
LD_seurat = Merge_seurat_precursor_final[,which(Merge_seurat_precursor_final$Condition == "LD")]
NMDA_seurat = Merge_seurat_precursor_final[,which(Merge_seurat_precursor_final$Condition == "NMDA")]

LD_RNA_mat = LD_seurat[['RNA']]@data
NMDA_RNA_mat = NMDA_seurat[['RNA']]@data

LD_Motif_mat = LD_M_seurat[['RNA']]@data
NMDA_Motif_mat = NMDA_M_seurat[['RNA']]@data

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

load("LD_Corr_RNA_Motif_Res_July23")
load("NMDA_Corr_RNA_Motif_Res_July23")

colnames(LD_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(LD_Corr_RNA_Motif_Res)[c(3)],sep='_') 
colnames(NMDA_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(NMDA_Corr_RNA_Motif_Res)[c(3)],sep='_') 

LD_NMDA_TF_activity <- rbind(LD_Corr_RNA_Motif_Res,NMDA_Corr_RNA_Motif_Res)
save(LD_NMDA_TF_activity,file='LD_NMDA_TF_activity_July23')

```


### STEP2: Identifying Cis-regulatory elements
``` r
### call peaks using ArchR then calcualte PtoG correlations ###
### LD_project is the ArchR object ####
### NMDA_project is the ArchR object ##

load("All_Peaks_GRNs")

LD_PtoG = Get_PtoG_fun(LD_project,Peakmat_all=GRanges(All_Peaks_GRNs))
NMDA_PtoG = Get_PtoG_fun(NMDA_project,Peakmat_all=GRanges(All_Peaks_GRNs))

save(LD_PtoG,file='LD_PtoG_July23')
save(NMDA_PtoG,file='NMDA_PtoG_July23')

```



``` r
###
### Next merge PtoG with Gene annotations ####
###

### Project_merge_sub_pre_cl is the ArchR object ####

load("All_Peaks_GRNs")

GTF_TSS_table_res = Get_TSS_table(Gtfs_protein_trans_GR_TSS,All_Peaks_GRNs,extend=500)
GTF_Body_table_res = Get_Body_table(Gtfs_protein_trans_GR_Body,All_Peaks_GRNs)

load(file='LD_PtoG_July23')
load(file='NMDA_PtoG_July23')

LD_PtoG_Ann = Get_all_peak_Gene_tables(GTF_TSS_table_res,GTF_Body_table_res,PtoG=LD_PtoG,Gtfs_protein_trans_GR_TSS)
NMDA_PtoG_Ann = Get_all_peak_Gene_tables(GTF_TSS_table_res,GTF_Body_table_res,PtoG=NMDA_PtoG,Gtfs_protein_trans_GR_TSS)

save(LD_PtoG_Ann,file='LD_PtoG_Ann_July23')
save(NMDA_PtoG_Ann,file='NMDA_PtoG_Ann_July23')
```



### STEP3: Predicting TF Binding Sites

#### The pipline of TOBIAS please see: https://github.com/Pinlyu3/IReNA-v2 
``` r
#### select the activator and repressor motifs which passed the cutoff ####
#### cutoff = 0.05 #####
setwd("/zp1/data/plyu3/Fish_muti_omic_dev/merge_ATAC")

load("LD_NMDA_TF_activity_July23")
load("InjuryDev_TF_activity_July23")

Globle_need_Motif = Get_all_need_Motif_tag(LD_NMDA_TF_activity,InjuryDev_TF_activity,0.05)

####
load(file='Dev_PtoG_Ann_July23')
load(file='Injury_PtoG_Ann_July23')
load(file='LD_PtoG_Ann_July23')
load(file='NMDA_PtoG_Ann_July23')

x1 = Dev_PtoG_Ann
x2 = Injury_PtoG_Ann
x3 = LD_PtoG_Ann
x4 = NMDA_PtoG_Ann

######## merge all peak regions to match TF binding mtoifs ##########
######## Globle_need_Peak is the union regions of all these peaks ####
Globle_need_Peak = Get_all_need_Peak_tag(x1,x2,x3,x4)

#### Fish_combined_motifs is the PWMList which stored Fish PMW matrix #####
load(file='Fish_combined_motifs') ### Fish_combined ####

#### Fish_combined
#### PWMatrixList of length 9155
#### names(9155): M00001 M00002 M00004 ... M11389_2.00 M11390_2.00 M11491_2.00 

Globle_motif_footprint <- Match_motifs(Globle_need_Peak,Globle_need_Motif,Fish_combined)
save(Globle_motif_footprint,file='Globle_motif_footprint_July23')

#### Next using LD or NMDA footprint to filter motif binding sites ######
#### Next load the corrected ATACseq signal calculated by TOBIAS ########

#### The pipline of TOBIAS please see: https://github.com/Pinlyu3/IReNA-v2 
#### STEP4: Predicting cell-type specific TFs binding in cis-regulatory elements
 
#### loading corrected signal by rtracklayer #####
NMDA_injury_footprint <- rtracklayer::import.bw("NMDA_injury_MG_fragments_GR_bam_tab_s_corrected.bw")
LD_injury_footprint <- rtracklayer::import.bw("LD_injury_MG_fragments_GR_bam_tab_s_corrected.bw")

Globle_motif_cleanedByLD = Get_footprint_Score_all_para(Globle_motif_footprint,footprint_signal=LD_injury_footprint)
save(Globle_motif_cleanedByLD,file='Globle_motif_cleanedByLD')

Globle_motif_cleanedByNMDA = Get_footprint_Score_all_para(Globle_motif_footprint,footprint_signal=NMDA_injury_footprint)
save(Globle_motif_cleanedByNMDA,file='Globle_motif_cleanedByNMDA')

```


### STEP4: TF-target correlation 

``` r
### calculate the gene-gene correlation ###

library(reshape2)
LD_mat = read.table("injury.LD.mat.txt",sep='\t',header=T)
LD_mat = as.matrix(LD_mat)
k = which(colSums(LD_mat) == 0)
LD_mat_cl = LD_mat[,-k]
LD_Corr_res = sparse.cor3(LD_mat_cl)
LD_Corr_res = melt(LD_Corr_res)
save(LD_Corr_res,file="LD_Corr_res")

#### injury_LD_network_new.tsv is output from grnboost2 in arboreto package #####
LD_network <- read.table("injury_LD_network_new.tsv",sep="\t",header=F)
colnames(LD_network) <- c("TF","Gene","Score")
LD_network$index = paste(LD_network$TF,LD_network$Gene,sep="-->")
k1 = which(LD_Corr_res$Var1 %in% LD_network$TF == T & LD_Corr_res$Var2 %in% LD_network$Gene == T)
LD_Corr_res_cl = LD_Corr_res[k1,]
LD_Corr_res_index = paste(LD_Corr_res_cl$Var1,LD_Corr_res_cl$Var2,sep="-->")
m = match(LD_network$index,LD_Corr_res_index)
LD_network$Corr = LD_Corr_res_cl$value[m]

save(LD_network,file="LD_network")
```



### STEP5: Construction of TF-peak-target links
``` r
#### load footprint, motif activity, peak-gene, and gene-gene correlations

load("Globle_motif_cleanedByNMDA") ### from Step3 ####
load("LD_NMDA_TF_activity_July23") ### from Step1 ####
load("NMDA_PtoG_Ann_July23") ### from Step2 ####
load("NMDA_network") ### from Step4 ####

GG_network = NMDA_network
PG_cor = NMDA_PtoG_Ann
TP_cor = LD_NMDA_TF_activity[which(LD_NMDA_TF_activity$TAG=='NMDA'),]
Foot = Globle_motif_cleanedByNMDA

NMDA_GRNs = Combined_all_things(GG_network,PG_cor,TP_cor,Foot)
save(NMDA_GRNs,file='NMDA_GRNs_July23')

load("Globle_motif_cleanedByLD") ### from Step3 ####
load("LD_NMDA_TF_activity_July23") ### from Step1 ####
load("LD_PtoG_Ann_July23") ### from Step2 ####
load("LD_network") ### from Step4 ####

GG_network = LD_network
PG_cor = LD_PtoG_Ann
TP_cor = LD_NMDA_TF_activity[which(LD_NMDA_TF_activity$TAG=='LD'),]
Foot = Globle_motif_cleanedByLD

LD_GRNs = Combined_all_things(GG_network,PG_cor,TP_cor,Foot)
save(LD_GRNs,file='LD_GRNs_July23')
```


### STEP6: Identification of enriched gene regulatory sub-networks
``` r
#### All_GRNs_need_Genes_July23 is the genes which expressed in MG group cells #####
load("All_GRNs_need_Genes_July23")
Genes_need = All_GRNs_need_Genes

#### load the LD_GRNs and NMDA_GRNs (total GRNs)
load("LD_GRNs_July23")
load("NMDA_GRNs_July23")

#### Next convert to triple pairs ######
#### only kept the genes expressed in MG group cells ####

LD_GRNs_MG_triple = Convert_to_triple_and_filterbygenes(LD_GRNs,Genes_need)
NMDA_GRNs_MG_triple = Convert_to_triple_and_filterbygenes(NMDA_GRNs,Genes_need)

#### Next convert to double pairs ######
#### only kept the genes expressed in MG group cells ####

LD_GRNs_MG_double = Convert_to_double_and_filterbygenes(LD_GRNs,Genes_need)
NMDA_GRNs_MG_double = Convert_to_double_and_filterbygenes(NMDA_GRNs,Genes_need)

##### Next load the log2 fold change files between LD and NMDA #######

##### Differential DARs between LD and NMDA
load("DARs_Combined_LDvsNMDA_July23")

##### Differential DEGs between LD and NMDA
load('DEGs_Combined_LDNMDA_July23')

##### call specific GRNs by the fold change of Genes and Peaks ######

LD_GRNs_MG_triple_specific <- Filter_triple_network_according_to_foldchange(LD_GRNs_MG_triple,DEGs,DARs,tag='pos')
NMDA_GRNs_MG_triple_specific <- Filter_triple_network_according_to_foldchange(NMDA_GRNs_MG_triple,DEGs,DARs,tag='neg')

##### 
library(writexl)
write_xlsx(LD_GRNs_MG_triple_specific, path = "LD_GRNs_MG_triple_specific.xlsx")
write_xlsx(NMDA_GRNs_MG_triple_specific, path = "NMDA_GRNs_MG_triple_specific.xlsx")
```





