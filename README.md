# Zebrafish retina GRNs

The test codes for constructing GRNs of MG cell groups between LD and NMDA conditions.

### System requirements ###
1.The code can be run on any operating system that has R installed.
2.ArchR (v1.0.2) and Seurat (v4.0.6) package are required to test the code.
3.No non-standard hardware is required.

### Installation guide ###
1.No need to install in an R environment.

### Expected run time ###
1.The expected construction time for GRNs is approximately 30 minutes.

### STEP0: Example datasets

We use LD and NMDA scRNAseq/scATACseq datasets as example datasets.

All example files can be downloaded in the Zenodo (https://doi.org/10.5281/zenodo.8317611)
or in the following link: [Datasets](https://drive.google.com/drive/folders/1yYuWGWyFog8xhMxbpK26uhdEOh620sz3?usp=sharing)

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

Get_corr_gene_motif <- function(motif_matrix,gene_matrix,Fish_All_motifs_table_short){
    ########
    colnames(motif_matrix) <- gsub("#",'_',colnames(motif_matrix))
    rownames(motif_matrix) <- gsub("-",'_',rownames(motif_matrix))
    ########
    print(head(colnames(motif_matrix)))
    print(head(rownames(motif_matrix)))
    print(head(colnames(gene_matrix)))
    print(head(rownames(gene_matrix)))
    ########
    k = which(Matrix::rowSums(gene_matrix) < 10)
    if(length(k) > 1){
        gene_matrix = gene_matrix[-k,]
    }
    rownames(gene_matrix) <- sapply(strsplit(rownames(gene_matrix),split='~~'),function(x) x[[2]])
    ######## clean cols #####
    cells = colnames(gene_matrix)[which(colnames(gene_matrix) %in% colnames(motif_matrix) == T)]
    gene_matrix_cl = gene_matrix[,which(colnames(gene_matrix) %in% cells == T)]
    motif_matrix_cl = motif_matrix[,which(colnames(motif_matrix) %in% cells == T)]
    #######
    gM = match(cells,colnames(gene_matrix_cl))
    gene_matrix_cl = gene_matrix_cl[,gM]
    mM = match(cells,colnames(motif_matrix_cl))
    motif_matrix_cl = motif_matrix_cl[,mM]
    ########
    all.equal(colnames(gene_matrix_cl),colnames(motif_matrix_cl))
    ########
    Fish_All_motifs_table_short$Cor = 0
    ########
    Fish_All_motifs_table_short_cl = Fish_All_motifs_table_short[which(Fish_All_motifs_table_short$Gene %in% rownames(gene_matrix_cl) == T),]
    Fish_All_motifs_table_short_cl$Cor = 0
    ########
    ########
    library(parallel)
    cl <- makeCluster(30)
    cor_res = parApply(cl,Fish_All_motifs_table_short_cl,1,cor_functon,motif_matrix_cl,gene_matrix_cl)
    stopCluster(cl)
    #########
    Fish_All_motifs_table_short_cl$Cor = cor_res
    return(Fish_All_motifs_table_short_cl)
}


cor_functon <- function(x,motif_matrix_cl,gene_matrix_cl){
    tmp1 = x[1]
    tmp2 = x[2]
    #######
    vector1 = which(rownames(motif_matrix_cl) == tmp1)
    vector2 = which(rownames(gene_matrix_cl) == tmp2)
    if(length(vector1) == 1 & length(vector2) == 1){
            v1 = motif_matrix_cl[vector1,]
            v2 = gene_matrix_cl[vector2,]
            ### all.equal(names(v1),names(v2))
            Cor = cor(v1,v2,method = c("spearman"))
            ###
            return(Cor)
    }else{
        Cor = 0
        return(Cor)
    }
}
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


#######
load("LD_Corr_RNA_Motif_Res_July23")
load("NMDA_Corr_RNA_Motif_Res_July23")

colnames(LD_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(LD_Corr_RNA_Motif_Res)[c(3)],sep='_') 
colnames(NMDA_Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(NMDA_Corr_RNA_Motif_Res)[c(3)],sep='_') 


LD_NMDA_TF_activity <- rbind(LD_Corr_RNA_Motif_Res,NMDA_Corr_RNA_Motif_Res)
save(LD_NMDA_TF_activity,file='LD_NMDA_TF_activity_July23')

####
```


### STEP2: Identifying Cis-regulatory elements
``` r
### call peaks using ArchR then calcualte PtoG correlations ###
### LD_project is the ArchR object ####
### NMDA_project is the ArchR object ##

Get_p2g_fun <- function(x){
    corCutOff = 0.25
    FDRCutOff = 1e-6
    varCutOffATAC = 0.7
    varCutOffRNA = 0.3
    p2g <- metadata(x@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
    if(!is.null(varCutOffATAC)){
        p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    if(!is.null(varCutOffRNA)){
        p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }
    mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
    mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
    p2g$peak <- paste0(rowRanges(mATAC))
    p2g$gene <- rowData(mRNA)$name
    return(p2g)
}

Get_PtoG_fun <- function(All_project_cl,Peakmat_all){
    All_project_cl <- addPeakSet(
        ArchRProj = All_project_cl,
        peakSet = Peakmat_all,
        force = TRUE
    )
    #####
    All_project_cl <- addPeakMatrix(
        ArchRProj = All_project_cl,
        ceiling = 4,
        binarize = FALSE,
        verbose = TRUE,
        threads = getArchRThreads(),
        force = TRUE
    )
    #####
    All_project_cl <- addPeak2GeneLinks(
        ArchRProj = All_project_cl,
        useMatrix = "GeneExpressionMatrix",
        reducedDims = "LSI_Combined"
    )
    ####
    PtoG <- Get_p2g_fun(All_project_cl)
    ####
    return(PtoG)
}


LD_PtoG = Get_PtoG_fun(LD_project,Peakmat_all=GRanges(All_Peaks_GRNs))
NMDA_PtoG = Get_PtoG_fun(NMDA_project,Peakmat_all=GRanges(All_Peaks_GRNs))

save(LD_PtoG,file='LD_PtoG_July23')
save(NMDA_PtoG,file='NMDA_PtoG_July23')

### LD_PtoG_July23 and NMDA_PtoG_July23 are in the shared links ####

```



``` r
###
### Next merge PtoG with Gene annotations ####
###

### Project_merge_sub_pre_cl is the ArchR object ####

load("All_Peaks_GRNs")"


###
Get_TSS_table <- function(GR_list,GR_peak,extend=2000){
    ##########
    outtable_list <- list()
    for(i in 1:length(GR_list)){
        tmp_GR = GR_list[[i]]
        ######
        start(tmp_GR) = start(tmp_GR) - extend
        end(tmp_GR) = end(tmp_GR) + extend
        ######
        k = which(countOverlaps(GR_peak,tmp_GR) > 0)
        if(length(k) > 0){
            GR_peak_sub = GR_peak[k]
            #####
            outtable_tmp = data.frame(Gene=tmp_GR$genes[1],Peak=as.character(GR_peak_sub),Class="TSS")
            outtable_list = c(outtable_list,list(outtable_tmp))
        }
        if(i %% 1000 == 0){
            print(i)
        }
    }
    outtable_list_merge = do.call("rbind",outtable_list)
    ##########
    return(outtable_list_merge)
}

Get_Body_table <- function(GR_list,GR_peak){
    ##########
    outtable_list <- list()
    for(i in 1:length(GR_list)){
        tmp_GR = GR_list[[i]]
        ######
        ######
        k = which(countOverlaps(GR_peak,tmp_GR) > 0)
        if(length(k) > 0){
            GR_peak_sub = GR_peak[k]
            #####
            outtable_tmp = data.frame(Gene=tmp_GR$genes[1],Peak=as.character(GR_peak_sub),Class="Body")
            outtable_list = c(outtable_list,list(outtable_tmp))
        }
        if(i %% 1000 == 0){
            print(i)
        }
    }
    outtable_list_merge = do.call("rbind",outtable_list)
    ##########
    return(outtable_list_merge)
}

GTF_TSS_table_res = Get_TSS_table(Gtfs_protein_trans_GR_TSS,All_Peaks_GRNs,extend=500)
GTF_Body_table_res = Get_Body_table(Gtfs_protein_trans_GR_Body,All_Peaks_GRNs)


Get_all_peak_Gene_tables <- function(GTF_TSS_table_res,GTF_Body_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
    ####
    GTF_TSS_table_res_1 = Process_TSS_peaks(GTF_TSS_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
    GTF_Body_table_res_1 = Process_Body_peaks(GTF_Body_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
    ####
    GTF_TSSBody_table_res = Merge_TSS_and_Body(GTF_TSS_table_res_1,GTF_Body_table_res_1)
    ####
    GTF_Inter_table_res_1 = Process_Inter_peaks(GTF_TSSBody_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
    ####
    ALL_peak_table = rbind(GTF_TSS_table_res_1,GTF_Body_table_res_1,GTF_Inter_table_res_1)
    ####
    ALL_peak_table_anno = Process_peak_table_to_annotation(ALL_peak_table)
    ####
    return(ALL_peak_table_anno)
}

Process_Body_peaks <- function(GTF_Body_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
        ###### First add PtoG link correaltions !!! ###########
        GTF_Body_table_res$Index = paste0(GTF_Body_table_res$Gene,"#",GTF_Body_table_res$Peak)
        ######
        GTF_Body_table_res$PtoG_cor = "ND"
        GTF_Body_table_res$PtoG_FDR = "ND"
        ######
        PtoG_index = paste0(PtoG$gene,"#",PtoG$peak)
        ######
        m = match(GTF_Body_table_res$Index,PtoG_index)
        GTF_Body_table_res$PtoG_cor = PtoG$Correlation[m]
        GTF_Body_table_res$PtoG_FDR = PtoG$FDR[m]
        #######
        k = which(is.na(GTF_Body_table_res$PtoG_FDR) == T)
        GTF_Body_table_res_cl = GTF_Body_table_res[-k,]
        #######
        ####### Next see the distance to the nearnest TSS gene ########
        #######
        GTF_Body_table_res_cl_GR = GRanges(GTF_Body_table_res_cl$Peak)
        #######
        dis_all = c()
        #######
        for(i in 1:length(GTF_Body_table_res_cl_GR)){
                tmp = GTF_Body_table_res_cl_GR[i]
                tmp_gene = GTF_Body_table_res_cl$Gene[i]
                ######
                k = which(names(Gtfs_protein_trans_GR_TSS) == tmp_gene)
                tmp_gene_tss = Gtfs_protein_trans_GR_TSS[[k]]
                ######
                tmp_dis = distance(tmp,tmp_gene_tss)
                out = min(na.omit(tmp_dis))
                dis_all = c(dis_all,out)
                if(i %% 1000 == 0){
                        print(i)
                }
        }
        #######
        #######
        GTF_Body_table_res_cl$DistanceTSS = dis_all
        #######
        return(GTF_Body_table_res_cl)
}

Process_TSS_peaks <- function(GTF_TSS_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
        GTF_TSS_table_res$Index = paste0(GTF_TSS_table_res$Gene,"#",GTF_TSS_table_res$Peak)
        ###### First add PtoG link correaltions !!! ###########
        ######
        GTF_TSS_table_res$PtoG_cor = "ND"
        GTF_TSS_table_res$PtoG_FDR = "ND"
        ######
        PtoG_index = paste0(PtoG$gene,"#",PtoG$peak)
        ######
        m = match(GTF_TSS_table_res$Index,PtoG_index)
        GTF_TSS_table_res$PtoG_cor = PtoG$Correlation[m]
        GTF_TSS_table_res$PtoG_FDR = PtoG$FDR[m]
        #######
        k = which(is.na(GTF_TSS_table_res$PtoG_FDR) == T)
        GTF_TSS_table_res$PtoG_FDR[k] = "ND"
        GTF_TSS_table_res$PtoG_cor[k] = "ND"
        GTF_TSS_table_res_cl = GTF_TSS_table_res
        #######
        ####### Next see the distance to the nearnest TSS gene ########
        #######
        #######
        #######
        GTF_TSS_table_res_cl$DistanceTSS = 0
        #######
        return(GTF_TSS_table_res_cl)
}

Process_Inter_peaks <- function(GTF_TSSBody_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
        New_table = data.frame(Gene = PtoG$gene, Peak = PtoG$peak, Class="Inter",Index=paste0(PtoG$gene,"#",PtoG$peak),PtoG_cor=PtoG$Correlation,PtoG_FDR=PtoG$FDR)
        ###### 
        ###### next remove peaks from GTF_TSSBody ######
        k = which(New_table$Peak %in% GTF_TSSBody_table_res$Peak == T)
        ######
        New_table_cl = New_table[-k,]
        ######
        ###### First add PtoG link correaltions !!! ###########
        ######
        #######
        ####### Next see the distance to the nearnest TSS gene ########
        #######
        New_table_cl_GR = GRanges(New_table_cl$Peak)
        #######
        dis_all = c()
        #######
        for(i in 1:length(New_table_cl_GR)){
                tmp = New_table_cl_GR[i]
                tmp_gene = New_table_cl$Gene[i]
                ######
                #print(tmp_gene)
                k = which(names(Gtfs_protein_trans_GR_TSS) == tmp_gene)
                if(length(k) > 0){
                    tmp_gene_tss = Gtfs_protein_trans_GR_TSS[[k]]
                    ######
                    tmp_dis = distance(tmp,tmp_gene_tss)
                    out = min(na.omit(tmp_dis))
                    dis_all = c(dis_all,out)
                }
                if(length(k) == 0){
                    dis_all = c(dis_all,"ND")
                }
                if(i %% 1000 == 0){
                        print(i)
                }
        }
        #######
        #######
        New_table_cl$DistanceTSS = dis_all
        #######
        k = which(New_table_cl$DistanceTSS == "ND")
        #######
        New_table_cl$DistanceTSS[k] = 10000000
        return(New_table_cl)
}

Merge_TSS_and_Body <- function(GTF_TSS_table_res,GTF_Body_table_res){
    ########
    GTF_TSS_table_res$Index = paste0(GTF_TSS_table_res$Gene,"#",GTF_TSS_table_res$Peak)
    GTF_Body_table_res$Index = paste0(GTF_Body_table_res$Gene,"#",GTF_Body_table_res$Peak)
    ########
    ######## rm Body with overlapped with tss ########
    ########
    k = which(GTF_Body_table_res$Index %in% GTF_TSS_table_res$Index == T)
    GTF_Body_table_res_cl = GTF_Body_table_res[-k,]
    ########
    GTF_TSSBody_table_res = rbind(GTF_TSS_table_res,GTF_Body_table_res_cl)
    ########
    return(GTF_TSSBody_table_res)
}

Process_peak_table_to_annotation <- function(ALL_peak_table){
    #######
    k1 = which(ALL_peak_table$Class %in% c("Body","Inter") == T & ALL_peak_table$DistanceTSS > 200000)
    ALL_peak_table_cl = ALL_peak_table[-k1,]
    #######
    k2 = which(ALL_peak_table_cl$Class == 'TSS' & ALL_peak_table_cl$PtoG_cor < 0)
    #######
    ALL_peak_table_cl2 = ALL_peak_table_cl[-k2,]
    ######
    k3 = which(ALL_peak_table_cl2$Class %in% c("Body","Inter") == T & ALL_peak_table_cl2$PtoG_cor < 0)
    ALL_peak_table_cl3 = ALL_peak_table_cl2[-k3,]
    ######
    ######
    ALL_peak_table_cl3$Index_Peak = paste0(ALL_peak_table_cl3$Gene,"::",ALL_peak_table_cl3$Class)
    ##### ALL_peak_table_cl3[which(ALL_peak_table_cl3$Gene == 'mmp9'),]
    #####
    ##### add tags for the peak #####
    ALL_peak_table_cl4 = split(ALL_peak_table_cl3,ALL_peak_table_cl3$Index_Peak)
    #####
    add_index <- function(x){
        x$Index_Peak2 = paste0(x$Index_Peak,"_",1:dim(x)[1])
        return(x)
    }
    #####
    ALL_peak_table_cl4 = lapply(ALL_peak_table_cl4,add_index)
    ALL_peak_table_cl4 = do.call("rbind",ALL_peak_table_cl4)
    #####
    return(ALL_peak_table_cl4)
}

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


Get_all_need_Motif_tag <- function(x,y,cutoff=0.05){
    x = rbind(x,y)
    ########
    k1 = which(x$Gene_Motif_Cor > cutoff | x$Gene_Motif_Cor < -cutoff)
    #########
    x_cl = x[k1,]
    #########
    k3 = which(x_cl$Gene_Motif_Cor > cutoff)
    k4 = which(x_cl$Gene_Motif_Cor < -cutoff)
    x_cl$class = 'Unknown'
    x_cl$class[k3] = 'pos'
    x_cl$class[k4] = 'neg'
    #########
    print(length(x_cl$ID[!duplicated(x_cl$ID)]))
    x_cl[which(x_cl$Gene == 'ascl1a'),]
    return(x_cl)
}


Globle_need_Motif = Get_all_need_Motif_tag(LD_NMDA_TF_activity,InjuryDev_TF_activity,0.05)

####
setwd("/zp1/data/plyu3/Fish_muti_omic_dev/merge_ATAC")

load(file='Dev_PtoG_Ann_July23')
load(file='Injury_PtoG_Ann_July23')
load(file='LD_PtoG_Ann_July23')
load(file='NMDA_PtoG_Ann_July23')

x1 = Dev_PtoG_Ann
x2 = Injury_PtoG_Ann
x3 = LD_PtoG_Ann
x4 = NMDA_PtoG_Ann

Get_all_need_Peak_tag <- function(x1,x2,x3,x4){
    all_tab = rbind(x1,x2,x3,x4)
    ######
    all_peak = all_tab$Peak
    all_peak = all_peak[!duplicated(all_peak)]
    ######
    return(all_peak)
}

########
######## merge all peak regions to match TF binding mtoifs ##########
######## Globle_need_Peak is the union regions of all these peaks ####
Globle_need_Peak = Get_all_need_Peak_tag(x1,x2,x3,x4)
######## 


setwd("/zp1/data/plyu3/Fish_muti_omic_dev/merge_ATAC")

#### Fish_combined_motifs is the PWMList which stored Fish PMW matrix #####
load(file='Fish_combined_motifs') ### Fish_combined ####

#### Fish_combined
#### PWMatrixList of length 9155
#### names(9155): M00001 M00002 M00004 ... M11389_2.00 M11390_2.00 M11491_2.00 

Merge_peaks <- function(x,extend = 150){
    library(GenomicRanges)
    x = GRanges(x)
    x = split(x,seqnames(x))
    ######
    for(i in 1:length(x)){
        tmp_x = x[[i]]
        tmp_x = sort(tmp_x)
        start(tmp_x) = start(tmp_x) - extend
        end(tmp_x) = end(tmp_x) + extend
        #####
        start(tmp_x[1]) = start(tmp_x[1]) + extend
        end(tmp_x[length(tmp_x)]) = end(tmp_x[length(tmp_x)]) - extend
        ######
        x[[i]] = tmp_x
    }
    return(x)
}

Merge_GRlist <- function(x_list){
    out = x_list[[1]]
    for(i in 2:length(x_list)){
        out = c(out,x_list[[i]])
    }
    return(out)
}


Match_motifs <- function(Globle_need_Peak,Globle_need_Motif,Fish_combined){
    ############
    Globle_need_Peak_1 = Merge_peaks(Globle_need_Peak,extend=150)
    ############
    ####Globle_need_Peak_2 = Merge_GRlist(Globle_need_Peak_1)
    ############ filter the Fish_combined #####
    k = which(names(Fish_combined) %in% Globle_need_Motif$ID == T)
    Fish_combined_cl = Fish_combined[k]
    ############
    library(parallel)
    cl <- makeCluster(35)
    ############
    footprint_Motif_list = parLapply(cl,Fish_combined_cl,Find_Motif,Globle_need_Peak_2=Globle_need_Peak_2)
    stopCluster(cl)
    ############
    head(names(footprint_Motif_list))
    ###########
    return(footprint_Motif_list)
    ###########
}

Find_Motif <- function(x,Globle_need_Peak_2){
    library(motifmatchr)
    library("BSgenome.Drerio.UCSC.danRer11")
    footprint = matchMotifs(x,Globle_need_Peak_2,genome = BSgenome.Drerio.UCSC.danRer11,out='positions',p.cutoff = 5e-05)
    ######
    footprint_tab = data.frame(footprint[[1]])
    ######
    return(footprint_tab)
}

#### 
Globle_motif_footprint <- Match_motifs(Globle_need_Peak,Globle_need_Motif,Fish_combined)
####
save(Globle_motif_footprint,file='Globle_motif_footprint_July23')

#### Next using LD or NMDA footprint to filter motif binding sites ######
#### Next load the corrected ATACseq signal calculated by TOBIAS ########

#### The pipline of TOBIAS please see: https://github.com/Pinlyu3/IReNA-v2 
#### STEP4: Predicting cell-type specific TFs binding in cis-regulatory elements
 
#### loading corrected signal by rtracklayer #####
NMDA_injury_footprint <- rtracklayer::import.bw("NMDA_injury_MG_fragments_GR_bam_tab_s_corrected.bw")
LD_injury_footprint <- rtracklayer::import.bw("LD_injury_MG_fragments_GR_bam_tab_s_corrected.bw")

####
Get_footprint_Score_all_para <- function(Globle_motif_footprint,footprint_signal){
    library(GenomicRanges)
    chr = sapply(Globle_motif_footprint,function(x) x$seqnames)
    start = sapply(Globle_motif_footprint,function(x) x$start)
    end = sapply(Globle_motif_footprint,function(x) x$end)
    score = sapply(Globle_motif_footprint,function(x) x$score)
    #########
    len = sapply(chr,function(x) length(x))
    names_tag = rep(names(len),len)
    #########
    chr =do.call(c,chr)
    start=do.call(c,start)
    end=do.call(c,end)
    score=do.call(c,score)
    #########
    new_GR = GRanges(seqnames=chr,IRanges(start,end),score=score,motif=names_tag)
    #########
    x_GRlist <- Get_left_mid_right_region_list(new_GR)
    x_left = cal_footprint_score_sub1(x_GRlist[[1]],footprint_signal)
    x_mid = cal_footprint_score_sub1(x_GRlist[[2]],footprint_signal)
    x_right = cal_footprint_score_sub1(x_GRlist[[3]],footprint_signal)
    #########
    x_res = Merge_footprint_score(x_left,x_mid,x_right)
    #########
    return(x_res)
}

cal_footprint_score_sub1 <- function(tmp_GR,tmp_F){
    library(GenomicRanges)
    ####
    res = findOverlaps(tmp_GR,tmp_F)
    ####
    score = tmp_F$score[subjectHits(res)]
    ####
    res_score_merge = tapply(score,queryHits(res),sum)
    ####
    tmp_GR$footprint_score = 0
    ####
    tmp_GR$footprint_score[as.numeric(names(res_score_merge))] = as.numeric(res_score_merge)
    ####
    tmp_GR$wid = width(tmp_GR)
    ####
    tmp_GR$footprint_score_norm = tmp_GR$footprint_score / tmp_GR$wid
    ####
    return(tmp_GR)
}

Get_left_mid_right_region_list <- function(x){
    library(GenomicRanges)
    #######
    x_L = x
    x_R = x
    x_M = x
    #######
    width = width(x)
    ######
    start(x_L) = start(x_L) - width*3
    end(x_L) = end(x_L) - width
    ########
    start(x_R) = start(x_R) + width
    end(x_R) = end(x_R) + width*3
    ########
    #######
    return(list(x_L,x,x_R))
}

Merge_footprint_score <- function(x_L,x_M,x_R){
    ###
    #############
    x = x_M
    #############
    x$left_score = 0
    x$mid_score = 0
    x$right_score = 0
    ##############
    x$left_delta = 0
    x$right_delta = 0
    #############
    x$left_score = x_L$footprint_score_norm
    x$mid_score = x_M$footprint_score_norm
    x$right_score = x_R$footprint_score_norm
    #############
    x$left_delta = x$left_score - x$mid_score
    x$right_delta = x$right_score - x$mid_score
    #############
    ## summary(x$mid_score)
    return(x)
    ###
}

#####
#####
#####

Globle_motif_cleanedByLD = Get_footprint_Score_all_para(Globle_motif_footprint,footprint_signal=LD_injury_footprint)
save(Globle_motif_cleanedByLD,file='Globle_motif_cleanedByLD')

Globle_motif_cleanedByNMDA = Get_footprint_Score_all_para(Globle_motif_footprint,footprint_signal=NMDA_injury_footprint)
save(Globle_motif_cleanedByNMDA,file='Globle_motif_cleanedByNMDA')

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

#### LD_network and NMDA_network are in the shared links ####
```



### STEP5: Construction of TF-peak-target links
``` r
#### load footprint, motif activity, peak-gene, and gene-gene correlations
Combined_all_things <- function(GG_network,PG_cor,TP_cor,Foot){
    library(GenomicRanges)
    ##### first filtr the GG_network #######
    ##### keep 95% intervals ###############
    GG_network = GG_network[order(GG_network$Score,decreasing=T),]
    GG_cutoff = quantile(GG_network$Score,0.9)
    #####
    GG_network_cl = GG_network[which(GG_network$Score > GG_cutoff),]
    #####
    k1 = which(GG_network_cl$Corr > 0.03)
    k2 = which(GG_network_cl$Corr < -0.03)
    #####
    GG_network_cl$Class = 'unknown'
    GG_network_cl$Class[k1] = 'pos'
    GG_network_cl$Class[k2] = 'neg'
    #####
    GG_network_clcl = GG_network_cl[which(GG_network_cl$Class != 'unknown'),]
    GG_network_clcl$TF_Gene = paste0(GG_network_clcl$TF,'->',GG_network_clcl$Gene)
    GG_network_clclcl = GG_network_clcl[,c('TF_Gene','Score','Corr')]
    colnames(GG_network_clclcl) = c('TF_Gene','GG_Score','GG_Corr')
    #####
    ## head(GG_network_clcl[which(GG_network_clcl$TF == 'foxj1a'),])
    #####
    ##### Next filter the footprint #####
    #####
    k_foot <- which(c(Foot$left_delta + Foot$right_delta) > 0.1)
    ######
    Foot_cl = Foot[k_foot]
    ######
    ###### Next link footprint to Peaks !!! ###########
    ######
    Footprint_Peak_table = All_footprint_to_peaks(Foot_cl,PG_cor)
    ###### Next merge the footpint peak to motif gene list !!! ##########
    ######
    TP_cor_cl = TP_cor[which(TP_cor$Gene_Motif_Cor > 0.05 | TP_cor$Gene_Motif_Cor < -0.05),]
    ######
    TP_cor_cl = TP_cor_cl[,c('Gene','ID','Gene_Motif_Cor')]
    colnames(TP_cor_cl) <- c("TF","Motif","TF_Motif_Cor")
    ###### Next we will merge big table ########
    library(data.table)
    library(dplyr)
    Footprint_Peak_table = data.table(Footprint_Peak_table)
    TP_cor_cl = data.table(TP_cor_cl)
    merged_table <- inner_join(Footprint_Peak_table, TP_cor_cl, by = "Motif")
    ###### test merge ##### correct !! #
    #Footprint_Peak_table = data.frame(Footprint_Peak_table)
    #TP_cor_cl = data.frame(TP_cor_cl)
    #merged_table_test <- merge(Footprint_Peak_table, TP_cor_cl)
    #dim(merged_table_test)
    ####### then next we will merge TF_gene corrections ######
    merged_table$TF_Gene = paste0(merged_table$TF,"->",merged_table$Target)
    ######
    GG_network_clclcl <- data.table(GG_network_clclcl)
    ######
    merged_table2 = inner_join(merged_table, GG_network_clclcl, by = "TF_Gene")
    ######
    k_filter1 = which(merged_table2$TF_Motif_Cor > 0 & merged_table2$GG_Corr > 0)
    k_filter2 = which(merged_table2$TF_Motif_Cor < 0 & merged_table2$GG_Corr < 0)
    #######
    merged_table3 = merged_table2[c(k_filter1,k_filter2),]
    #######
    print(dim(merged_table3))
    return(merged_table3)
}


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
Convert_to_triple_and_filterbygenes <- function(GRNs,Genes_need){
    ###########
    GRNs_cl = GRNs[,c('TF','Motif','Peak','Target','TF_Motif_Cor','PtoG_Class','GG_Corr')]
    ###########
    GRNs_cl_Triple = paste0(GRNs_cl$TF,"->",GRNs_cl$Peak,"->",GRNs_cl$Target)
    ###########
    GRNs_cl$Triple = GRNs_cl_Triple
    GRNs_cl_Triple_cl = GRNs_cl[!duplicated(GRNs_cl$Triple),]
    ###########
    k = which(GRNs_cl_Triple_cl$TF %in% Genes_need$gene_short == T & GRNs_cl_Triple_cl$Target %in% Genes_need$gene_short == T)
    ###########
    GRNs_cl_Triple_clcl = GRNs_cl_Triple_cl[k,]
    ############
    print(dim(GRNs_cl_Triple_clcl))
    return(GRNs_cl_Triple_clcl)
}

LD_GRNs_MG_triple = Convert_to_triple_and_filterbygenes(LD_GRNs,Genes_need)
NMDA_GRNs_MG_triple = Convert_to_triple_and_filterbygenes(NMDA_GRNs,Genes_need)



#### Next convert to double pairs ######
#### only kept the genes expressed in MG group cells ####

Convert_to_double_and_filterbygenes <- function(GRNs,Genes_need){
    ###########
    GRNs_cl = GRNs[,c('TF','Motif','Peak','Target','TF_Motif_Cor','PtoG_Class','GG_Corr')]
    ###########
    GRNs_cl_Double = paste0(GRNs_cl$TF,"->",GRNs_cl$Target)
    ###########
    GRNs_cl$Double = GRNs_cl_Double
    GRNs_cl_Double_cl = GRNs_cl[!duplicated(GRNs_cl$Double),]
    ###########
    k = which(GRNs_cl_Double_cl$TF %in% Genes_need$gene_short == T & GRNs_cl_Double_cl$Target %in% Genes_need$gene_short == T)
    ###########
    GRNs_cl_Double_clcl = GRNs_cl_Double_cl[k,]
    ############
    print(dim(GRNs_cl_Double_clcl))
    return(GRNs_cl_Double_clcl)
}

LD_GRNs_MG_double = Convert_to_double_and_filterbygenes(LD_GRNs,Genes_need)
NMDA_GRNs_MG_double = Convert_to_double_and_filterbygenes(NMDA_GRNs,Genes_need)


#####
##### Next load the log2 fold change files between LD and NMDA #######
#####

##### Differential DARs between LD and NMDA
load("DARs_Combined_LDvsNMDA_July23")

##### Differential DEGs between LD and NMDA
load('DEGs_Combined_LDNMDA_July23')

##### call specific GRNs by the fold change of Genes and Peaks ######

Filter_triple_network_according_to_foldchange <- function(GRNs_MG_triple,DEGs,DARs,tag='pos'){
    if(tag == 'pos'){
        DEGs=DEGs
        DARs=DARs
    }
    if(tag == 'neg'){
        DEGs$avg_log2FC = -DEGs$avg_log2FC
        DARs$Log2FC = -DARs$Log2FC
    }
    ######## ad activator tag ######
    GRNs_MG_triple$type = 'Act'
    k = which(GRNs_MG_triple$GG_Corr < 0)
    GRNs_MG_triple$type[k] = 'Rep'
    ########
    ######## for Act #######
    ########
    GRNs_MG_triple_Act = GRNs_MG_triple[which(GRNs_MG_triple$type == 'Act'),]
    GRNs_MG_triple_Rep = GRNs_MG_triple[which(GRNs_MG_triple$type == 'Rep'),]
    ######## for Act #######
    genes_1 = DEGs$gene_short[which(DEGs$avg_log2FC < 0)]
    genes_2 = DEGs$gene_short[which(DEGs$avg_log2FC > 0)]
    peaks = DARs$peaks[which(DARs$Log2FC < 0)]
    k = which(GRNs_MG_triple_Act$TF %in% genes_2 == T & GRNs_MG_triple_Act$Target %in% genes_2 == T & GRNs_MG_triple_Act$Peak %in% peaks == F)
    GRNs_MG_triple_Act_cl = GRNs_MG_triple_Act[k,]
    ######## for Rep #######
    genes_1 = DEGs$gene_short[which(DEGs$avg_log2FC < 0)]
    genes_2 = DEGs$gene_short[which(DEGs$avg_log2FC > 0)]
    peaks = DARs$peaks[which(DARs$Log2FC < 0)]
    k = which(GRNs_MG_triple_Rep$TF %in% genes_1 == T & GRNs_MG_triple_Rep$Target %in% genes_2 == T & GRNs_MG_triple_Rep$Peak %in% peaks == F)
    GRNs_MG_triple_Rep_cl = GRNs_MG_triple_Rep[k,]
    #######
    #######
    DEGs_pos = DEGs[which(DEGs$avg_log2FC > 0),]
    DEGs_neg = DEGs[which(DEGs$avg_log2FC < 0),]
    Peak_pos = DARs[which(DARs$Log2FC > 0),]
    ########
    ######## for GRNs_MG_triple_Act_cl #####
    GRNs_MG_triple_Act_cl$TF_fold = sapply(as.list(GRNs_MG_triple_Act_cl$TF),Get_highest_score_G,tab=DEGs_pos)
    GRNs_MG_triple_Act_cl$Target_fold = sapply(as.list(GRNs_MG_triple_Act_cl$Target),Get_highest_score_G,tab=DEGs_pos)
    GRNs_MG_triple_Act_cl$Peak_fold = sapply(as.list(GRNs_MG_triple_Act_cl$Peak),Get_highest_score_P,tab=Peak_pos)
    ########
    GRNs_MG_triple_Rep_cl$TF_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$TF),Get_lowest_score_G,tab=DEGs_neg)
    GRNs_MG_triple_Rep_cl$Target_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$Target),Get_highest_score_G,tab=DEGs_pos)
    GRNs_MG_triple_Rep_cl$Peak_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$Peak),Get_highest_score_P,tab=Peak_pos)
    ########
    GRNs_MG_triple_cl_out = rbind(GRNs_MG_triple_Act_cl,GRNs_MG_triple_Rep_cl)
    ########
    return(GRNs_MG_triple_cl_out)
}

Get_highest_score_G <- function(x,tab){
    k = which(tab$gene_short == x)
    #####
    if(length(k) == 0){
        out = 0
        return(out)
    }
    if(length(k) > 0){
        out = max(tab$avg_log2FC[k])
        return(out)
    }
}


Get_highest_score_P <- function(x,tab){
    k = which(tab$peaks == x)
    #####
    if(length(k) == 0){
        out = 0
        return(out)
    }
    if(length(k) > 0){
        out = max(tab$Log2FC[k])
        return(out)
    }
}


Get_lowest_score_G <- function(x,tab){
    k = which(tab$gene_short == x)
    #####
    if(length(k) == 0){
        out = 0
        return(out)
    }
    if(length(k) > 0){
        out = min(tab$avg_log2FC[k])
        return(out)
    }
}
LD_GRNs_MG_triple_specific <- Filter_triple_network_according_to_foldchange(LD_GRNs_MG_triple,DEGs,DARs,tag='pos')
NMDA_GRNs_MG_triple_specific <- Filter_triple_network_according_to_foldchange(NMDA_GRNs_MG_triple,DEGs,DARs,tag='neg')

##### 
library(writexl)
write_xlsx(LD_GRNs_MG_triple_specific, path = "LD_GRNs_MG_triple_specific.xlsx")
write_xlsx(NMDA_GRNs_MG_triple_specific, path = "NMDA_GRNs_MG_triple_specific.xlsx")
```




### STEP7: Identication of key activator TFs
``` r
#### 
load("DEG_Combined_LDNMDA_plot_July23")

Gene_list = data.frame(Gene=rownames(DEG_Combined_LDNMDA_plot[[1]]),Cluster=DEG_Combined_LDNMDA_plot[[3]])
Gene_list$Gene_short = sapply(strsplit(Gene_list$Gene,split='~~'),function(x) x[[2]])
Gene_list$ID = sapply(strsplit(Gene_list$Gene,split='~~'),function(x) x[[1]])

table(Gene_list$Cluster)

LDNMDA_DEG_list = Gene_list

Get_Key_Activate_Factors <- function(all_network,specific_network,gene_groups){
    ###########
    specific_network_pos = specific_network[which(specific_network$type=='Act'),]
    all_observed_TF = names(table(specific_network$TF))
    ###########
    res = list()
    ############
    for(i in 1:length(all_observed_TF)){
        tmp_TF = all_observed_TF[i]
        #######
        tmp_TF_all_network = all_network[which(all_network$TF == tmp_TF),]
        tmp_TF_sp_network = specific_network[which(specific_network$TF == tmp_TF),]
        ########
        tmp_TF_all_target = tmp_TF_all_network$Target[!duplicated(tmp_TF_all_network$Target)]
        tmp_TF_sp_target = tmp_TF_sp_network$Target[!duplicated(tmp_TF_sp_network$Target)]
        ########
        ######## OK!! for each groups !!! #####
        ########
        Cluster = as.numeric(names(table(as.numeric(gene_groups$Cluster))))
        ########
        for(j in 1:length(Cluster)){
            tmp_Cluster = gene_groups[which(gene_groups$Cluster == Cluster[j]),]
            #####
            tmp_Cluster_gene = tmp_Cluster$Gene_short
            #####
            overlap = which(tmp_TF_sp_target %in% tmp_Cluster_gene == T)
            #####
            overlap_all = which(tmp_TF_all_target %in% tmp_Cluster_gene == T)
            #####
            tab = data.frame(TF=tmp_TF,All_Target=length(tmp_TF_all_target),Sp_Target=length(tmp_TF_sp_target),Cluster=Cluster[j],Sp_Cluster_overlap=length(overlap),All_Cluster_overlap=length(overlap_all),Cluster_size=length(tmp_Cluster_gene))
            ####
            res = c(res,list(tab))
        }
        #########
        #########
    }
    ###########
    res = do.call('rbind',res)
    res$coverage = res$Sp_Cluster_overlap / res$Cluster_size
    ###########
    res[which(res$TF=='foxj1a'),]
    ########### Next we will calculate the p-value ########
    ########### speicific target to that cluster ##########
    ###########
    res$pvalue = apply(res,1,perform_phyper)
    ###########
    return(res)
}

perform_phyper <- function(x){
    ####
    population_size = as.numeric(x[2])
    successes_in_population = as.numeric(x[6])
    successes_in_sample = as.numeric(x[5])
    sample_size = as.numeric(x[3])
    ####
    p_value <- phyper(successes_in_sample - 1, successes_in_population, population_size - successes_in_population, sample_size, lower.tail = FALSE)
    ####
    return(p_value)
}

all_network = LD_GRNs_MG_triple
specific_network = LD_GRNs_MG_triple_specific
gene_groups = LDNMDA_DEG_list[which(LDNMDA_DEG_list$Cluster %in% c(1,2,3) == T),]


LD_key_TF = Get_Key_Activate_Factors(all_network,specific_network,gene_groups)

prepare_plot <- function(res){
    res$log10_P = -log10(res$pvalue)
    res[which(res$TF=='foxj1a'),]
    #####
    res_list = split(res,res$TF)
    #####
    res_list$foxj1a
    #####
    new_list = list()
    for(i in 1:length(res_list)){
        res_list_1 = res_list[[i]]
        if(max(res_list_1$log10_P > 3 & max(res_list_1$coverage > 0.01))){
            new_list <- c(new_list,list(res_list_1))
        }
    }
    #####
    new_list = do.call('rbind',new_list)
    #####
    ##### next we will rank these TFs by cluster ########
    #####
    k = which(is.infinite(new_list$log10_P) == TRUE)
    new_list$log10_P[k] = 100
    ##### rank !!! #######
    new_list_2 = new_list[order(new_list$log10_P,decreasing=T),]
    rm_dup = which(duplicated(new_list_2$TF) == T)
    new_list_2 = new_list_2[-rm_dup,]
    ####
    GO_res_sub2_sp = split(new_list_2,new_list_2$Cluster)
    ####
    final_order = c()
    for(j in 1:length(GO_res_sub2_sp)){
        final_order = c(final_order,GO_res_sub2_sp[[j]]$TF)
    }
    ####
    new_list$TF = factor(new_list$TF,levels=final_order)
    ####
    return(new_list) 
}

LD_key_TF_plot = prepare_plot(LD_key_TF)

LD_key_TF_plot$Cluster = factor(LD_key_TF_plot$Cluster,levels=c(3,2,1))

k = which(LD_key_TF_plot$coverage > 0.2)
LD_key_TF_plot$coverage[k] = 0.2

#######
library(ggplot2)

ggplot(LD_key_TF_plot,aes(x=TF,y=Cluster)) + geom_point(aes(size=log10_P,color=coverage)) + theme_classic() +
xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),panel.background = element_rect(colour = "black",size = 1, linetype = "solid")) + scale_size(range = c(0,8),breaks = c(0,2,4,6)) +
scale_color_gradient2(low="grey",mid = "#A74997",high="darkred",limits = c(0,0.2),midpoint=0.1) + guides(size = guide_legend(order = 2))

ggsave("LD.pdf",height=3,width=10)

```





