#### all functions ####

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


Get_all_need_Peak_tag <- function(x1,x2,x3,x4){
  all_tab = rbind(x1,x2,x3,x4)
  ######
  all_peak = all_tab$Peak
  all_peak = all_peak[!duplicated(all_peak)]
  ######
  return(all_peak)
}

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
  Globle_need_Peak_2 = Merge_GRlist(Globle_need_Peak_1)
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

All_footprint_to_peaks <- function(footprint_GR,peak_table){
  ####### extend to 801bp #################
  peak_table_GR = GRanges(peak_table$Peak)
  start(peak_table_GR) = start(peak_table_GR)-150
  end(peak_table_GR) = end(peak_table_GR)+150
  #######
  #######
  ####### Next we extend the peak_table #######
  ####### find the overlaps with peaks and footprint ####
  ####### left is peak right is the footprint ###
  countIndex = findOverlaps(peak_table_GR,footprint_GR)
  ####### 
  left_table = peak_table[queryHits(countIndex),]
  #######
  right_table = data.frame(footprint_GR)[subjectHits(countIndex),]
  ####### Clean_left_table ########
  left_table_cl = left_table[,c("Index","Peak","Class","PtoG_cor","Gene")]
  colnames(left_table_cl) <- c("PtoG_Index","Peak","PtoG_Class","PtoG_cor","Target")
  ### 
  k = which(left_table_cl$PtoG_cor == 'ND')
  left_table_cl$PtoG_cor[k] = 0
  left_table_cl$PtoG_cor = round(as.numeric(left_table_cl$PtoG_cor),3)
  ####### head(right_table)
  right_table$footprint = paste0(right_table$seqnames,':',right_table$start,'-',right_table$end)
  right_table_cl = right_table[,c('motif','footprint','mid_score','left_delta','right_delta')]
  colnames(right_table_cl) <- c("Motif",'Footprint_region','Footprint_Mid',"Footprint_Left_del","Footprint_Right_del")
  #######
  merge_table = cbind(left_table_cl,right_table_cl)
  #######
  return(merge_table)
}


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

