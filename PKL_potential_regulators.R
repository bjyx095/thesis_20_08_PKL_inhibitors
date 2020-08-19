# potential drug repositioning
# PKL co-related gene signauture with cmap l1000 data

# PKL inhibitors drug repositioning

library(psych)
library(cmapR)

cellLine_ID <- as.data.frame(read.csv('cell_id.txt',
                                      header = F
))


# PKL_ovz_num.csv, PKL_siz_num.csv, PKL_corz_num.csv, PKL_corGene_num.csv


PKLR_gene_temp <- as.data.frame(read.csv("PKL_statistical_overlap_genes_strict.csv",
                                         sep = ",",
                                         header = T,
                                         row.names = 2))

PKLR_gene <- PKLR_gene_temp[,-1]

sum_cor_r = data.frame()

sum_cor_p = data.frame()

path_name = "F:\\sweden\\master_thesis\\data_analysis\\l1000\\cmap_drug_discovery\\cell_line_data\\"

for (i in cellLine_ID) {
  file_name = paste(path_name,i,".gctx",sep="")
}


corCal <- function(i,PKLR_gene_list, type){
  
  df <- parse_gctx(i, rid=row.names(PKLR_gene_list))
  
  cor <- corr.test(x =PKLR_gene, y=df@mat, use= "pairwise",method= "spearman" , adjust = "BH")
  
  r_file = paste(i,  "_", type, "_cor.csv", sep="")
  r= t(cor$r)
  p = t(cor$p)
  df_2 = merge(r, p, by = "row.names", sort = F)
  write.csv(df_2, file = r_file)
  
  return(r_file)
}

# use all cell line data

for (i in file_name){corCal(i, PKLR_gene, "combined")}


# except the drugs which r <-0.4 and not appear in other cell line


hepg2_com_strict <- read.csv("HEPG2.gctx_combined_cor_2.csv",
                             header = T,
                             row.names = 2,
                             stringsAsFactors = F)


sum_strict <- read.csv("sum_2.csv",
                       header = T,
                       row.names = 1,
                       stringsAsFactors = F)

col_meta <- as.matrix(read.csv("GSE70138_Broad_LINCS_sig_info.txt",
                               sep = "\t",
                               header =T,
                               stringsAsFactors = F,
                               row.names = 1))

pert_meta <- as.matrix (read.csv("GSE70138_Broad_LINCS_pert_info.txt",
                                 sep= "\t",
                                 header = T,
                                 stringsAsFactors = F,
                                 row.names = 1))


hepg2_com_strict$rank_sum <- rank(hepg2_com_strict$neg_Zscore.x)

col_name <- row.names(subset(hepg2_com_strict, hepg2_com_strict$res_ov_down_z.x< -0.4 | hepg2_com_strict$res_si_up_z.x < -0.4 | hepg2_com_strict$cor_down_z.x< -0.4 | hepg2_com_strict$neg_Zscore.x < -0.4))
col_meta_hepg2 <- col_meta[which (col_meta[,"cell_id"]=="HEPG2"),]

col_meta_hepg2_sub <- data.frame()
hepg2_drug <- c()

for (i in col_name) {
  hepg2_drug_temp <- col_meta_hepg2[which(row.names(col_meta_hepg2)==i),"pert_id"]
  if (length(hepg2_drug)<1){
    hepg2_drug <- hepg2_drug_temp
  } else{
    hepg2_drug <- c(hepg2_drug,hepg2_drug_temp)
  }
  
}

for (i in hepg2_drug){
  col_meta_hepg2_sub_temp <- as.data.frame(col_meta_hepg2[which(col_meta_hepg2[,"pert_id"]==i),])
  if (nrow(col_meta_hepg2_sub)<=1){
    col_meta_hepg2_sub <- col_meta_hepg2_sub_temp
  } else{
    col_meta_hepg2_sub <- rbind(col_meta_hepg2_sub, col_meta_hepg2_sub_temp)
    
  }
}


PKL_drug_sel_strict <- merge(col_meta_hepg2_sub,
                             hepg2_com_strict,
                             by= "row.names",
                             sort=F)

write.csv(PKL_drug_sel_strict, file="PKL_drug_hepg2_strict.csv", row.names = F)




