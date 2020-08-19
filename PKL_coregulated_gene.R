# PKL co-regulated genes detection
# PKLR inhibited/over-expressed hepg2 RNA-seq data


# z score transferation and combination

library(ggplot2)
library(psych)

res_ov <- as.data.frame(read.table("res_PKLR_ov.csv",
                                   header = T,
                                   sep = ",",
                                   stringsAsFactors = F))

res_si <- as.data.frame(read.table("res_PKLR_si.csv",
                                   header = T,
                                   sep = ",",
                                   stringsAsFactors = F))

res_ov_up <- na.omit(subset(res_ov,res_ov$log2FoldChange > 0))
res_ov_down <- na.omit(subset(res_ov,res_ov$log2FoldChange < 0))
res_si_up <- na.omit(subset(res_si,res_si$log2FoldChange > 0))
res_si_down <- na.omit(subset(res_si,res_si$log2FoldChange < -0))

res_ovz <- qnorm(res_ov_up$pvalue/2)
res_ovz [which (res_ovz == "Inf")] = -100
res_ovz [which (res_ovz == -100)] = max(res_ovz)
res_ovz [which (res_ovz == "-Inf")] = 100
res_ovz [which (res_ovz == 100)] = min(res_ovz)
res_ovz <- -res_ovz

res_ov_up <- cbind(res_ov_up,res_ovz)

res_ovz <- qnorm(res_ov_down$pvalue/2)
res_ovz [which (res_ovz == "Inf")] = -100
res_ovz [which (res_ovz == -100)] = max(res_ovz)
res_ovz [which (res_ovz == "-Inf")] = 100
res_ovz [which (res_ovz == 100)] = min(res_ovz)

res_ov_down <- cbind(res_ov_down,res_ovz)
# inverse the value of negative correlated genes for ov

res_siz <- qnorm(res_si_up$pvalue/2)
res_siz [which (res_siz == "Inf")] = -100
res_siz [which (res_siz == -100)] = max(res_siz)
res_siz [which (res_siz == "-Inf")] = 100
res_siz [which (res_siz == 100)] = min(res_siz)

res_si_up <- cbind(res_si_up,res_siz)
# inverse the value of positive correlated genes for si

res_siz <- qnorm(res_si_down$pvalue/2)
res_siz [which (res_siz == "Inf")] = -100
res_siz [which (res_siz == -100)] = max(res_siz)
res_siz [which (res_siz == "-Inf")] = 100
res_siz [which (res_siz == 100)] = min(res_siz)
res_siz <- -res_siz
res_si_down <- cbind(res_si_down,res_siz)


res_ov_z <- rbind(res_ov_up, res_ov_down)
res_si_z <- rbind(res_si_up, res_si_down)



rownames(res_ov_z) <- 1:nrow(res_ov_z)
rownames(res_si_z) <- 1:nrow(res_si_z)


res_ov_z <- res_ov_z[order(res_ov_z$baseMean),]
res_si_z <- res_si_z[order(res_si_z$baseMean),]

res_ov_z <- res_ov_z[ !duplicated(res_ov_z$gene_symbol, fromLast = TRUE),]
res_si_z <- res_si_z[ !duplicated(res_si_z$gene_symbol, fromLast = TRUE),]



#spearman results with p value
cor_up <- as.data.frame(read.table("F:/sweden/master_thesis/data_analysis/GTExCluster/PKL_upCor_v7.csv",
                                   header = T,
                                   stringsAsFactors = F,
                                   sep =","))

cor_down <- as.data.frame(read.table("F:/sweden/master_thesis/data_analysis/GTExCluster/PKL_downCor_v7.csv",
                                     header = T,
                                     stringsAsFactors = F,
                                     sep =","))



total_gene_temp <- as.data.frame(read.csv("com_codeGene_est.csv",
                                          row.names = 1,
                                          stringsAsFactors = F,
                                          header = T))

total_gene <- as.data.frame(total_gene_temp$gene_symbol)
row.names(total_gene) <- row.names(total_gene_temp)


corz <- qnorm(cor_up$PKL_pval/2)
corz [which (corz == "Inf")] = -100
corz [which (corz == -100)] = max(corz)
corz [which (corz == "-Inf")] = 100
corz [which (corz == 100)] = min(corz)
corz <- -corz
corz <- data.frame(corz)
cor_up <- cbind(cor_up, corz)

corz <- qnorm(cor_down$PKL_pval/2)
corz [which (corz == "Inf")] = -100
corz [which (corz == -100)] = max(corz)
corz [which (corz == "-Inf")] = 100
corz [which (corz == 100)] = min(corz)

corz <- data.frame(corz)
cor_down <- cbind(cor_down, corz)
#inverse the negative correlated genes value for spearman

cor_z <- rbind(cor_up, cor_down)


row.names( cor_z) <- cor_z[,1]
cor_z <- merge(cor_z,
               total_gene,
               by ="row.names",
               sort = F)

cor_z <- cor_z[,-1]


cor_z <- cor_z[order(cor_z$mean_counts),]
cor_z <- cor_z[!duplicated(cor_z$`total_gene_temp$gene_symbol`, fromLast = TRUE),]

write.csv(res_ov_z, file = "PKLR_ovz.csv", row.names = F)
write.csv(res_si_z, file = "PKLR_siz.csv", row.names = F)
write.csv(cor_z, file= "PKLR_corz.csv", row.names = F)


path_cor = "F:\\sweden\\master_thesis\\data_analysis\\RNAquantify\\DEG\\differential_expressed_analysis\\"

# PKLR_ovz.csv, PKLR_siz.csv, PKLR_corz.csv, PKL_corGene.csv
for (i in c("ovz","siz","corz")){
  PKL_file <- paste(path_cor,"PKLR_",i,".csv", sep="")
  
  PKLR_gene <- as.data.frame(read.table(PKL_file,
                                        sep = ",",
                                        header =T,
                                        stringsAsFactors = F,
                                        row.names = 6))
  
  
  
  cmap_gene <- as.data.frame(read.csv("GSE92742_Broad_LINCS_gene_info.txt",
                                      sep = "\t",
                                      row.names = 2,
                                      header = T))
  
  cmap_gene_num <- as.data.frame(cmap_gene$pr_gene_id)
  row.names(cmap_gene_num) <- row.names(cmap_gene)
  
  PKLR_gene <- merge(PKLR_gene, 
                     cmap_gene_num,
                     by = "row.names",
                     sort=F)
  
  PKLR_gene = PKLR_gene[,-2]
  write.csv(PKLR_gene, file="PKL_",i,"_num.csv", row.names = F)
}



# merge the data with gene number for cmap
PKL_si_num <- as.data.frame(read.csv("PKL_siz_num.csv",
                                     header =T,
                                     row.names = 1,
                                     stringsAsFactors = F))


PKL_ov_num <- as.data.frame(read.csv("PKL_ovz_num.csv",
                                     header =T,
                                     row.names = 1,
                                     stringsAsFactors = F))


PKL_cor_num <- as.data.frame(read.csv("PKL_corz_num.csv",
                                      header=T,
                                      row.names = 1,
                                      stringsAsFactors = F))


PKL_ovsi_num <- merge(PKL_si_num,
                      PKL_ov_num,
                      by="row.names",
                      sort=F)

gene_name <- PKL_ovsi_num$Row.names
gene_num <- PKL_ovsi_num$cmap_gene.pr_gene_id.y
ov_z <- PKL_ovsi_num$res_ovz
si_z <- PKL_ovsi_num$res_siz

PKL_ovsi_num <- data.frame(gene_num, ov_z, si_z)
row.names(PKL_ovsi_num) <- gene_name



PKL_num <- merge(PKL_ovsi_num,
                 PKL_cor_num,
                 by="row.names",
                 sort = F)

cor_z <- PKL_num$corz
ov_z <- PKL_num$ov_z
si_z <- PKL_num$si_z

sum_z <- (ov_z+si_z+cor_z)/sqrt(3)
gene_name <- PKL_num$Row.names
gene_num <- PKL_num$gene_num

PKL_comCor <- data.frame(gene_name, gene_num, ov_z, si_z, cor_z, sum_z)
PKL_comCor_up <- subset(PKL_comCor, PKL_comCor$ov_z>=0 & PKL_comCor$si_z>=0 & PKL_comCor$cor_z>=0)
PKL_comCor_down <- subset(PKL_comCor, PKL_comCor$ov_z<0 & PKL_comCor$si_z<0 & PKL_comCor$cor_z<0)

PKL_comCor_overlap <- rbind(PKL_comCor_down, PKL_comCor_up)

write.csv(PKL_comCor_overlap, file="PKL_statistical_overlap_genes_strict.csv", row.names=FALSE, quote=FALSE)






