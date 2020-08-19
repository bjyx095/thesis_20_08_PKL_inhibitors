# PKL co-related gene signature
# meta-analysis to merge the gene list and combine probability

# strict combined data without/with gene num

corz <- read.csv("PKLR_corz.csv",
                 header= T,
                 row.names = 1,
                 stringsAsFactors = F)

siz<- read.csv("PKLR_siz.csv",
               header= T,
               row.names = 1,
               stringsAsFactors = F)

ovz<- read.csv("PKLR_ovz.csv",
               header= T,
               row.names = 1,
               stringsAsFactors = F)


# pkl_cor_v5 is calculated by bh methods, with alpha=0.1
corq <- read.csv("F:\\sweden\\master_thesis\\data_analysis\\GTExCluster\\PKL_cor_v5.csv",
                 header = T,
                 row.names = 1,
                 stringsAsFactors = F)

cor <- merge(corz, corq, by = "row.names", sort=F)
rownames(cor) <- cor$Row.names
ovsi_temp <- merge(siz,ovz, by="row.names", sort=F)
rownames(ovsi_temp)<- ovsi_temp$Row.names
ovsicor_temp <- merge(ovsi_temp, cor, by = "row.names", sort=F)
write.csv(ovsicor_temp,file="sparse_rawdata.csv", row.names = F)


sparse_com <- data.frame(ovsicor_temp$gene_symbol.x, ovsicor_temp$padj.x, ovsicor_temp$res_siz, ovsicor_temp$padj.y, ovsicor_temp$res_ovz, ovsicor_temp$PKL_pval.y,ovsicor_temp$corz)
colnames(sparse_com)<-c("gene_symbol","si_qval","siz","ov_qval","ovz","cor_qval","corz")
row.names(sparse_com)<- sparse_com$gene_symbol
sparse_com <- sparse_com[,-1]


sparse_com$sumz <- (sparse_com$siz+sparse_com$ovz+sparse_com$corz)/sqrt(3)
write.csv(sparse_com, file ="sparse_qzdata.csv")



# Hypergeometric distribution

# For each two interacted group, the hypergeometric distribution was tested to check if the shared genes make statistical sense or is by chance.


# overlapped genes for each two groups
overlap_posCor_ovsi <- intersect(res_ov_up_gene, res_si_down_gene)
overlap_posCor_ovcor <- intersect(res_ov_up_gene, row.names(cor_up))
overlap_posCor_sicor <- intersect(res_si_down_gene, row.names(cor_up))

overlap_negCor_ovsi <- intersect(res_ov_down_gene, res_si_up_gene)
overlap_negCor_ovcor <- intersect(res_ov_down_gene, row.names(cor_down))
overlap_negCor_sicor <- intersect(res_si_up_gene, row.names(cor_down))



overlap_posCor_ovsi_x <- as.numeric(length(overlap_posCor_ovsi)-1)
overlap_posCor_ovsi_m <- as.numeric(min(length(res_ov_up_gene),length(res_si_down_gene)))
overlap_posCor_ovsi_n <- as.numeric(totalGene_num - overlap_posCor_ovsi_m)
overlap_posCor_ovsi_k <- as.numeric(max(length(res_ov_up_gene),length(res_si_down_gene)))

overlap_posCor_ovsi_hyper <- phyper(overlap_posCor_ovsi_x,overlap_posCor_ovsi_m,overlap_posCor_ovsi_n,overlap_posCor_ovsi_k, lower.tail = FALSE)
overlap_posCor_ovsi_hyper_log <- phyper(overlap_posCor_ovsi_x,overlap_posCor_ovsi_m,overlap_posCor_ovsi_n,overlap_posCor_ovsi_k, lower.tail = FALSE,log.p = TRUE)



overlap_posCor_ovcor_x <- as.numeric(length(overlap_posCor_ovcor)-1)
overlap_posCor_ovcor_m <- as.numeric(min(length(res_ov_up_gene),length(rownames(cor_up))))
overlap_posCor_ovcor_n <- as.numeric(totalGene_num  - overlap_posCor_ovcor_m)
overlap_posCor_ovcor_k <- as.numeric(max(length(res_ov_up_gene),length(rownames(cor_up))))

overlap_posCor_ovcor_hyper <- phyper(overlap_posCor_ovcor_x,overlap_posCor_ovcor_m,overlap_posCor_ovcor_n,overlap_posCor_ovcor_k, lower.tail = FALSE)
overlap_posCor_ovcor_hyper_log <- phyper(overlap_posCor_ovcor_x,overlap_posCor_ovcor_m,overlap_posCor_ovcor_n,overlap_posCor_ovcor_k, lower.tail = FALSE,log.p = TRUE)


overlap_posCor_sicor_x <- as.numeric(length(overlap_posCor_sicor)-1)
overlap_posCor_sicor_m <- as.numeric(min(length(res_si_down_gene),length(rownames(cor_up))))
overlap_posCor_sicor_n <- as.numeric(totalGene_num - overlap_posCor_sicor_m)
overlap_posCor_sicor_k <- as.numeric(max(length(res_si_down_gene),length(rownames(cor_up))))

overlap_posCor_sicor_hyper <- phyper(overlap_posCor_sicor_x,overlap_posCor_sicor_m,overlap_posCor_sicor_n,overlap_posCor_sicor_k, lower.tail = FALSE)
overlap_posCor_sicor_hyper_log <- phyper(overlap_posCor_sicor_x,overlap_posCor_sicor_m,overlap_posCor_sicor_n,overlap_posCor_sicor_k, lower.tail = FALSE,log.p = TRUE)


overlap_negCor_ovsi_x <- as.numeric(length(overlap_negCor_ovsi)-1)
overlap_negCor_ovsi_m <- as.numeric(min(length(res_ov_down_gene),length(res_si_up_gene)))
overlap_negCor_ovsi_n <- as.numeric(totalGene_num  - overlap_negCor_ovsi_m)
overlap_negCor_ovsi_k <- as.numeric(max(length(res_ov_down_gene),length(res_si_up_gene)))


overlap_negCor_ovsi_hyper <- phyper(overlap_negCor_ovsi_x,overlap_negCor_ovsi_m,overlap_negCor_ovsi_n,overlap_negCor_ovsi_k, lower.tail = FALSE)
overlap_negCor_ovsi_hyper_log <- phyper(overlap_negCor_ovsi_x,overlap_negCor_ovsi_m,overlap_negCor_ovsi_n,overlap_negCor_ovsi_k, lower.tail = FALSE,log.p = TRUE)


overlap_negCor_ovcor_x <- as.numeric(length(overlap_negCor_ovcor)-1)
overlap_negCor_ovcor_m <- as.numeric(min(length(res_ov_down_gene),length(rownames(cor_down))))
overlap_negCor_ovcor_n <- as.numeric(totalGene_num  - overlap_negCor_ovcor_m)
overlap_negCor_ovcor_k <- as.numeric(max(length(res_ov_down_gene),length(rownames(cor_down))))

overlap_negCor_ovcor_hyper <- phyper(overlap_negCor_ovcor_x,overlap_negCor_ovcor_m,overlap_negCor_ovcor_n,overlap_negCor_ovcor_k, lower.tail = FALSE)
overlap_negCor_ovcor_hyper_log <- phyper(overlap_negCor_ovcor_x,overlap_negCor_ovcor_m,overlap_negCor_ovcor_n,overlap_negCor_ovcor_k, lower.tail = FALSE,log.p = TRUE)


overlap_negCor_sicor_x <- as.numeric(length(overlap_negCor_sicor)-1)
overlap_negCor_sicor_m <- as.numeric(min(length(res_si_up_gene),length(rownames(cor_down))))
overlap_negCor_sicor_n <- as.numeric(totalGene_num  - overlap_negCor_sicor_m)
overlap_negCor_sicor_k <- as.numeric(max(length(res_si_up_gene),length(rownames(cor_down))))

overlap_negCor_sicor_hyper <- phyper(overlap_negCor_sicor_x,overlap_negCor_sicor_m,overlap_negCor_sicor_n,overlap_negCor_sicor_k, lower.tail = FALSE)
overlap_negCor_sicor_hyper_log <- phyper(overlap_negCor_sicor_x,overlap_negCor_sicor_m,overlap_negCor_sicor_n,overlap_negCor_sicor_k, lower.tail = FALSE,log.p = TRUE)


hyperP <- data.frame(
  x = c(overlap_posCor_ovsi_x,overlap_posCor_sicor_x,overlap_posCor_ovcor_x,overlap_negCor_ovsi_x,overlap_negCor_sicor_x,overlap_negCor_ovcor_x),
  m = c(overlap_posCor_ovsi_m,overlap_posCor_sicor_m,overlap_posCor_ovcor_m,overlap_negCor_ovsi_m,overlap_negCor_sicor_m,overlap_negCor_ovcor_m),
  n = c(overlap_posCor_ovsi_n,overlap_posCor_sicor_n,overlap_posCor_ovcor_n,overlap_negCor_ovsi_n,overlap_negCor_sicor_n,overlap_negCor_ovcor_n),
  k = c(overlap_posCor_ovsi_k,overlap_posCor_sicor_k,overlap_posCor_ovcor_k,overlap_negCor_ovsi_k,overlap_negCor_sicor_k,overlap_negCor_ovcor_k),
  p = c(overlap_posCor_ovsi_hyper,overlap_posCor_sicor_hyper,overlap_posCor_ovcor_hyper,overlap_negCor_ovsi_hyper,overlap_negCor_sicor_hyper,overlap_negCor_ovcor_hyper),
  log_p = c(overlap_posCor_ovsi_hyper_log,overlap_posCor_sicor_hyper_log,overlap_posCor_ovcor_hyper_log,overlap_negCor_ovsi_hyper_log,overlap_negCor_sicor_hyper_log,overlap_negCor_ovcor_hyper_log),
  row.names = c("posCor_ovsi", "posCor_sicor","posCor_ovcor","negCor_ovsi", "negCor_sicor","negCor_ovcor")
)


write.csv(hyperP, file = "hyperP_PKL.csv")


