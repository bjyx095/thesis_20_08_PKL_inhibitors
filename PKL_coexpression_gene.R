# PKL co-expression genes detection in GTEx dataset (110 human liver)
library(psych)


# load data
df1 <- as.data.frame(read.table('Liver_exp_gene_level_rENSG.txt',
                                sep = '\t',
                                row.names = 1,
                                header = T))

df2 <- as.data.frame(read.table('PKLR_liver_trans_tpm_Gtex_v23.txt',
                                sep = '\t',
                                row.names = 1,
                                header =T))

# data preprocessing
df1_sub <- subset(df1, rowMeans(df1) >= 1 )
df2t <- as.data.frame(t(df2))
df1_subt <- as.data.frame(t(df1_sub))
df3 <- as.data.frame(t(rbind(df1_sub, df2)))

cor_pkl_v4 <- corr.test(x = df2t, y=df1_subt, method = "spearman", adjust = "none")


corMatrix_PKLR_v4_r <- as.data.frame(cor_pkl_v4$r)
corMatrix_PKLR_v4_p <- as.data.frame(cor_pkl_v4$p)
corMatrix_PKLR_v4_r <- as.data.frame(t(corMatrix_PKLR_v4_r))
corMatrix_PKLR_v4_p <- as.data.frame(t(corMatrix_PKLR_v4_p))


PKL_r <- corMatrix_PKLR_v4_r$`ENST00000392414(PKL)`
PKL_pval <- corMatrix_PKLR_v4_p$`ENST00000392414(PKL)`
corMatrix_PKL_v4 <- data.frame(PKL_r,PKL_pval)
row.names(corMatrix_PKL_v4) <- row.names(corMatrix_PKLR_v4_r)

write.csv(corMatrix_PKL_v4, file = "PKL_cor_v6.csv")

# add the level of transcript counts for further delete of duplicate genes 
mean_counts <- apply(df1_sub,1,mean)
df_1_mean <- as.data.frame(mean_counts)
row.names(df_1_mean) <- row.names(df1_sub)
corMatrix_PKL_v5 <- merge(corMatrix_PKL_v4,
                          df_1_mean,
                          by= "row.names",
                          sort = F)

write.csv(corMatrix_PKL_v5, file = "PKL_cor_v7.csv")

corMat_PKL_up_v5 <- subset(corMatrix_PKL_v5, corMatrix_PKL_v5$PKL_r>0)
corMat_PKL_down_v5 <- subset(corMatrix_PKL_v5, corMatrix_PKL_v5$PKL_r<0)

write.csv(corMat_PKL_down_v5, file = "PKL_downCor_v7.csv", row.names=F)
write.csv(corMat_PKL_up_v5, file = "PKL_upCor_v7.csv", row.names = F)


