

codeGene_ID <- as.data.frame(read.table("codeGene_ID.csv",
                                        sep =",",
                                        row.names = 1,
                                        header = T,
                                        stringsAsFactors = F))

codeGene_symbol_temp <- as.data.frame(read.table("sep_code_ID.csv",
                                                 sep = ",",
                                                 row.names = 2,
                                                 header =T,
                                                 stringsAsFactors = F,
                                                 encoding = 'UTF-8'))

# create a dataframe containing matched gene ID and symbol
codeGene_symbol <- data.frame (gene_symbol = codeGene_symbol_temp$X.U.FEFF.gene_symbol)
row.names(codeGene_symbol) <- rownames(codeGene_symbol_temp)

codeGene_name <- merge(codeGene_ID,
                       codeGene_symbol,
                       by = "row.names",
                       sort = F)


sample_101 <- as.data.frame(read.table("./sep_est/P8717_101_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_102 <- as.data.frame(read.table("./sep_est/P8717_102_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_103 <- as.data.frame(read.table("./sep_est/P8717_103_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_104 <- as.data.frame(read.table("./sep_est/P8717_104_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_105 <- as.data.frame(read.table("./sep_est/P8717_105_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_106 <- as.data.frame(read.table("./sep_est/P8717_106_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_107 <- as.data.frame(read.table("./sep_est/P8717_107_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_108 <- as.data.frame(read.table("./sep_est/P8717_108_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_109 <- as.data.frame(read.table("./sep_est/P8717_109_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))


sample_110 <- as.data.frame(read.table("./sep_est/P8717_110_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_111 <- as.data.frame(read.table("./sep_est/P8717_111_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_112 <- as.data.frame(read.table("./sep_est/P8717_112_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_113 <- as.data.frame(read.table("./sep_est/P8717_113_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_114 <- as.data.frame(read.table("./sep_est/P8717_114_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_115 <- as.data.frame(read.table("./sep_est/P8717_115_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_116 <- as.data.frame(read.table("./sep_est/P8717_116_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_117 <- as.data.frame(read.table("./sep_est/P8717_117_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_118 <- as.data.frame(read.table("./sep_est/P8717_118_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_119 <- as.data.frame(read.table("./sep_est/P8717_119_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))


sample_120 <- as.data.frame(read.table("./sep_est/P8717_120_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_121 <- as.data.frame(read.table("./sep_est/P8717_121_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_122 <- as.data.frame(read.table("./sep_est/P8717_122_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_123 <- as.data.frame(read.table("./sep_est/P8717_123_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))

sample_124 <- as.data.frame(read.table("./sep_est/P8717_124_codeGene_est.csv",
                                       sep = ",",
                                       row.names=1,
                                       header =T,
                                       stringsAsFactors = F))


gene_ID <- codeGene_name$Row.names
gene_symbol <- codeGene_name$gene_symbol
P8717_101_est <- round(sample_101$est)
P8717_102_est <- round(sample_102$est)
P8717_103_est <- round(sample_103$est)
P8717_104_est <- round(sample_104$est)
P8717_105_est <- round(sample_105$est)
P8717_106_est <- round(sample_106$est)
P8717_107_est <- round(sample_107$est)
P8717_108_est <- round(sample_108$est)
P8717_109_est <- round(sample_109$est)
P8717_110_est <- round(sample_110$est)
P8717_111_est <- round(sample_111$est)
P8717_112_est <- round(sample_112$est)
P8717_113_est <- round(sample_113$est)
P8717_114_est <- round(sample_114$est)
P8717_115_est <- round(sample_115$est)
P8717_116_est <- round(sample_116$est)
P8717_117_est <- round(sample_117$est)
P8717_118_est <- round(sample_118$est)
P8717_119_est <- round(sample_119$est)
P8717_120_est <- round(sample_120$est)
P8717_121_est <- round(sample_121$est)
P8717_122_est <- round(sample_122$est)
P8717_123_est <- round(sample_123$est)
P8717_124_est <- round(sample_124$est)



# Combine all the data together in one dataframe and save
# The new data frame containing ensembol gene ID, gene names and rounding gene counts for 24 samples.

com_codeGene_est <- data.frame(gene_ID, gene_symbol,
                               P8717_101_est,P8717_102_est,P8717_103_est,P8717_104_est,P8717_105_est,P8717_106_est,P8717_107_est,P8717_108_est,P8717_109_est,P8717_110_est,
                               P8717_111_est,P8717_112_est,P8717_113_est,P8717_114_est,P8717_115_est,P8717_116_est,P8717_117_est,P8717_118_est,P8717_119_est,P8717_120_est,
                               P8717_121_est,P8717_122_est,P8717_123_est,P8717_124_est)


write.csv(com_codeGene_est, file = "./sep_est/com_codeGene_est.csv", quote = FALSE, row.names = FALSE)


# * Combine PKLR-siRNA groups data and save
# The data frame containing ensembol gene ID, gene names and rounding gene counts for PKLR-siRNA experimental groups as well as control groups.


PKLR.siRNA_EX_1 <- P8717_101_est
PKLR.siRNA_EX_2 <- P8717_102_est
PKLR.siRNA_EX_3 <- P8717_103_est
siRNA_NC_1 <- P8717_110_est
siRNA_NC_2 <- P8717_111_est
siRNA_NC_3 <- P8717_112_est

PKLR_si_est <- data.frame(gene_ID,
                          gene_symbol,
                          PKLR.siRNA_EX_1,
                          PKLR.siRNA_EX_2,
                          PKLR.siRNA_EX_3,
                          siRNA_NC_1,
                          siRNA_NC_2,
                          siRNA_NC_3)



# * Combine PKLR-plasmid groups data and save
# The data frame containing ensembol gene ID, gene names and rounding gene counts for PKLR-plasmid experimental groups as well as control groups.


PKLR.plasmid_EX_1 <- P8717_113_est
PKLR.plasmid_EX_2 <- P8717_114_est
PKLR.plasmid_EX_3 <- P8717_115_est
plasmid_NC_1 <- P8717_122_est
plasmid_NC_2 <- P8717_123_est
plasmid_NC_3 <- P8717_124_est

PKLR_ov_est <- data.frame(gene_ID,
                          gene_symbol,
                          PKLR.plasmid_EX_1,
                          PKLR.plasmid_EX_2,
                          PKLR.plasmid_EX_3,
                          plasmid_NC_1,
                          plasmid_NC_2,
                          plasmid_NC_3)

write.csv(PKLR_si_est, file = "PKLR_si_est_sym.csv",quote = FALSE, row.names = FALSE)
write.csv(PKLR_ov_est, file = "PKLR_ov_est_sym.csv",quote = FALSE, row.names = FALSE)


# * Simplized PKLR data group
# Compared with the last data, this one does not have gene symbol. The main reason is to avoid time consuming.

PKLR_si_est <- data.frame(gene_ID,
                          PKLR.siRNA_EX_1,
                          PKLR.siRNA_EX_2,
                          PKLR.siRNA_EX_3,
                          siRNA_NC_1,
                          siRNA_NC_2,
                          siRNA_NC_3)
PKLR_ov_est <- data.frame(gene_ID,
                          PKLR.plasmid_EX_1,
                          PKLR.plasmid_EX_2,
                          PKLR.plasmid_EX_3,
                          plasmid_NC_1,
                          plasmid_NC_2,
                          plasmid_NC_3)
write.csv(PKLR_si_est, file = "PKLR_si_est.csv",quote = FALSE, row.names = FALSE)
write.csv(PKLR_ov_est, file = "PKLR_ov_est.csv",quote = FALSE, row.names = FALSE)



# Differential expressed analysis (overexpressed-PKLR plasmid groups)

library(SummarizedExperiment)  
library(DESeq2) 

ov_counts_sym = as.matrix(read.table("PKLR_ov_est_sym.csv",
                                     sep = ",",
                                     header = T,
                                     stringsAsFactors = F,
                                     row.names = 1
))
ov_counts = as.matrix(read.table("PKLR_ov_est.csv",
                                 sep = ",",
                                 header = T,
                                 stringsAsFactors = F,
                                 row.names = 1
))


grep = grepl(colnames(ov_counts_sym), pattern = "gene_symbol")


treatment <- sapply(strsplit(colnames(ov_counts),"_"),"[[",1)
group <-sapply(strsplit(colnames(ov_counts),"_"),"[[",2)
sample <- sapply(strsplit(colnames(ov_counts),"_"),"[[",3)

metadata <- as.data.frame(cbind(treatment, group, sample))


se_ov <- SummarizedExperiment(assays = ov_counts, colData = metadata)
se_ov <- SummarizedExperiment(list(counts = ov_counts), colData=metadata)
rowData(se_ov) $symbol <- ov_counts_sym[,grep]


# DESeq2 analysis
# Differential expressed analysis by groups


dds_ov <- DESeqDataSetFromMatrix(countData = assays(se_ov)$counts, colData = metadata,
                                 design = ~1)
dds_ov$group <- factor(paste0(dds_ov$group))
design(dds_ov)<- ~group

dds_ov <-DESeq(dds_ov)
resultsNames(dds_ov)


# Combine result data with gene symbol and save

res_ov_deseq<- results(dds_ov,contrast = c('group','EX','NC'))
gene_ID <- rownames(se_ov)
gene_symbol <- rowData(se_ov)$symbol
res_ov_gene <- data.frame(gene_ID, gene_symbol)
row.names(res_ov_gene) <- res_ov_gene$gene_ID

res_ov_counts <- merge(res_ov_gene,
                       assays(se_ov)$counts,
                       by = "row.names",
                       sort = FALSE)
row.names(res_ov_counts)<- res_ov_counts[,1]

res_ov <- merge(res_ov_counts,
                as.data.frame(res_ov_deseq),
                by = "row.names",
                sort = FALSE)
res_ov <- res_ov [,-1]

write.csv(res_ov, file = "res_PKLR_ov.csv", quote = FALSE, row.names = FALSE)



# Differential expressed analysis (PKLR-siRNA groups)
# Almost the same procedure with plasmid group


si_counts_sym = as.matrix(read.table("PKLR_si_est_sym.csv",
                                     sep = ",",
                                     header = T,
                                     stringsAsFactors = F,
                                     row.names = 1
))
si_counts = as.matrix(read.table("PKLR_si_est.csv",
                                 sep = ",",
                                 header = T,
                                 stringsAsFactors = F,
                                 row.names = 1
))

grep = grepl(colnames(si_counts_sym), pattern = "gene_symbol")


treatment <- sapply(strsplit(colnames(si_counts),"_"),"[[",1)
group <-sapply(strsplit(colnames(si_counts),"_"),"[[",2)
sample <- sapply(strsplit(colnames(si_counts),"_"),"[[",3)

metadata <- as.data.frame(cbind(treatment, group, sample))


se_si <- SummarizedExperiment(assays = si_counts, colData = metadata)
se_si <- SummarizedExperiment(list(counts = si_counts), colData=metadata)
rowData(se_si) $symbol <- si_counts_sym[,grep]

dds_si <- DESeqDataSetFromMatrix(countData = assays(se_si)$counts, colData = metadata,
                                 design = ~1)
dds_si$group <- factor(paste0(dds_si$group))
design(dds_si)<- ~group

dds_si <-DESeq(dds_si)
resultsNames(dds_si)

res_si_deseq<- results(dds_si,contrast = c('group','EX','NC'))
gene_ID <- rownames(se_si)
gene_symbol <- rowData(se_si)$symbol
res_si_gene <- data.frame(gene_ID, gene_symbol)
row.names(res_si_gene) <- res_si_gene$gene_ID

res_si_counts <- merge(res_si_gene,
                       assays(se_si)$counts,
                       by = "row.names",
                       sort = FALSE)
res_si_counts <- res_si_counts[,-1]
row.names(res_si_counts)<- res_si_counts[,1]

res_si <- merge(res_si_counts,
                as.data.frame(res_si_deseq),
                by = "row.names",
                sort = FALSE)
res_si <- res_si [,-1]

write.csv(res_si, file = "res_PKLR_si.csv", quote = FALSE, row.names = FALSE)



# Detect the relationships between over-expressed group and inhibition group

# With the two gorups of data, we can detect if there is any correlation between them. For example, if they have a correlation coefficient which absolute value larger than 0.7, or if we can find a clearly linage relationship between them.

library (ggplot2)


deseq_siRNA <- as.data.frame(read.table('res_PKLR_si.csv',
                                        sep = ',',
                                        stringsAsFactors = F,
                                        row.names = 1,
                                        header = T))

deseq_plasmid <- as.data.frame(read.table('res_PKLR_ov.csv',
                                          sep = ',',
                                          stringsAsFactors = F,
                                          row.names = 1,
                                          header = T))


# Seperate the input log2foldchange data and save

cor_login <- data.frame(gene_ID = row.names(deseq_plasmid),
                        gene_symbol = deseq_plasmid$gene_symbol,
                        inhibit_log = deseq_siRNA$log2FoldChange,
                        ov_log = deseq_plasmid$log2FoldChange)
write.csv(cor_login, file="res_PKLR_merge_log.csv", row.names = F)





cor_login_data <- cor_login[,-2]
cor_login_data <- cor_login_data[,-1]
row.names(cor_login_data) <- cor_login$gene_ID

# correlation detection
#The data with value "NA" were ignored, and only the genes which have value in both groups were calculated. Spearman method was used here, which is more powerful for two different groups of data.


cor_logout <- cor(cor_login_data, use ="pairwise.complete.obs", method = c("spearman"))
print(cor_logout)

# Draw the scatter plot between overexpressed group and inhibition group and compare



my_log <- data.frame(my_ov_log,
                      my_si_log)
ggplot(my_log, aes(x = my_ov_log, y = my_si_log))+
  geom_point()


