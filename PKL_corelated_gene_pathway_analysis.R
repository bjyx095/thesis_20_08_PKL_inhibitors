# pathway analysis
# functional analysis

library(clusterProfiler)
library(org.Hs.eg.db)

df_3 <- read.csv("PKL_statistical_overlap_genes_strict.csv",
                 header=T,
                 stringsAsFactors = F)
df_3 <- subset(df_3, df_3[,"sum_z"]>=0)

gene_list_3 <- df_3[,"sum_z"]

names(gene_list_3) <- as.character(df_3[,"gene_num"])
gene_list_3 <- sort(gene_list_3, decreasing=TRUE)

gene_3 <- as.data.frame(gene_3)
ptw_gsea_kegg <- gseKEGG(
  geneList = gene_list_3,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  by="fgsea"
)
ptw_gsea_kegg_res <- ptw_gsea_kegg@result

tiff(filename ="ptw_gsea_kegg_graph_up.tiff", width=720, height=720)
dotplot(ptw_gsea_kegg, color="pvalue",  showCategory = 40, font.size=15)
dev.off()