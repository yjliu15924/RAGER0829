library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genekitr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
geneIDs <- readLines(file)
expData <- read.csv(args[2],row.names = 1)
exprIDs <- rownames(expData)
exprIDs <- unlist(lapply(strsplit(exprIDs,"\\|"),function(x) x[[1]]))
exprIDs <- gsub("\\.\\d+","",exprIDs)
rownames(expData) <- exprIDs
matchIndexes <- match(geneIDs,rownames(expData))
expData <- expData[matchIndexes,]
expData <- expData[, -1]
entrezIDs <- transId(
   id = rownames(expData),
   transTo = "entrez", org = "human", keepNA = FALSE
)

entrezIDs <- entrezIDs[!duplicated(entrezIDs$input_id), ]

expData$input_id <- rownames(expData)

mergedData <- merge(entrezIDs, expData, by = "input_id")
sortedData <- mergedData[order(mergedData$log2FoldChange, decreasing = TRUE), ]
geneList <- sortedData$log2FoldChange
names(geneList) <- sortedData$entrezid
res <- gseGO(
   geneList,
   ont = "ALL",
   OrgDb = org.Hs.eg.db,
   keyType = "ENTREZID",
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "pos",
)
gsea_result = res@result
write.csv(gsea_result,args[3])

res1 <- gseKEGG(
   geneList = geneList,
   organism = "hsa",  # Human KEGG code
   keyType = "ncbi-geneid",  # If using EntrezID, otherwise adjust as needed
   pvalueCutoff = 1,
   pAdjustMethod = "fdr",
   scoreType = "pos"
)
gsea_result1 = res1@result
write.csv(gsea_result1,args[4])

# #这一部分作为单独的，不在snakemake中
# # 指定需要绘制的通路ID
# target_pathway_id <- 'GO:0070372'
# 
# # 查找该 ID 所对应的行号
# geneSetIDs <- which(gsea_result$ID == target_pathway_id)
# 
# # 设置画图输出
# pdf("/home/yjliu/RAGER/case2/datasets//GSEA_plot_for_GO:0070372.pdf", width = 10, height = 8)
# 
# # 绘制 GSEA 图
# gseaplot2(res, 
#           geneSetID = geneSetIDs, 
#           title = gsea_result$Description[geneSetIDs])
# 
# # 关闭绘图设备
# dev.off()