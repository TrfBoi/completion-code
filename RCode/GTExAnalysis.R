#
### V8数据处理
# 
options(stringsAsFactors = F)
GTExExpr <- read.table('D:/Users/Lenovo/Desktop/毕设/code/RCode/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct'
                , header = T, sep = '\t', skip = 2)
GTEXSimpleInfo <- read.table("D:\\Users\\Lenovo\\Desktop\\毕设\\GTEx\\GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
                             , header = T, sep = '\t')
GTEXSimpleInfo <- GTEXSimpleInfo[, c("SAMPID", "SMTS", "SMTSD")]
GTEXSimpleInfo <- GTEXSimpleInfo[GTEXSimpleInfo$SMTSD == "Whole Blood", ]
GTExExpr <- GTExExpr[, c(T, T, gsub('[.]','-', colnames(GTExExpr)[3:ncol(GTExExpr)]) %in% GTEXSimpleInfo$SAMPID)]
rm(list = c("GTEXSimpleInfo", "GTExSimpleIds"))
write.table(GTExExpr, file = "GTEx/GTExExpr.txt", sep = "\t", row.names = F)
#
### V7数据处理
#
options(stringsAsFactors = F)
GTExExpr <- read.table("D:\\Users\\Lenovo\\Desktop\\毕设\\GTEx\\v7\\All_Tissue_Site_Details.combined.reads.gct"
                       , header = T, sep = '\t', skip = 2)
GTEXSimpleInfo <- read.table("D:\\Users\\Lenovo\\Desktop\\毕设\\GTEx\\v7\\GTEx_v7_Annotations_SampleAttributesDS.txt"
                             , header = T, sep = '\t')
GTEXSimpleInfo <- GTEXSimpleInfo[, c("SAMPID", "SMTS", "SMTSD")]
GTEXSimpleInfo <- GTEXSimpleInfo[GTEXSimpleInfo$SMTS == "Blood", ]
GTExExpr <- GTExExpr[, c(T, T, gsub('[.]','-', colnames(GTExExpr)[3:ncol(GTExExpr)]) %in% GTEXSimpleInfo$SAMPID)]
rm(list = c("GTEXSimpleInfo", "GTExSimpleIds"))
write.table(GTExExpr, file = "GTEx/GTExExpr_v7.txt", sep = "\t", row.names = F)
#
### 表达矩阵处理
#
GTExExpr <- read.table(file = "GTEx/GTExExpr.txt", header = T, sep = '\t')
for (i in 1:length(strsplit(GTExExpr$Name, split = "\\."))) {
    GTExExpr[i, "Name"] <- tmpENSGIdList[[i]][1]
}
library(org.Hs.eg.db)
tmp <- mapIds(org.Hs.eg.db, keys = GTExExpr[, 1], keytype = "ENSEMBL", column="SYMBOL")
GTExExpr[, 1] <- tmp
GTExExpr <- GTExExpr[, -2]
GTExExpr <- aggregate(x = GTExExpr[, 2:ncol(GTExExpr)], by = list(GTExExpr[, 1]), FUN = mean)
rownames(GTExExpr) <- GTExExpr[,1]
GTExExpr <- GTExExpr[,-1]
GTExExpr <- GTExExpr[rowSums(GTExExpr) > 500 & apply(GTExExpr, 1, function(x){ all(x > 0) }),]
write.table(GTExExpr, file = "GTEx/GTExExprFilter.txt", sep = "\t")


##########################################################################################################################

#
### 与TCGA批次效应处理及差异分析
#
GTExExpr <- read.table(file = "GTEx/GTExExprFilter.txt", sep = "\t", row.names = 1, header = T)
sameGene <- intersect(row.names(GTExExpr), row.names(filterGeneExpressMatrix))
GTExExpr <- GTExExpr[sameGene, ]
filterGeneExpressMatrix <- filterGeneExpressMatrix[sameGene, ]
rm(sameGene)
degExp <- cbind(GTExExpr, filterGeneExpressMatrix)
# rm(GTExExpr)

library(edgeR)
group_info = c(rep("ctrl", 51), rep("KD", 151))
edgeRObj <- DGEList(counts = degExp, group = group_info)
edgeRObj <- calcNormFactors(edgeRObj, method = "TMM")
plotMDS(edgeRObj, col = c(rep("red", 51), rep("blue", 151)), pch = 16)
design.mat <- model.matrix(~group_info)
# estimate dispersion
edgeRObj <- estimateDisp(edgeRObj, design.mat)
edgeRObj$common.dispersion
edgeRObj$tagwise.dispersion
# 1st common dispersion
edgeRObj <- estimateCommonDisp(edgeRObj)
# 2nd tagwise dispersion
edgeRObj <- estimateTagwiseDisp(edgeRObj)
# plotBCV(edgeRObj, cex = 0.8)
# plot var and mean
# plotMeanVar(edgeRObj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
# test with likelihood ratio test
dgeRes1 <- exactTest(edgeRObj)
dgeRes <- as.data.frame(topTags(dgeRes1, n=nrow(degExp), sort.by = "logFC"))
# MA plot
# select.sign.gene = decideTestsDGE(dgeRes1, p.value = 0.001) 
# select.sign.gene_id = rownames(dgeRes1)[as.logical(select.sign.gene)]
# plotSmear(dgeRes1, de.tags = select.sign.gene_id, cex = 0.5, ylim=c(-4,4)) 
# abline(h = c(-2, 2), col = "blue")

rm(list = c("group_info", "design.mat", "select.sign.gene", "select.sign.gene_id"))
rm(list = c("edgeRObj", "dgeRes1"))
library(ggplot2)
huoshanRes <- dgeRes[rownames(dgeRes) %in% cellAgeInfo[, "gene_name"], ]
huoshanRes$Sig = ifelse(huoshanRes$PValue < 0.05 & 
                        abs(huoshanRes$logFC) >= 0.585, 
                    ifelse(huoshanRes$logFC > 0.585 ,'Up','Down'),'None')
table(huoshanRes$Sig)
ggplot(huoshanRes, aes(x = logFC, y = -log10(PValue), colour=Sig)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
    # 辅助线
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),
               lty=4,col="black",lwd=0.8) +
    # 坐标轴
    labs(x="log2(Fold Change)",
         y="-log10 (P-value)")+
    theme_bw()+
    ggtitle("Volcano Plot")+
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank()
    )
dgeRes <- dgeRes[dgeRes$PValue < 0.05, ]
# 差异1.5倍, 设置为1就是差异2倍
dgeRes <- dgeRes[abs(dgeRes$logFC) >= 0.585, ]
degExp <- degExp[rownames(dgeRes), ]
diffCellAgeGene <- rownames(degExp[rownames(degExp) %in% cellAgeInfo[, "gene_name"], ])
# 热图，这里只绘制一下细胞衰老基因的热图(差异度前8，在所有基因中前500)
# heatMapIn <- degExp[rownames(degExp) %in% cellAgeInfo[, "gene_name"], ][1:8, ]
# annotation <- data.frame(type=c(rep("normal", 51), rep("tumor", 151)))
# rownames(annotation) <- colnames(heatMapIn)
# library(pheatmap)
# pheatmap(as.matrix(heatMapIn), 
#          annotation = annotation, 
#          cluster_cols = FALSE,
#          fontsize_row = 11,
#          show_colnames = F,
#          fontsize_col = 3,
#          color = colorRampPalette(c("green", "black", "red"))(50) )

# 火山图
rm(list = c("annotation", "heatMapIn", "degExp", "dgeRes", "huoshanRes"))
filterGeneExpressMatrix <- filterGeneExpressMatrix[diffCellAgeGene, ]
write.table(filterGeneExpressMatrix, file = "GTEx/diffCellAgeGene.txt", sep = "\t", row.names = T, col.names = T)



##### GO与KEGG分析
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05           
qvalueFilter=0.05 

colorSel="qvalue"
if(qvalueFilter>0.05){
    colorSel="pvalue"
}
# 基因向量
genes=diffCellAgeGene
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"] 

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter), ]

write.table(GO, file="GTEx/GO.txt", sep="\t", quote=F, row.names = F)

showNum = 10
if (nrow(GO)<30) {
    showNum=nrow(GO)
}

pdf(file="GTEx/go_barplot.pdf",width = 10,height = 10)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

pdf(file="GTEx/go_bubble.pdf",width = 10,height = 10)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###
### kegg
###
pvalueFilter=0.05          
qvalueFilter=0.05  

colorSel="qvalue"
if(qvalueFilter>0.05){
    colorSel="pvalue"
}

genes=diffCellAgeGene
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="GTEx/KEGG.txt",sep="\t",quote=F,row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
    showNum=nrow(KEGG)
}

pdf(file="GTEx/kegg_barplot.pdf",width = 9,height = 9)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file="GTEx/kegg_bubble.pdf",width = 9,height = 9)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

