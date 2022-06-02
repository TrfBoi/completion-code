#
### 读取RNA-Seq测序数据
#
geneExpressMatrix <- read.table("D:/Users/Lenovo/Desktop/毕设/matrix.txt", sep=" ", row.names = 1, encoding = "UTF-8", header = TRUE)
# Ensembl gene ID转成gene symbol。同名基因合并, 表达量取平均值
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
geneExpressMatrix <- geneExpressMatrix[1:(length(rownames(geneExpressMatrix))-5),]
geneExpressMatrix <- cbind(rownames(geneExpressMatrix), geneExpressMatrix)
geneExpressMatrix[,1] <- mapIds(org.Hs.eg.db, keys = geneExpressMatrix[,1], keytype = "ENSEMBL", column="SYMBOL")
length(unique(geneExpressMatrix[,1]))
uniqueGeneExpressMatrix <- aggregate(x = geneExpressMatrix[, 2:ncol(geneExpressMatrix)], by = list(geneExpressMatrix[, 1]), FUN = mean)
rownames(uniqueGeneExpressMatrix) <- uniqueGeneExpressMatrix[,1]
uniqueGeneExpressMatrix <- uniqueGeneExpressMatrix[,-1]
filterGeneExpressMatrix <- uniqueGeneExpressMatrix[rowSums(uniqueGeneExpressMatrix) > 500 & apply(uniqueGeneExpressMatrix, 1, function(x){ all(x > 0) }),]
#
# 从这里去跑差异分析和批次校验
#
allGeneExp <- filterGeneExpressMatrix

################################
# 筛选出细胞衰老基因表达矩阵   #
################################
cellAgeInfo <- read.table("D:/Users/Lenovo/Desktop/毕设/cellAge/cellAge1.csv", sep = ";", row.names = 1, header = T, encoding = "UTF-8")
cellAgeInfo <- cellAgeInfo[cellAgeInfo[, "organism"] == "Human",]
inhibitsCellAgeGene <- cellAgeInfo[cellAgeInfo[,"senescence_effect"] == "Inhibits", ][, "gene_name"]
inducesCellAgeGene <- cellAgeInfo[cellAgeInfo[,"senescence_effect"] == "Induces", ][, "gene_name"]
unClearCellAgeGene <- cellAgeInfo[cellAgeInfo[,"senescence_effect"] == "Unclear", ][, "gene_name"]
###
### 到GTEx中找出差异细胞衰老基因
###
filterGeneExpressMatrix <- filterGeneExpressMatrix[rownames(filterGeneExpressMatrix) %in% cellAgeInfo[, "gene_name"], ]
# diffCellAgeGene <- rownames(read.table(file = "GTEx/diffCellAgeGene.txt", sep = "\t", row.names = 1, header = T))
# diffCellAgeGene <- diffCellAgeGene[diffCellAgeGene %in% rownames(geoExprMatrix)]
filterGeneExpressMatrix <- filterGeneExpressMatrix[rownames(filterGeneExpressMatrix) %in% rownames(geoExprMatrix), ]
rm(list = c("inducesCellAgeGene", "inhibitsCellAgeGene", "unClearCellAgeGene", "geneExpressMatrix", "uniqueGeneExpressMatrix"))
###
### 到GEO数据库脚本做批次校验，但是做不了，不知道GEO数据是啥格式
###

###################################
# 整理出生存分析所需input matrix  #
###################################
survivalInputMatrix <- read.table("D:/Users/Lenovo/Desktop/毕设/survivalInputMatrix.txt", sep = "\t", encoding = "UTF-8")
survivalInputMatrix <- t(survivalInputMatrix)
colnames(survivalInputMatrix) <- survivalInputMatrix[1, ]
survivalInputMatrix <- survivalInputMatrix[-1, ]
survivalInputMatrix[, "age_at_index"] <- as.integer(survivalInputMatrix[, "age_at_index"])
survivalInputMatrix[, "alive"] <- as.integer(survivalInputMatrix[, "alive"])
survivalInputMatrix[, "time"] <- as.integer(survivalInputMatrix[, "time"])
survivalInputMatrix <- survivalInputMatrix[survivalInputMatrix[, "time"] > -1, ] # 去掉缺失生存时间的样本
survivalInputMatrix <- as.data.frame(survivalInputMatrix)
############################
# 构建COX模型所需输入矩阵  #
############################
coxMatrix <- survivalInputMatrix[, c("id", "time", "alive")]
coxMatrix[, "id"] <- gsub("-", ".", coxMatrix[, "id"])
colNames <- strsplit(colnames(filterGeneExpressMatrix), split = "\\.")
tmpIds <- rep(" ", ncol(filterGeneExpressMatrix))
for (i in 1:length(colNames)) {
    tmpIds[i] <- paste(colNames[[i]][1:3], collapse = ".")
}
colnames(filterGeneExpressMatrix) <- tmpIds
colNames <- strsplit(colnames(allGeneExp), split = "\\.")
tmpIds <- rep(" ", ncol(allGeneExp))
for (i in 1:length(colNames)) {
    tmpIds[i] <- paste(colNames[[i]][1:3], collapse = ".")
}
colnames(allGeneExp) <- tmpIds
tmpGeneExp <- filterGeneExpressMatrix[, colnames(filterGeneExpressMatrix) %in% coxMatrix[, "id"]]
rm(colNames)
rm(tmpIds)
tmpGeneExp <- t(tmpGeneExp)
tmpGeneExp <- cbind(id=rownames(tmpGeneExp), tmpGeneExp)
coxMatrix <- merge(coxMatrix, tmpGeneExp, by = "id")
rm(tmpGeneExp)



###############################################
# 单因素COX分析得到预后相关的细胞衰老相关基因 #
###############################################
# install.packages('survival')
library(survival)
for (i in 2:ncol(coxMatrix)) {
    coxMatrix[, i] <- as.integer(coxMatrix[, i])
}
coxMatrix[, 4:ncol(coxMatrix)] <- log2(coxMatrix[, 4:ncol(coxMatrix)] + 1) # 取log2

coxOut <- data.frame()
importantGenes <- colnames(coxMatrix)[2:3]
for (i in colnames(coxMatrix[, 4:ncol(coxMatrix)])) {
    cox <- coxph(Surv(time, alive) ~ coxMatrix[, i], data = coxMatrix) # , ties = "exact"
    coxSummary <- summary(cox)
    coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
    if (coxP < 0.05) {
        importantGenes <- c(importantGenes, i)
        coxOut <- rbind(coxOut,
                     cbind(id = i, HR = coxSummary$conf.int[, "exp(coef)"],
                           HR.95L = coxSummary$conf.int[, "lower .95"],
                           HR.95H = coxSummary$conf.int[, "upper .95"],
                           pvalue = coxSummary$coefficients[, "Pr(>|z|)"]) )
    }
}
rm(coxP)
rm(coxSummary)
rm(cox)
importantGenesExp <- coxMatrix[, importantGenes]
importantGenesExp <- cbind(id = coxMatrix[, 1], importantGenesExp)
rownames(importantGenesExp) <- importantGenesExp[, 1]
importantGenesExp <- importantGenesExp[, -1]
for (i in 2:ncol(coxOut)) {
    coxOut[, i] <- as.double(coxOut[, i])
}
coxOut
coxOut <- coxOut[coxOut$pvalue < 0.02, ] # 由于单因素基因结果过多，绘图时选取较为明显的几个结果 即p < 0.02的18个基因
gene <- coxOut[, 1]
hr <- sprintf("%.3f", coxOut$"HR")
hrLow  <- sprintf("%.3f", coxOut$"HR.95L")
hrHigh <- sprintf("%.3f", coxOut$"HR.95H")
Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
pVal <- ifelse(coxOut$pvalue < 0.001, "<0.001", sprintf("%.3f", coxOut$pvalue))
# pdf(file="singleFactorForest.pdf", width = 8, height = 4.5)
n <- nrow(coxOut)
nRow <- n + 1
ylim <- c(1, nRow)
layout(matrix(c(1,2), nc=2), width=c(2,2))
xlim = c(0, 4)
par(mar=c(2, 1, 1, 0))
plot(1, xlim = xlim, ylim = ylim, type="n", axes=F, xlab="", ylab="")
text.cex <- 0.8
text(0, n:1, gene,adj = 0, cex = text.cex)
text(1.5 - 0.5*0.2, n:1, pVal, adj = 1, cex = text.cex)
text(1.5 - 0.5*0.2, n+1,'pvalue', cex = text.cex, font = 2, adj = 1)
text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(3, n+1, 'Hazard ratio', cex=text.cex, font=2, adj=1)
par(mar = c(4, 2, 0, 1), mgp = c(2, 0.5, 0))
xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)) + 1)
plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, ylab="", xaxs="i", xlab="Hazard ratio")
arrows(as.numeric(hrLow), n:1,as.numeric(hrHigh), n:1, angle=90, code=3, length=0.05, col="darkblue", lwd=2.5)
abline(v=1, col="black", lty=2, lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
# dev.off()
rm(list=c("gene", "hr", "hrHigh", "hrLow", "Hazard.ratio", "pVal", "n", "nRow", "ylim", "xlim", "boxcolor", "text.cex"))



##################
# Go 与 KEGG分析 #
##################
par()
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
genes=importantGenes[3:length(importantGenes)]
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"] 

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter), ]

write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

showNum = 10
if (nrow(GO)<30) {
    showNum=nrow(GO)
}

pdf(file="go_barplot.pdf",width = 10,height = 7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

pdf(file="go_bubble.pdf",width = 10,height = 7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


pvalueFilter=0.05          
qvalueFilter=0.05  

colorSel="qvalue"
if(qvalueFilter>0.05){
    colorSel="pvalue"
}

genes=importantGenes[3:length(importantGenes)]
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
    showNum=nrow(KEGG)
}

pdf(file="kegg_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file="kegg_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()


##########
# 分型   #
##########
# install.packages("survival")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ConsensusClusterPlus")
# library(limma)
# library(survival)
# library(ConsensusClusterPlus)
# maxK <- 9
# clusterIn <- importantGenesExp[, c(-1, -2)]
# clusterIn <- t(as.matrix(clusterIn))
# clusterResults <- ConsensusClusterPlus(clusterIn,
#                              maxK = maxK,
#                              reps = 50,
#                              pItem = 0.8,
#                              pFeature = 1,
#                              clusterAlg = "km",
#                              title = "clusterOut",
#                              distance = "euclidean",
#                              seed = 123456,
#                              plot = "png")
# 
# clusterNum <- 2
# clusterOut <- clusterResults[[clusterNum]][["consensusClass"]]
# clusterOut <- as.data.frame(clusterOut)
# clusterOut[, 1] <- paste0("C", clusterOut[, 1])
# clusterOut <- cbind(id = rownames(clusterOut), clusterOut)
# rm(list = c("maxK", "clusterResults", "clusterNum", "clusterIn"))

################
# 分型生存分析 #
################
# install.packages("survival")
# install.packages("survminer")
# library(survival)
# library(survminer)
# clusterSurvivalIn <- coxMatrix[, 1:3]
# clusterSurvivalIn$time <- clusterSurvivalIn$time/365
# clusterSurvivalIn <- merge(clusterSurvivalIn, clusterOut, by = "id")
# rownames(clusterSurvivalIn) <- clusterSurvivalIn[, 1]
# colnames(clusterSurvivalIn) <- c(colnames(clusterSurvivalIn)[-ncol(clusterSurvivalIn)], "clusterType")
# clusterSurvivalIn <- clusterSurvivalIn[, -1]
# length <- length(levels(factor(clusterSurvivalIn$clusterType)))
# diff <- survdiff(Surv(time, alive) ~ clusterType, data = clusterSurvivalIn)
# pValue <- 1-pchisq(diff$chisq, df=length-1)
# if (pValue<0.001) {
#     pValue="p<0.001"
# } else {
#     pValue=paste0("p=",sprintf("%.03f",pValue))
# }
# fit <- survfit(Surv(time, alive) ~ clusterType, data = clusterSurvivalIn)
# bioCol <- c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol <- bioCol[1:length]
# surPlot <- ggsurvplot(fit, 
#                    data=clusterSurvivalIn,
#                    conf.int=F,
#                    pval=pValue,
#                    pval.size=6,
#                    legend.title="Cluster",
#                    legend.labs=levels(factor(clusterSurvivalIn[,"clusterType"])),
#                    legend = c(0.8, 0.8),
#                    font.legend=10,
#                    xlab="Time(years)",
#                    break.time.by = 1,
#                    palette = bioCol,
#                    surv.median.line = "hv",
#                    risk.table=T,
#                    cumevents=F,
#                    risk.table.height=.25)
# pdf(file="clusterSurvival.pdf", onefile = FALSE, width=7, height=5.5)
# print(surPlot)
# dev.off()
# rm(list = c("clusterSurvivalIn", "length", "diff", "pValue", "fit", "bioCol", "surPlot"))

################
# 分型差异分析 #
################
# BiocManager::install("edgeR")
# library(edgeR)
# degExp <- allGeneExp[, clusterOut$id]
# colnames(degExp) <- clusterOut[colnames(degExp), "clusterOut"]
# group_info = c(rep("ctrl", length(degExp[, colnames(degExp) == "C1"])), rep("KD", length(degExp[, colnames(degExp) == "C2"])))
# degExp <- cbind(degExp[, colnames(degExp) == "C1"], degExp[, colnames(degExp) == "C2"])
# 
# edgeRObj <- DGEList(counts = degExp, group = group_info)
# edgeRObj <- calcNormFactors(edgeRObj, method = "TMM")
# design.mat <- model.matrix(~group_info)
# # estimate dispersion
# edgeRObj <- estimateDisp(edgeRObj, design.mat)
# edgeRObj$common.dispersion
# edgeRObj$tagwise.dispersion
# # 1st common dispersion
# edgeRObj <- estimateCommonDisp(edgeRObj)
# # 2nd tagwise dispersion
# edgeRObj <- estimateTagwiseDisp(edgeRObj)
# # test with likelihood ratio test
# dgeRes <- exactTest(edgeRObj)
# dgeRes <- as.data.frame(topTags(dgeRes, n=nrow(degExp), sort.by = "logFC"))
# rm(list = c("group_info", "degExp", "design.mat"))
# dgeRes <- dgeRes[dgeRes$PValue < 0.05, ]
# # 差异1.5倍, 设置为1就是差异2倍
# dgeRes <- dgeRes[abs(dgeRes$logFC) >= 0.585, ]
# degExp <- allGeneExp[rownames(allGeneExp) %in% rownames(dgeRes), clusterOut$id] # 1273
# 
# # 热图，这里只绘制一下细胞衰老基因的热图
# degCellAgeExp <- degExp[rownames(degExp) %in% cellAgeInfo[, "gene_name"], ]
# rm(edgeRObj)


###################################
# TODO 预后相关基因GO与KEGG分析   #
###################################

#######################
# Lasso Cox 回归分析  #
#######################
# install.packages("glmnet")
# install.packages("survival")
library("glmnet")
library("survival")
lassoGeneExp <- importantGenesExp
lassoGeneExp$time[lassoGeneExp$time <= 0] = 1 # lasso不允许参数为0

x <- as.matrix(lassoGeneExp[, c(3:ncol(lassoGeneExp))])
y <- data.matrix(Surv(lassoGeneExp$time, lassoGeneExp$alive))

fit <- glmnet(x, y, family = "cox", alpha = 0.7)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox")
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene = row.names(coef)[index]  # 14个基因
lassoGene = c("time","alive",lassoGene)
lassoSigExp = lassoGeneExp[, lassoGene]
lassoSigExp = cbind(id=row.names(lassoSigExp),lassoSigExp)
rm(list  = c("x", "y", "fit", "actCoef", "coef", "index", "cvfit", "lassoGeneExp"))


#############################################
# 多因素COX分析构建预后模型 TODO 输入是lass #
#############################################
lassoSigExp$time <- lassoSigExp$time / 365
lassoSigExp <- lassoSigExp[, -1]

multiCox <- coxph(Surv(time, alive) ~ ., data = lassoSigExp) # .代表单因素分析中的所有基因共同作用
# multiCox <- step(multiCox, direction = "both") # 减少基因数量
multiCoxSum <- summary(multiCox)
coxOut <- data.frame()
coxOut <- cbind(
    coef = multiCoxSum$coefficients[,"coef"],
    HR = multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L = multiCoxSum$conf.int[,"lower .95"],
    HR.95H = multiCoxSum$conf.int[,"upper .95"],
    pvalue = multiCoxSum$coefficients[,"Pr(>|z|)"])
coefOut <- coxOut
coxOut <- cbind(id = rownames(coxOut), coxOut)
riskScore <- predict(multiCox, type="risk", newdata=lassoSigExp)
outCol <- c("time", "alive", rownames(multiCoxSum$coefficients))
riskType <- as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
riskOut <- cbind(id = rownames(lassoSigExp), cbind(lassoSigExp[, outCol], riskScore, riskType))
rm(list=c("riskScore", "outCol", "riskType"))

################
# 绘制生存曲线 #
################
# install.packages("survminer")
library(survminer)
diff <- survdiff(Surv(time, alive) ~ riskType, data = riskOut)
pValue <- 1-pchisq(diff$chisq, df = 1)
pValue <- signif(pValue,4)
pValue <- format(pValue, scientific = TRUE)
fit <- survfit(Surv(time, alive) ~ riskType, data = riskOut)
pdf(file="survivorshipCurve.pdf", onefile = FALSE, width = 5.5, height =5)
ggsurvplot(fit, 
           data=riskOut,
           conf.int=TRUE,
           pval=paste0("p=", pValue),
           pval.size=4,
           risk.table=TRUE, # 是否绘制每一阶段生存人数
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()
summary(fit) # 查看生存率
rm(list = c("pValue", "diff", "fit"))


# TODO 高低风险差异分析，绘制差异热图，


#################################
# 绘制风险曲线（风险得分相关图）#
#################################
# install.packages("pheatmap")
par()
library(pheatmap)
riskOut <- riskOut[order(riskOut$riskScore), ]
riskClass <- riskOut[, "riskType"]
lowLength <- length(riskClass[riskClass == "low"])
highLength <- length(riskClass[riskClass == "high"])
line <- riskOut[, "riskScore"]
line[line > 10] <- 10 # 降低x轴长度
pdf(file="riskScore.pdf", width = 10, height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green", lowLength),
           rep("red", highLength)))
abline(h = median(riskOut$riskScore), v=lowLength, lty=2)
legend("topleft", c("High risk", "low Risk"), bty="n", pch=19, col=c("red","green"), cex=1.2)
rm(list = c("riskClass", "highLength", "line"))
dev.off()
color <- as.vector(riskOut$alive)
color[color == 1] <- "red"
color[color == 0] <- "green"
pdf(file="survStat.pdf", width = 10, height = 4)
plot(riskOut$time,
     pch = 19,
     xlab = "Patients (increasing risk socre)",
     ylab = "Survival time (years)",
     col = color)
legend("topright", c("Dead", "Alive"), bty="n", pch=19, col=c("red","green"), cex=1.2)
abline(v=lowLength, lty=2)
rm(list = c("color", "lowLength"))
dev.off()
heatMapIn <- riskOut[c(4:(ncol(riskOut)-2))]
heatMapIn <- t(heatMapIn)
annotation <- data.frame(type=riskOut[, ncol(riskOut)])
rownames(annotation) <- rownames(riskOut)
pdf(file="heatmap.pdf", width = 10, height = 4)
pheatmap(heatMapIn, 
         annotation = annotation, 
         cluster_cols = FALSE,
         fontsize_row = 11,
         show_colnames = F,
         fontsize_col = 3,
         color = colorRampPalette(c("green", "red"))(50) )
dev.off()
rm(list = c("heatMapIn", "annotation"))

############################################################
# 独立分析 单因素和多因素, 如果p值都较为客观, 则是影响因素 #
############################################################
library(survival)
survivalInputMatrix[, 1] <- gsub("-", ".", survivalInputMatrix[, 1])
oneAnalysisIn <- merge(riskOut[, c(1, ncol(riskOut)-1)], survivalInputMatrix, by = "id")
oneAnalysisIn[, "gender"] <- ifelse(oneAnalysisIn[, "gender"] == "male", 1, 0)
oneAnalysisIn <- oneAnalysisIn[, colnames(oneAnalysisIn) != "race"]
for (i in 2:ncol(oneAnalysisIn)) {
    oneAnalysisIn[, i] <- as.numeric(oneAnalysisIn[, i])
}
oneAnalysisIn <- cbind(oneAnalysisIn[, colnames(oneAnalysisIn) %in% c("id", "time", "alive")],
                       oneAnalysisIn[, !(colnames(oneAnalysisIn) %in% c("id", "time", "alive"))])
rownames(oneAnalysisIn) <- oneAnalysisIn$id
oneAnalysisIn <- oneAnalysisIn[, !(colnames(oneAnalysisIn) %in% c("id"))]
coxOut <- data.frame()
for(i in colnames(oneAnalysisIn[, 3:ncol(oneAnalysisIn)])){
    cox <- coxph(Surv(time, alive) ~ oneAnalysisIn[, i], data = oneAnalysisIn)
    coxSummary <- summary(cox)
    coxOut <- rbind(coxOut,
                 cbind(id = i,
                       HR = coxSummary$conf.int[,"exp(coef)"],
                       HR.95L = coxSummary$conf.int[,"lower .95"],
                       HR.95H = coxSummary$conf.int[,"upper .95"],
                       pvalue = coxSummary$coefficients[,"Pr(>|z|)"]))
}
multiCox <- coxph(Surv(time, alive) ~ ., data = oneAnalysisIn)
multiCoxSum <- summary(multiCox)
multiCoxOut <- data.frame()
multiCoxOut <- cbind(
    HR = multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L = multiCoxSum$conf.int[,"lower .95"],
    HR.95H = multiCoxSum$conf.int[,"upper .95"],
    pvalue = multiCoxSum$coefficients[,"Pr(>|z|)"])
multiCoxOut <- cbind(id = row.names(multiCoxOut), multiCoxOut)
bioForest <- function(rt = data.frame(), fileName) {
    rt <- as.data.frame(rt)
    gene <- rt[, "id"]
    for (i in colnames(rt)) {
        if (i == "id") {
            next
        }
        rt[, i] <- as.numeric(rt[, i])
    }
    hr <- sprintf("%.3f", rt$"HR")
    hrLow  <- sprintf("%.3f", rt$"HR.95L")
    hrHigh <- sprintf("%.3f", rt$"HR.95H")
    Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
    pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    n <- nrow(rt)
    nRow <- n + 1
    pdf(file = fileName, width = 10, height = 4.5)
    ylim <- c(1, nRow)
    layout(matrix(c(1,2), nc=2), width=c(2, 2))
    xlim = c(0,3)
    par(mar=c(4,2.5,2,1))
    plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
    text.cex=0.8
    text(0,n:1,gene,adj=0,cex=text.cex)
    text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
    text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
    par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
    xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
    arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
    abline(v=1,col="black",lty=2,lwd=2)
    boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
    axis(1)
    dev.off()
}
bioForest(rt = coxOut, "independentSingle.pdf")
bioForest(rt = multiCoxOut, "independentMulti.pdf")

############
# ROC 分析 #
############
# install.packages("survivalROC")
library(survivalROC)
oneAnalysisIn$time <- oneAnalysisIn$time / 365
rocCol <- rainbow(ncol(oneAnalysisIn) - 2)
aucText <- c()
par()
pdf("ROC single.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc <- survivalROC(Stime=oneAnalysisIn$time, status=oneAnalysisIn$alive, marker = oneAnalysisIn$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText, paste0("risk score"," (AUC=", sprintf("%.3f", roc$AUC), ")"))
abline(0,1)
j = 1
oneAnalysisIn <- as.data.frame(oneAnalysisIn)
for (i in colnames(oneAnalysisIn[, 4:ncol(oneAnalysisIn)])) {
    roc=survivalROC(Stime=oneAnalysisIn$time, status=oneAnalysisIn$alive, marker = oneAnalysisIn[, i], predict.time =1, method="KM")
    j = j+1
    aucText=c(aucText,paste0(i, " (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol)
dev.off()
rm(list = c("aucText", "rocCol", "j", "bioForest", "multiCoxOut", 
            "multiCoxSum", "multiCox", "cox", "coxSummary", "coxOut", "roc"))
# 时间ROC曲线
#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")
library(survival)
library(survminer)
library(timeROC)
timeRoc <- function(input=null, rocFile=null){
    ROC_rt <- timeROC(T=input$time, delta=input$alive,
                   marker=input$riskScore, cause=1,
                   weighting='aalen',
                   times=c(1,3,5), ROC=TRUE)
    pdf(file=rocFile, width=5, height=5)
    plot(ROC_rt,time=1, col='green', title=FALSE, lwd=2)
    plot(ROC_rt,time=3, col='blue', add=TRUE, title=FALSE, lwd=2)
    plot(ROC_rt,time=5, col='red', add=TRUE, title=FALSE, lwd=2)
    legend('bottomright',
           c(paste0('AUC at 1 years: ', sprintf("%.03f",ROC_rt$AUC[1])),
             paste0('AUC at 3 years: ', sprintf("%.03f",ROC_rt$AUC[2])),
             paste0('AUC at 5 years: ', sprintf("%.03f",ROC_rt$AUC[3]))),
           col=c("green","blue","red"), lwd=2, bty = 'n')
    dev.off()
}
timeRoc(input=riskOut, rocFile="train.ROC.pdf")
rm(list = c("timeRoc"))

###########################
#  PCA和t-SNE分析         #
###########################
# 查看预后模型是否可以准确分出高低风险的病人
# install.packages("Rtsne")
# install.packages("ggplot2")
library(Rtsne)
library(ggplot2)
doPca <- function(input = null, pcaFile = null, tsneFile = null) {
    rt = input
    data=rt[c(3:(ncol(rt)-2))]
    risk=rt[, "riskType"]
    # PCA
    data.pca=prcomp(data, scale. = TRUE)
    pcaPredict=predict(data.pca)
    PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2], risk = risk)	
    pdf(file=pcaFile, height=4.5, width=5.5)   
    p = ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
        scale_colour_manual(name="Risk",  values =c("red", "blue")) +
        theme_bw() +
        theme(plot.margin=unit(rep(1.5,4),'lines')) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
    # t-SNE
    tsneOut = Rtsne(data, dims=2, pca = F, perplexity=(nrow(data) - 1 )/4, verbose=F, max_iter= 500,check_duplicates=F)
    tsneOut
    tsne = data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2], risk=risk)	
    pdf(file=tsneFile, height=4.5, width=5.5)
    p = ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk)) +
        scale_colour_manual(name="Risk",  values =c("red", "blue"))+
        theme_bw()+
        theme(plot.margin=unit(rep(1.5,4),'lines'))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
doPca(input = riskOut[, -1], pcaFile="PCA.pdf", tsneFile="t-SNE.pdf")
rm(doPca)

###################################################################
# 临床相关性分析, 即分析模型得分，基因表达 和 临床信息之间的关系  #
###################################################################
clinicalAnalysisGeneExp <- lassoSigExp[, !(colnames(lassoSigExp) %in% c("time", "alive"))]
clinicalAnalysisGeneExp <- cbind(id = rownames(clinicalAnalysisGeneExp), clinicalAnalysisGeneExp)
oneAnalysisIn <- oneAnalysisIn[, !(colnames(oneAnalysisIn) %in% c("time", "alive"))]
oneAnalysisIn <- cbind(id = rownames(oneAnalysisIn), oneAnalysisIn)
clinicalAnalysisIn <- merge(oneAnalysisIn, clinicalAnalysisGeneExp, by = "id")
clinicalAnalysisIn <- cbind(clinicalAnalysisIn[, !(colnames(clinicalAnalysisIn) %in% c("riskScore"))], riskScore = clinicalAnalysisIn$riskScore)
clinicalAnalysisIn$age_at_index <- ifelse(clinicalAnalysisIn$age_at_index <= 65, "<=65", ">65")
# install.packages("beeswarm")
library(beeswarm)
clinicalNum <- ncol(oneAnalysisIn) -  2                                                       
pFilter <- 0.05                                                         
clinicalOut <- data.frame(gene=colnames(clinicalAnalysisIn[, (clinicalNum+2):ncol(clinicalAnalysisIn)]))
for(clinical in colnames(clinicalAnalysisIn[, 2:(clinicalNum+1)])) {
    xlabel = vector()
    tab1 = table(clinicalAnalysisIn[, clinical])
    labelNum = length(tab1)
    dotCol = c("blue","red")
    if (labelNum == 3) {
        dotCol = c(2,3,4)
    }
    if (labelNum == 4) {
        dotCol = c(2,3,4,5)
    }
    if (labelNum > 4) {
        dotCol = rainbow(labelNum)
    }
    for(i in 1:labelNum) {
        xlabel=c(xlabel, names(tab1[i]) )
    }
    clinicalPvalVector = c()
    for(i in colnames(clinicalAnalysisIn[, (clinicalNum+2):ncol(clinicalAnalysisIn)])) {
        rt1 = rbind(expression = clinicalAnalysisIn[, i], clinical = clinicalAnalysisIn[, clinical])
        rt1 = as.matrix(t(rt1))
        if (labelNum == 2) {
            cliTest <- t.test(expression ~ clinical, data = rt1)
        } else {
            cliTest <- kruskal.test(expression ~ clinical, data = rt1)
        }
        pValue = cliTest$p.value
        stat=round(cliTest$statistic, 3)
        pval = 0
        if (pValue < 0.001) {
            pval = signif(pValue,4)
            pval = format(pval, scientific = TRUE)
        } else {
            pval = sprintf("%.03f", pValue)
        }
        clinicalPvalVector = c(clinicalPvalVector, paste0(stat, "(", pval, ")"))
        if (pValue < pFilter) {
            b = boxplot(expression ~ clinical, data = rt1, outline = FALSE, plot=F)
            yMin=min(b$stats)
            yMax = max(b$stats/5+b$stats)
            n = ncol(b$stats)
            outPdf=paste0(i, ".", clinical, ".pdf")
            pdf(file=outPdf, width = 7, height = 5)
            par(mar = c(4.5,6,3,3))
            ylab = ifelse(i == "riskScore", "Risk score", "Gene expression")
            boxplot(expression ~ clinical, data = rt1, names=xlabel,
                    ylab = ylab, main=paste0(i," (p=",pval,")"), xlab=clinical,
                    cex.main=1.4, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax), outline = FALSE)
            beeswarm(expression ~ clinical, data = rt1, col = dotCol, lwd = 0.1,
                     pch = 16, add = TRUE, corral = "wrap")
            dev.off()
        }
    }
    clinicalOut <- cbind(clinicalOut, clinicalPvalVector)
}
colnames(clinicalOut) = c("id", colnames(clinicalAnalysisIn[, 2:(clinicalNum+1)]))
rm(list = c("clinicalAnalysisGeneExp", "xlabel", "tab1", "labelNum", 
            "pFilter", "clinicalOut", "clinicalNum", "oneAnalysisIn",
            "dotCol", "clinicalPvalVector", "rt1", "pValue", "stat",
            "pval", "b", "yMin", "yMax", "n", "outPdf",
            "ylab", "cliTest", "clinical"))



#####################################
#  ssGSEA进行免疫打分并差异分析     #
#####################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("GSVA")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("GSEABase")
library(GSVA)
library(limma)
library(GSEABase)
# ssGSEA免疫打分
getImmuneScore <- function(expIn=null, gmtFile=null, project=null) {
    rt=expIn
    rt=as.matrix(rt)
    exp=rt
    dimnames=list(rownames(exp), colnames(exp))
    mat=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
    mat=avereps(mat)
    mat=mat[rowMeans(mat)>0,]
    geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
    
    ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, min.sz = 2)
    
    normalize=function(x){
        return((x-min(x))/(max(x)-min(x)))}
    ssgseaOut <<- normalize(ssgseaScore)
}
getImmuneScore(expIn = allGeneExp, gmtFile="D:\\Users\\Lenovo\\Desktop\\毕设\\immune.gmt", project="TCGA")
rm(getImmuneScore)
# 免疫差异分析(免疫细胞和免疫功能)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
# install.packages("reshape2")
# install.packages("ggpubr")
library(limma)
library(reshape2)
library(ggpubr)
scoreCor <- function(riskFile=null, scoreFile=null, project=null){
    data=t(scoreFile)
    risk=riskFile[, -1]
    data=data[rownames(risk), ]
    rt=cbind(data, risk[, c("riskScore","riskType")])
    rt=rt[, -(ncol(rt)-1)]
    
    immCell=c("B_cells","Macrophages","Neutrophils","pDCs",
              "Th1_cells","Th2_cells","TIL","Treg")
    rt1=rt[, c(immCell, "riskType")]
    data=melt(rt1, id.vars=c("riskType"))
    colnames(data)=c("Risk","Type","Score")
    data$Risk=factor(data$Risk, levels=c("low","high"))
    p=ggboxplot(data, x="Type", y="Score", color = "Risk",
                xlab="",ylab="Score",add = "none",palette = c("blue","red") )
    p=p+rotate_x_text(50)
    p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
    
    pdf(file=paste0(project,".immCell.pdf"), width=7, height=6)
    print(p)
    dev.off()
    
    immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
                  "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
                  "MHC_class_I","Parainflammation","T_cell_co-inhibition",
                  "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
    rt1=rt[,c(immFunction,"riskType")]
    data=melt(rt1,id.vars=c("riskType"))
    colnames(data)=c("Risk","Type","Score")
    data$Risk=factor(data$Risk, levels=c("low","high"))
    p=ggboxplot(data, x="Type", y="Score", color = "Risk",
                xlab="",ylab="Score",add = "none",palette = c("blue","red") )
    p=p+rotate_x_text(50)
    p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
    
    pdf(file=paste0(project,".immFunction.pdf"), width=7, height=6)
    print(p)
    dev.off()
}
scoreCor(riskFile=riskOut, scoreFile=ssgseaOut, project="TCGA")

write.table(coefOut, file = "GEO/coef.txt", row.names = T, col.names = T)

