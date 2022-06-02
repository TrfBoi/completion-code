geoPath <- "D:\\Users\\Lenovo\\Desktop\\毕设\\code\\RCode\\GEO\\GPL"
# 
# # 基因id转换文件
# annoMatrix <- read.table(file = paste(geoPath, "570", ".txt", sep = ""), sep = "\t", header = T, fill = T, quote = "")
# annoMatrix <- annoMatrix[, c("ID", "Gene.Symbol")]
# annoMatrix <- annoMatrix[!duplicated(annoMatrix$ID, fromLast=TRUE), ]
# rownames(annoMatrix) <- annoMatrix[, 1]
# annoMatrix <- annoMatrix[annoMatrix$Gene.Symbol != "" & unlist(lapply(strsplit(annoMatrix$Gene.Symbol, split = "///"), function(x){ length(as.array(x)) == 1 })), ]
# # 表达矩阵文件
# geoAllExprMatrix <- read.table(file = paste(geoPath, "570", "_matrix.txt", sep = ""), sep = "\t", header = T)
# geoAllExprMatrix[, "ID_REF"] <- annoMatrix[geoAllExprMatrix[, "ID_REF"], "Gene.Symbol"]
# geoAllExprMatrix <- aggregate(x = geoAllExprMatrix[, 2:ncol(geoAllExprMatrix)], by = list(geoAllExprMatrix[, 1]), FUN = mean)
# rownames(geoAllExprMatrix) <- geoAllExprMatrix[, 1]
# geoAllExprMatrix <- geoAllExprMatrix[, -1]
# # 生存信息文件
# geoAllSurvMatrix <- read.table(file = paste(geoPath, "570", "_surv.txt", sep = ""), sep = "\t", header = T, row.names = 1)
# nrow(geoAllExprMatrix)

# 基因id转换文件
# annoMatrix <- read.table(file = paste(geoPath, 96, ".txt", sep = ""), sep = "\t", header = T, fill = TRUE, quote = "")
# annoMatrix <- annoMatrix[, c("ID", "Gene.Symbol")]
# annoMatrix <- annoMatrix[!duplicated(annoMatrix$ID, fromLast=TRUE), ]
# rownames(annoMatrix) <- annoMatrix[, 1]
# annoMatrix <- annoMatrix[annoMatrix$Gene.Symbol != "" & unlist(lapply(strsplit(annoMatrix$Gene.Symbol, split = "///"), function(x){ length(as.array(x)) == 1 })), ]
# # 表达矩阵文件
# geoExprMatrix <- read.table(file = paste(geoPath, 96, "_matrix.txt", sep = ""), sep = "\t", header = T)
# geoExprMatrix[, "ID_REF"] <- annoMatrix[geoExprMatrix[, "ID_REF"], "Gene.Symbol"]
# geoExprMatrix <- aggregate(x = geoExprMatrix[, 2:ncol(geoExprMatrix)], by = list(geoExprMatrix[, 1]), FUN = mean)
# rownames(geoExprMatrix) <- geoExprMatrix[, 1]
# geoExprMatrix <- geoExprMatrix[, -1]
# print(nrow(geoExprMatrix))
# uniqueGene <- rownames(geoExprMatrix)[duplicated(c(rownames(geoExprMatrix), rownames(geoAllExprMatrix)), fromLast = T)]
# print(length(uniqueGene))
# geoAllExprMatrix <- geoAllExprMatrix[uniqueGene, ]
# geoExprMatrix <- geoExprMatrix[uniqueGene, ]
# geoAllExprMatrix <- cbind(geoAllExprMatrix, geoExprMatrix)
# 生存信息文件
# geoSurvMatrix <- read.table(file = paste(geoPath, 96, "_surv.txt", sep = ""), sep = "\t", header = T, row.names = 1)
# geoAllSurvMatrix <- cbind(geoAllSurvMatrix, geoSurvMatrix)

# 基因id转换文件
annoMatrix <- read.table(file = paste(geoPath, 570, ".txt", sep = ""), sep = "\t", header = T, fill = TRUE, quote = "")
annoMatrix <- annoMatrix[, c("ID", "Gene.Symbol")]
annoMatrix <- annoMatrix[!duplicated(annoMatrix$ID, fromLast=TRUE), ]
rownames(annoMatrix) <- annoMatrix[, 1]
annoMatrix <- annoMatrix[annoMatrix$Gene.Symbol != "" & unlist(lapply(strsplit(annoMatrix$Gene.Symbol, split = "///"), function(x){ length(as.array(x)) == 1 })), ]
# 表达矩阵文件
geoExprMatrix <- read.table(file = paste(geoPath, 570, "_matrix.txt", sep = ""), sep = "\t", header = T)
geoExprMatrix[, "ID_REF"] <- annoMatrix[geoExprMatrix[, "ID_REF"], "Gene.Symbol"]
geoExprMatrix <- aggregate(x = geoExprMatrix[, 2:ncol(geoExprMatrix)], by = list(geoExprMatrix[, 1]), FUN = mean)
rownames(geoExprMatrix) <- geoExprMatrix[, 1]
geoExprMatrix <- geoExprMatrix[, -1]
print(nrow(geoExprMatrix))
# uniqueGene <- rownames(geoExprMatrix)[duplicated(c(rownames(geoExprMatrix), rownames(geoAllExprMatrix)), fromLast = T)]
# print(length(uniqueGene))
# geoAllExprMatrix <- geoAllExprMatrix[uniqueGene, ]
# geoExprMatrix <- geoExprMatrix[uniqueGene, ]
# geoAllExprMatrix <- cbind(geoAllExprMatrix, geoExprMatrix)
# 生存信息文件
geoSurvMatrix <- read.table(file = paste(geoPath, 570, "_surv.txt", sep = ""), sep = "\t", header = T, row.names = 1)
# geoAllSurvMatrix <- cbind(geoAllSurvMatrix, geoSurvMatrix)
rm(list = c("annoMatrix", "uniqueGene", "geoPath"))
tmp <- geoSurvMatrix
tmp <- apply(tmp, 2, function(x){
    unlist(strsplit(x, split = ";"))[3:4]
})
tmp <- rbind(apply(as.matrix(tmp[1, ]), 1, function(x){
    unlist(strsplit(x, split = " "))[4]
}), apply(as.matrix(tmp[2, ]), 1, function(x){
    unlist(strsplit(x, split = ":"))[2]
}))
tmp[2, ] <- as.integer(tmp[2, ])
tmp[1, ] <- as.integer(tmp[1, ])
geoSurvMatrix <- tmp
rownames(geoSurvMatrix) <- c("time", "alive")
geoSurvMatrix <- t(geoSurvMatrix)
rm(tmp)
write.table(x = geoExprMatrix, file = "GEO/geoExprMatrix.txt", sep = "\t", row.names = T, col.names = T)
write.table(x = geoSurvMatrix, file = "GEO/geoSurvMatrix.txt", sep = "\t", row.names = T, col.names = T)


#######################################################################################################################################


geoExprMatrix <- read.table(file = "GEO/geoExprMatrix.txt", sep = "\t", row.names = 1, header = T)
geoSurvMatrix <- read.table(file = "GEO/geoSurvMatrix.txt", sep = "\t", row.names = 1, header = T)
coefOut <- read.table(file = "GEO/coef.txt", row.names = 1, header = T)

geoSurvMatrix$time <- (geoSurvMatrix$time) / 365
geoExprMatrix <- geoExprMatrix[rownames(coefOut), ]
geoExprMatrix <- t(geoExprMatrix)
geoExprMatrix <- geoExprMatrix[rownames(geoSurvMatrix), ]
riskScore <- c()
for (i in 1:nrow(geoExprMatrix)) {
    riskScore <- c(riskScore, sum(geoExprMatrix[i, ] * coefOut$coef))
}
riskType <- as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
riskOut <- cbind(geoSurvMatrix, riskScore = riskScore, riskType = riskType)
# riskOut <- riskOut[riskOut$time > 0.1,]
#############
## 生存分析 #
#############
library(survminer)
diff <- survdiff(Surv(time, alive) ~ riskType, data = riskOut)
pValue <- 1-pchisq(diff$chisq, df = 1)
pValue <- signif(pValue,4)
pValue <- format(pValue, scientific = TRUE)
fit <- survfit(Surv(time, alive) ~ riskType, data = riskOut)
pdf(file="GEO/survivorshipCurve.pdf", onefile = FALSE, width = 5.5, height =5)
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

#################################
# 绘制风险曲线（风险得分相关图）#
#################################
# install.packages("pheatmap")
par()
# library(pheatmap)
# riskOut <- riskOut[order(riskOut$riskScore), ]
# riskClass <- riskOut[, "riskType"]
lowLength <- length(riskType[riskType == "low"])
# highLength <- length(riskClass[riskClass == "high"])
# line <- riskOut[, "riskScore"]
# line[line > 10] <- 10 # 降低x轴长度
# pdf(file="riskScore.pdf", width = 10, height = 4)
# plot(line,
#      type="p",
#      pch=20,
#      xlab="Patients (increasing risk socre)",
#      ylab="Risk score",
#      col=c(rep("green", lowLength),
#            rep("red", highLength)))
# abline(h = median(riskOut$riskScore), v=lowLength, lty=2)
# legend("topleft", c("High risk", "low Risk"), bty="n", pch=19, col=c("red","green"), cex=1.2)
# rm(list = c("riskClass", "highLength", "line"))
# dev.off()
color <- as.vector(riskOut$alive)
color[color == 1] <- "red"
color[color == 0] <- "green"
pdf(file="GEO/survStat.pdf", width = 10, height = 4)
plot(riskOut$time,
     pch = 19,
     xlab = "Patients (increasing risk socre)",
     ylab = "Survival time (years)",
     col = color)
legend("topright", c("Dead", "Alive"), bty="n", pch=19, col=c("red","green"), cex=1.2)
abline(v=lowLength, lty=2)
rm(list = c("color", "lowLength"))
dev.off()
# heatMapIn <- riskOut[c(4:(ncol(riskOut)-2))]
# heatMapIn <- t(heatMapIn)
# annotation <- data.frame(type=riskOut[, ncol(riskOut)])
# rownames(annotation) <- rownames(riskOut)
# pdf(file="heatmap.pdf", width = 10, height = 4)
# pheatmap(heatMapIn, 
#          annotation = annotation, 
#          cluster_cols = FALSE,
#          fontsize_row = 11,
#          show_colnames = F,
#          fontsize_col = 3,
#          color = colorRampPalette(c("green", "red"))(50) )
# dev.off()
# rm(list = c("heatMapIn", "annotation"))

############
# ROC 分析 #
############
# install.packages("survivalROC")
# library(survivalROC)
oneAnalysisIn <- riskOut
rocCol <- rainbow(1)
aucText <- c()
par()
pdf("GEO/ROC single.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc <- survivalROC(Stime=oneAnalysisIn$time, status=oneAnalysisIn$alive, marker = oneAnalysisIn$riskScore, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1],
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText, paste0("risk score"," (AUC=", sprintf("%.3f", roc$AUC), ")"))
abline(0,1)
j = 1
oneAnalysisIn <- as.data.frame(oneAnalysisIn)
for (i in colnames(oneAnalysisIn[, 3])) {
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
timeRoc(input=riskOut, rocFile="GEO/train.ROC.pdf")
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
    data=rt[, -ncol(rt)]
    print(ncol(data))
    for (i in 1:ncol(data)) {
        data[, i] <- as.numeric(data[, i])
    }
    risk=riskOut$riskType
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
doPca(input = geoExprMatrix, pcaFile="GEO/PCA.pdf", tsneFile="GEO/t-SNE.pdf")
rm(doPca)

