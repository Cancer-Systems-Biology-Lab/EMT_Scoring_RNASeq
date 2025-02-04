library(dplyr)
library(data.table)

source("./counts_to_TPM.R")
source("./EMT_score_func.R")

## Folder with all the raw read counts
setwd("./Data")

fileList = list.files(pattern = "*.tsv")

counts = lapply(fileList, read.delim, header = T)
gseIDs = sapply(strsplit(fileList, split = "_"), function(x) x[1])

corMat = matrix(0, nrow = length(gseIDs), ncol = 7)
colnames(corMat) = c("GSEID","76GS-MLR_Cor","76GS-MLR_Pval","76GS-KS_Cor","76GS-KS_Pval","KS-MLR_Cor","KS-MLR_Pval")
for(dataNum in 1:length(counts)){
	log2_TPM = countToTpm(counts[[dataNum]], gseIDs[dataNum])
	MA_val = rnaToMA(log2_TPM)
	EMT76GS_Score = EMT76GS(MA_val)
	KS_score = KSScore(MA_val)
	MLR_score = mlrEMTPred(MA_val, gseIDs[dataNum])
    writeEMTscore(gseIDs[dataNum], EMT76GS_Score, KS_score, MLR_score)
    corMat[dataNum, ] = c(gseIDs[dataNum],all_scoreCor(list(EMT76GS_Score[,1], MLR_score[,2], KS_score[,1])))
}

fwrite(corMat, "../Output/EMT_Score_Correlations.tsv",row.names=F,sep='\t')
