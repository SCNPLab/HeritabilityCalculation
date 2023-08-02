# -----------------------------------------------------------------------
# Program: UnivACE.R  
# Univariate Twin Analysis model to estimate causes of variation (ACE) for continuous data
# Matrix style model input - Raw data input
# -----------------------------------------------------------------------
# Yue Cui 2013.12

t1 <- Sys.time()
args <- commandArgs(TRUE)
#--------------------------------------------------------------------
# WD <- args[1]
WD <- "F:\\IPCAS_TWIN\\CBF\\SmoothedData\\20220620\\Info"
if (!is.null(WD)) setwd(WD)
(WD <- getwd())
#------------------------------------------------------
require(OpenMx)
require(psych)
source("GenEpiHelperFunctions.R")
source("FunR_Univ_sub.R")		
#-----------------------------------------------------

# Prepare Data
#----------------------------------------------------------------------------------------------------------
# Reads data from csv spreadsheet in which mis val were recoded to 'NA'
# Variabels: Famid ADHD1 IQ1 ADHD2 IQ2 zyg (1=MZ, 2=DZ)
# ----------------------------------------------------------------------------------------------------------
# FileName="twins_data_CBF_BA_forACE.csv"
# filelist <- list.files(path=WD,pattern="twins_data_CBF*",full.names=T,recursive=FALSE)
filelist <- list.files(path=WD,pattern="twins_data_CBF_HOA_whole_withSmoo*",full.names=T,recursive=FALSE)

for (FileName in filelist){
# ThickAreadata <<- read.table (args[2], header=T, sep=',')
ThickAreadata <<- read.table (FileName, header=T, sep=',')
namevar <<- names (ThickAreadata)
#namevar

numvar <- (length(namevar)-2)/2			# should input number of variables

# initialize output matrix
Allresults <- matrix(NA,numvar,72)
summary_name_ACE <- c('2ll','obs','ep','obsStat','df','albound','aubound','clbound', 'cubound','elbound','eubound','a','c','e','a^2','c^2','e^2','p_Sat_ACE','p_Sat_AE','p_Sat_CE','p_Sat_E')
summary_name_AE <- c('2ll','obs','ep','obsStat','df','albound','aubound','clbound','cubound','elbound','eubound','a','c','e','a^2','c^2','e^2','p_ACE_AE','p_ACE_CE','p_ACE_E','p_AE_E','ACE_AIC','ACE_BIC','AE_AIC','AE_BIC','CE_AIC','CE_BIC','E_AIC','E_BIC')
summary_name_CE <- c('2ll','obs','ep','obsStat','df','albound','aubound','clbound','cubound','elbound','eubound','a','c','e','a^2','c^2','e^2')
summary_name_E <- c('2ll','obs','ep','obsStat','df')
# call UnivACE fun
values <- 1:numvar
x <- lapply(values,UnivACE)

for (i in 1:numvar){
y <- unlist(x[i])
Allresults[i,] <- y
}

# output
ACEresults <- Allresults[,1:21]
AEresults <- Allresults[,22:50]
CEresults <- Allresults[,51:67]
Eresults <- Allresults[,68:72]
colnames(ACEresults) <- summary_name_ACE
colnames(AEresults) <- summary_name_AE
colnames(CEresults) <- summary_name_CE
colnames(Eresults) <- summary_name_E
rownames(ACEresults) <- namevar[2:(numvar+1)]
rownames(AEresults) <- namevar[2:(numvar+1)]
rownames(CEresults) <- namevar[2:(numvar+1)]
rownames(Eresults) <- namevar[2:(numvar+1)]
t2 <- Sys.time()
runtime <- difftime(t2,t1,units="mins")
print(runtime)		#Time difference of **mins
# save output
outfile1 <- gsub(".csv$", "_UnivACE.csv", FileName);
#outfile1 <- paste(strsplit(args[2],".csv"),"_ACE.csv",sep="")
write.csv(ACEresults, outfile1, row.names=TRUE)
outfile2 <- gsub(".csv$", "_UnivAE.csv", FileName);
#outfile2 <- paste(strsplit(args[2],".csv"),"_AE.csv",sep="")
write.csv(AEresults, outfile2, row.names=TRUE)
# outfile3 <- gsub(".csv$", "_Univruntime.csv",FileName );
#outfile3 <- paste(strsplit(args[2],".csv"),"_runtime.csv",sep="")
# write.csv(runtime, outfile3)
outfile5 <- gsub(".csv$", "_UnivCE.csv", FileName);
#outfile2 <- paste(strsplit(args[2],".csv"),"_CE.csv",sep="")
write.csv(CEresults, outfile5, row.names=TRUE)
outfile6 <- gsub(".csv$", "_UnivE.csv", FileName);
#outfile2 <- paste(strsplit(args[2],".csv"),"_E.csv",sep="")
write.csv(Eresults, outfile6, row.names=TRUE)
#warnings()
} #for (FileName in filelist){
