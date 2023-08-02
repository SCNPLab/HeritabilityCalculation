# -----------------------------------------------------------------------
# Program: UnivACE.R  
# Univariate Twin Analysis model to estimate causes of variation (ACE) for continuous data
# Matrix style model input - Raw data input
# -----------------------------------------------------------------------
# Yue 2013.12.22
# function: UnivACE: univariate ACE model
# input: 
# namevar1 -- THICK1
# namevar2 -- THICK2
# mydata -- ThickAreadata
# output: 
# matrix[1*27] 
# matrix[1,1:14]= ACE('2ll','obs','ep','obsStat','df','AIC','BIC','a','c','e','a^2','c^2','e^2','p')
# matrix[1,15:27]= AE('2ll','obs','ep','obsStat','df','AIC','BIC','a','c','e','a^2','c^2','e^2')

# --------------------------------------------------------------------------

UnivACE <- function(namevarindex){

nv <<- 1			# number of variables for a twin = 1 in Univariate
ntv <<- 2*nv			# number of variables for a pair = 2* 1 for Univariate

ACEsummary <- matrix(NA,1,11)			# save ('2ll','obs','ep','obsStat','df','albound','aubound',,'clbound''cubound','elbound','eubound')
AEsummary <- matrix(NA,1,11)			# save ('2ll','obs','ep','obsStat','df','albound','aubound',,'clbound''cubound','elbound','eubound')
CEsummary <- matrix(NA,1,11)			# save ('2ll','obs','ep','obsStat','df','albound','aubound',,'clbound''cubound','elbound','eubound')
Esummary <- matrix(NA,1,5)			# save ('2ll','obs','ep','obsStat','df')
ACEpara <- matrix(NA,1,6)			# save ('a','c','e')('a^2','c^2','e^2')
AEpara <- matrix(NA,1,6)			# save ('a','c','e')('a^2','c^2','e^2')
CEpara <- matrix(NA,1,6)			# save ('a','c','e')('a^2','c^2','e^2')
AICsummary <- matrix(NA,1,8)       # save ('ACE_AIC','ACE_BIC','AE_AIC','AE_BIC','CE_AIC','CE_BIC','E_AIC','E_BIC')
result <- matrix(NA,1,72)			# for all, including 'p' for ACE vs AE



# variables to be analyzed
namevar1 <- namevar[namevarindex+1]
namevar2 <- namevar[namevarindex+1+numvar]
selVars <- c(namevar1,namevar2)
mzData <- subset(ThickAreadata, zyg==1, selVars)
dzData <- subset(ThickAreadata, zyg==2, selVars)
# if the other stats should be output, run the following lines
#xx <- paste("m01_01_",namevar1,".Ro",sep="")
# save output
#--------------------------------------------------------------------
#sink(xx, append=FALSE, split=TRUE)
#--------------------------------------------------------------------
# Print Descriptive Statistics
# -----------------------------------------------------------------------

#print(describe(ThickAreadata))
#print(colMeans(mzData,na.rm=TRUE))
#print(cov(mzData,use="complete"))
#print(colMeans(dzData,na.rm=TRUE))
#print(cov(dzData,use="complete"))

#------------------------------------------------------------------------
# Fit Univariate Saturated Model
# -----------------------------------------------------------------------
univTwinSatModel <- mxModel("univTwinSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=10, name="CholMZ" ),
        mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=90, name="expMeanMZ" ),
        mxData( observed=mzData, type="raw" ),
        mxFIMLObjective( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars),
    # Algebra's needed for equality constraints    
        mxAlgebra( expression=expMeanMZ[1,1:nv], name="expMeanMZt1"),
        mxAlgebra( expression=expMeanMZ[1,(nv+1):ntv], name="expMeanMZt2"),
        mxAlgebra( expression=t(diag2vec(expCovMZ)), name="expVarMZ"),
        mxAlgebra( expression=expVarMZ[1,1:nv], name="expVarMZt1"),
        mxAlgebra( expression=expVarMZ[1,(nv+1):ntv], name="expVarMZt2")
    ),
    mxModel("DZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=10, name="CholDZ" ),
        mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=90, name="expMeanDZ" ),
        mxData( observed=dzData, type="raw" ),
        mxFIMLObjective( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars),
    # Algebra's needed for equality constraints    
        mxAlgebra( expression=expMeanDZ[1,1:nv], name="expMeanDZt1"),
        mxAlgebra( expression=expMeanDZ[1,(nv+1):ntv], name="expMeanDZt2"),
        mxAlgebra( expression=t(diag2vec(expCovDZ)), name="expVarDZ"),
        mxAlgebra( expression=expVarDZ[1,1:nv], name="expVarDZt1"),
        mxAlgebra( expression=expVarDZ[1,(nv+1):ntv], name="expVarDZt2")
    ),
    mxAlgebra( MZ.objective + DZ.objective, name="min2sumll" ),
    mxAlgebraObjective("min2sumll")
)

univTwinSatFit <- mxRun(univTwinSatModel, intervals=TRUE)
univTwinSatSumm <- summary(univTwinSatFit)
#print(univTwinSatSumm)

# Generate Saturated Output
# -----------------------------------------------------------------------
parameterSpecifications(univTwinSatFit)
expectedMeansCovariances(univTwinSatFit)
#tableFitStatistics(univTwinSatFit)

        
# Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
 univACEModel <- mxModel("univACE",
    mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=10, label="a11", name="a" ), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=10, label="c11", name="c" ), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=10, label="e11", name="e" ), 
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=a %*% t(a), name="A" ),
        mxAlgebra( expression=c %*% t(c), name="C" ),
        mxAlgebra( expression=e %*% t(e), name="E" ),
    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=A+C+E, name="V" ),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
    # Algebra to compute standardized path estimares and variance components
        mxAlgebra( expression=a%*%iSD, name="sta"),
        mxAlgebra( expression=c%*%iSD, name="stc"),
        mxAlgebra( expression=e%*%iSD, name="ste"),
        mxAlgebra( expression=A/V, name="h2"),
        mxAlgebra( expression=C/V, name="c2"),
        mxAlgebra( expression=E/V, name="e2"),
    # Note that the rest of the mxModel statements do not change for bivariate/multivariate case
    # Matrix & Algebra for expected means vector
        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= 90, label="mean", name="Mean" ),
        mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 
    ),
    mxModel("MZ",
        mxData( observed=mzData, type="raw" ),
        mxFIMLObjective( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selVars )
    ),
    mxModel("DZ", 
        mxData( observed=dzData, type="raw" ),
        mxFIMLObjective( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selVars ) 
    ),
    mxAlgebra( expression=MZ.objective + DZ.objective, name="m2ACEsumll" ),
    mxAlgebraObjective("m2ACEsumll"),
    mxCI(c('ACE.h2', 'ACE.c2', 'ACE.e2'))

)

univACEFit <- mxRun(univACEModel, intervals=TRUE)
univACESumm <- summary(univACEFit)
#univACESumm
#univACEFit$ACE.h2
#univACEFit$ACE.c2
#univACEFit$ACE.e2
#univACEFit$ACE.sta
#univACEFit$ACE.stc
#univACEFit$ACE.ste


# Generate ACE Output
# -----------------------------------------------------------------------
parameterSpecifications(univACEFit)
expectedMeansCovariances(univACEFit)
#tableFitStatistics(univACEFit)
#-----------------------------------------------------------
# Yue2013.12.19
# extract some values
#
# ACEsummary[,1] <- univACESumm[[2]]
# ACEsummary[,2] <- univACESumm[[5]]
# ACEsummary[,3] <- univACESumm[[6]]
# ACEsummary[,4] <- univACESumm[[7]]
# ACEsummary[,5] <- univACESumm[[8]]
#ACEsummary[,6] <- univACESumm[[13]]
#ACEsummary[,7] <- univACESumm[[14]]
# ACEsummary[,6] <- univACESumm[[23]][1]
# ACEsummary[,7] <- univACESumm[[23]][7]
# ACEsummary[,8] <- univACESumm[[23]][2]
# ACEsummary[,9] <- univACESumm[[23]][8]
# ACEsummary[,10] <- univACESumm[[23]][3]
# ACEsummary[,11] <- univACESumm[[23]][9]
ACEsummary[,1] <- univACESumm[[8]]
ACEsummary[,2] <- univACESumm[[11]]
ACEsummary[,3] <- univACESumm[[12]]
ACEsummary[,4] <- univACESumm[[13]]
ACEsummary[,5] <- univACESumm[[14]]
ACEsummary[,6] <- univACESumm[[36]][1,1]
ACEsummary[,7] <- univACESumm[[36]][1,3]
ACEsummary[,8] <- univACESumm[[36]][2,1]
ACEsummary[,9] <- univACESumm[[36]][2,3]
ACEsummary[,10] <- univACESumm[[36]][3,1]
ACEsummary[,11] <- univACESumm[[36]][3,3]

AICsummary[,1] <- univACESumm[[22]]
AICsummary[,2] <- univACESumm[[23]]
#------------------------------------------------------------

# Generate Table of Parameter Estimates using mxEval
pathEstimatesACE <- mxEval(cbind(ACE.sta,ACE.stc,ACE.ste), univACEFit)
#    rownames(pathEstimatesACE) <- 'pathEstimates'
#    colnames(pathEstimatesACE) <- c('a','c','e')
#print(pathEstimatesACE)

varComponentsACE <- mxEval(cbind(ACE.h2,ACE.c2,ACE.e2), univACEFit)
#    rownames(varComponentsACE) <- 'varComponents'
#    colnames(varComponentsACE) <- c('a^2','c^2','e^2')
#print(varComponentsACE)
ACEpara <- c(pathEstimatesACE,varComponentsACE)

# Fit AE model
# -----------------------------------------------------------------------
univAEModel <- mxModel(univACEFit, name="univAE",
    mxModel(univACEFit$ACE,
        mxMatrix( type="Lower", nrow=1, ncol=1, free=FALSE, values=0, label="c11", name="c" ) # drop c at 0
    )
)
univAEFit <- mxRun(univAEModel, intervals=TRUE)
univAESumm <- summary(univAEFit)
#univAESumm

# Generate AE Output
# -----------------------------------------------------------------------
parameterSpecifications(univAEFit)
expectedMeansCovariances(univAEFit)
#tableFitStatistics(univAEFit)
#-----------------------------------------------------------
# Yue2013.12.19
# extract some values
#
AEsummary[,1] <- univAESumm[[8]]
AEsummary[,2] <- univAESumm[[11]]
AEsummary[,3] <- univAESumm[[12]]
AEsummary[,4] <- univAESumm[[13]]
AEsummary[,5] <- univAESumm[[14]]
#AEsummary[,6] <- univAESumm[[13]]
#AEsummary[,7] <- univAESumm[[14]]
AEsummary[,6] <- univAESumm[[36]][1,1]
AEsummary[,7] <- univAESumm[[36]][1,3]
AEsummary[,8] <- univAESumm[[36]][2,1]
AEsummary[,9] <- univAESumm[[36]][2,3]
AEsummary[,10] <- univAESumm[[36]][3,1]
AEsummary[,11] <- univAESumm[[36]][3,3]


AICsummary[,3] <- univAESumm[[22]]
AICsummary[,4] <- univAESumm[[23]]
#------------------------------------------------------------

# Generate Table of Parameter Estimates using mxEval
pathEstimatesAE <- print(round(mxEval(cbind(ACE.sta,ACE.stc,ACE.ste), univAEFit),4))
varComponentsAE <- print(round(mxEval(cbind(ACE.h2,ACE.c2,ACE.e2), univAEFit),4))
#	rownames(pathEstimatesAE) <- 'pathEstimates'
#	colnames(pathEstimatesAE) <- c('a','c','e')
#	rownames(varComponentsAE) <- 'varComponents'
#	colnames(varComponentsAE) <- c('a^2','c^2','e^2')
#print(pathEstimatesAE)
#print(varComponentsAE)
AEpara <- c(pathEstimatesAE,varComponentsAE)

# Fit CE model
# -----------------------------------------------------------------------
univCEModel <- mxModel(univACEFit, name="univCE",
    mxModel(univACEFit$ACE,
        mxMatrix( type="Lower", nrow=1, ncol=1, free=FALSE, values=0, label="a11", name="a" ) # drop a at 0
    )
)
univCEFit <- mxRun(univCEModel, intervals=TRUE)
univCESumm <- summary(univCEFit)
#univCESumm

# Generate CE Output
# -----------------------------------------------------------------------
parameterSpecifications(univCEFit)
expectedMeansCovariances(univCEFit)
#tableFitStatistics(univCEFit)

CEsummary[,1] <- univCESumm[[8]]
CEsummary[,2] <- univCESumm[[11]]
CEsummary[,3] <- univCESumm[[12]]
CEsummary[,4] <- univCESumm[[13]]
CEsummary[,5] <- univCESumm[[14]]
#CEsummary[,6] <- univCESumm[[13]]
#CEsummary[,7] <- univCESumm[[14]]
CEsummary[,6] <- univCESumm[[36]][1,1]
CEsummary[,7] <- univCESumm[[36]][1,3]
CEsummary[,8] <- univCESumm[[36]][2,1]
CEsummary[,9] <- univCESumm[[36]][2,3]
CEsummary[,10] <- univCESumm[[36]][3,1]
CEsummary[,11] <- univCESumm[[36]][3,3]


AICsummary[,5] <- univCESumm[[22]]
AICsummary[,6] <- univCESumm[[23]]

# Generate Table of Parameter Estimates using mxEval
pathEstimatesCE <- print(round(mxEval(cbind(ACE.sta,ACE.stc,ACE.ste), univCEFit),4))
varComponentsCE <- print(round(mxEval(cbind(ACE.h2,ACE.c2,ACE.e2), univCEFit),4))
#	rownames(pathEstimatesCE) <- 'pathEstimates'
#	colnames(pathEstimatesCE) <- c('a','c','e')
#	rownames(varComponentsCE) <- 'varComponents'
#	colnames(varComponentsCE) <- c('a^2','c^2','e^2')
#print(pathEstimatesCE)
#print(varComponentsCE)
CEpara <- c(pathEstimatesCE,varComponentsCE)


# Fit E model
# -----------------------------------------------------------------------
# Note: we call the AE model and drop a, fix it to 0
univEModel <- mxModel(univAEFit, name="univE",
    mxModel(univAEFit$ACE,
        mxMatrix( type="Lower", nrow=1, ncol=1, free=FALSE, values=0, label="a11", name="a" ) # drop a at 0
    )
)
univEFit <- mxRun(univEModel, intervals=TRUE)
univESumm <- summary(univEFit)
#univESumm

Esummary[,1] <- univESumm[[8]]
Esummary[,2] <- univESumm[[11]]
Esummary[,3] <- univESumm[[12]]
Esummary[,4] <- univESumm[[13]]
Esummary[,5] <- univESumm[[14]]

AICsummary[,7] <- univESumm[[22]]
AICsummary[,8] <- univESumm[[23]]

# Generate E Output
# -----------------------------------------------------------------------
parameterSpecifications(univEFit)
expectedMeansCovariances(univEFit)
#tableFitStatistics(univEFit)

# Generates an output table of all submodels (in 'list') compared to Fully Saturated model 
# -----------------------------------------------------------------------
#univACENested <- list(univACEFit, univAEFit, univCEFit, univEFit)
#tableFitStatistics(univTwinSatFit,univACENested)
Sat.fit <- 	rbind(mxCompare(univTwinSatFit, univACEFit),
                                mxCompare(univTwinSatFit, univAEFit)[2,],
	   				  mxCompare(univTwinSatFit, univCEFit)[2,],
	  				  mxCompare(univTwinSatFit, univEFit)[2,])
Sat.fit


# Generates an output table of all submodels (in 'list') compared to previous (Nested models)
# --------------------------------------------------------------------------------------------------------------------------------

Nested.fit <- 	rbind(mxCompare(univACEFit, univAEFit),
	   				  mxCompare(univACEFit, univCEFit)[2,],
	  				  mxCompare(univACEFit, univEFit)[2,],
					  mxCompare(univAEFit, univEFit)[2,])
# print(Nested.fit)



result <- c(ACEsummary, ACEpara,Sat.fit[[9]][2],Sat.fit[[9]][3],Sat.fit[[9]][4],Sat.fit[[9]][5],AEsummary, AEpara, Nested.fit[[9]][2],Nested.fit[[9]][3],Nested.fit[[9]][4],Nested.fit[[9]][5],AICsummary,CEsummary,CEpara,Esummary)
return(result)

# restores it to the console once again
#----------------------------------------------------------------------
#sink()
#----------------------------------------------------------------------
}



