# ------------------------------------------------------------------------
# Script Name: awsom.R
# Description: This script implements the Adaptive Weights Smoothing with 
#              Optimized Metrics (AWSOM) algorithm. 
# 
# Author: Karsten Tabelow
# Contact: arsten.tabelow@wias-berlin.de
# 
# Contributor:
#   - Igor Tellez (igor.tellez@gmail.com) 
#
# Created: 2024-07-31
# Last Modified: 2024-07-31
# ------------------------------------------------------------------------


library(fmri)
library(oro.nifti)

awsom <- function(cbetaFile, # FSL cope1.nii.gz
                   resFile,   # FSL res4d.nii.gz
                   varFile,   # FSL varcope1.nii.gz
                   maskFile,  # mask file
                   npar,      # number of regressors
                   hmax = 4,  # AWSOM bandwidth
                   outputPath) {

  # Read the NIfTI files
  cbeta <- readNIfTI(cbetaFile, reorient = FALSE)
  res <- readNIfTI(resFile, reorient = FALSE)
  mask <- as.logical(readNIfTI(maskFile, reorient = FALSE))
  
  # Function definitions
  spm <- list()
  spm$var <- readNIfTI(varFile, reorient = FALSE)
  ddim <- dim(cbeta)
  dimt <- dim(res)[4]
  dim(mask) <- ddim
  res <- aperm(res, c(4, 1, 2, 3))
  lags <- c(5, 5, 3)
  corr <- aws::residualSpatialCorr(res, mask, lags = c(5, 5, 3), compact = FALSE)
  scale <- max(abs(range(res)))/32767
  lags <- c(5, 5, 3)
  bw <- optim(c(2, 2, 2),
              fmri:::corrrisk,
              method = "L-BFGS-B",
              lower = c(.59, .59, .59),
              upper = c(10, 10, 10),
              lag = lags,
              data = corr)$par
  bw[bw < .6] <- 0
  dim(bw) <- c(1, 3)
  rxyz <- voxeld / bw
  rxyz[rxyz > 1] <- 1
  dim(rxyz) <- c(1, 3) # nothing changes
  spm$beta <- array(c(rep(1, prod(ddim)), cbeta), dim = c(ddim, 2))
  spm$cbeta <- cbeta
  rsum <- apply(res[, , , ], 2:4, sum)
  mask[rsum==0] <- F
  rsd <- apply(res[, , , ], 2:4, sd)
  mask[rsd==0] <- F
  spm$mask <- mask
  spm$residuals <- writeBin(as.integer(res/scale), raw(), 2)
  spm$resscale <- scale
  spm$maskOnly <- FALSE
  spm$arfactor <- array(0.2, dim = ddim)
  spm$rxyz <- rxyz
  spm$scorr <- corr
  spm$weights <- c(1, 1, 1)
  spm$dim <- c(ddim, dimt)
  spm$hrf <- rep(1, dimt)
  spm$bw <- bw
  spm$df <- dimt - npar
  spm$call <- NULL
  spm$roixa <- 1
  spm$roixe <- ddim[1]
  spm$roiya <- 1
  spm$roiye <- ddim[2]
  spm$roiza <- 1
  spm$roize <- ddim[3]
  spm$roit <- 1
  spm$header <- list()
  spm$format <- "NIFTI"
  class(spm) <- c("fmridata", "fmrispm")
  object.size(spm)
  spm <- condensefMRI(spm, mask)
  object.size(spm)
  adjPar <- 0.8
  
  # AWSOM filtering
  spm.smooth <- fmri.smooth(spm, hmax = hmax, adaptation = "fullaws", ladjust = adjPar)
  
  # Writing AWSOM p-value map
  pvalueawsom <- fmri.pvalue(spm.smooth)
  tempData <- pvalueawsom$pvalue
  tempData[is.na(tempData)] <- 0
  pValAWSOM <- mask
  pValAWSOM@.Data <- tempData
  writeNifti(pValAWSOM, paste(outputPath, "/AWSOM_pVal", sep = ""))
  
  # Writing AWSOM COPE map
  cope <- mask
  cope@.Data <- spm.smooth$cbeta
  writeNifti(cope, paste(outputPath, "/AWSOM_cope", sep = ""))
  
}
