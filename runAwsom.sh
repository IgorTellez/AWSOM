#!/bin/bash

# --------------------------------------------------------------------------------------------------------------------
# Script Name: runAwsom.R
# Description: This script calls and runs the Adaptive Weights Smoothing with Optimized Metrics (AWSOM) algorithm (awsom.R). 
#              The results will be the p-value activity map, z-values and the filtered beta-value activity map.
#
#              It requires the path to the following files obtained by the GLM analysis:
#              - copePath: beta map, also called contrast parameter (COPE) in FSL
#              - residualsPath: residuals
#              - varianceCopePath: COPE variance
#              - maskPath: binary brain mask
#              - numRegressors: indicate with an integer the number of regressors used in the analysis
#              - outputPath: path to the output folder
#
#              If the GLM analysis has been performed with FSL, the necessary files are located in the folder named 'stats' 
#              within the output folder, normally named with the suffix '.feat'. The default names of the required files are: 
#              -cope*.nii.gz
#              -res4d.nii.gz
#              -varcope*.nii.gz
#
#              To run the awsom.R script you need to have R installed, 
#              as well as the following libraries:
#              - fmri: https://cran.r-project.org/web/packages/fmri/index.html
#              - RNifti: https://cran.r-project.org/web/packages/RNifti/index.html
#              - oro.nifti : https://cran.r-project.org/web/packages/oro.nifti/index.html
#              In addition you need to install the fmri_1.9.12.3.tar.gz file from:
#              https://github.com/IgorTellez/AWSOM/blob/e97913dd646b3d0ae9b66983f7e9228ad56d5736/fmri_1.9.12.3.tar.gz
#
#
#              The FSL toolbox needs to be installed to convert the map from p-values to z-values. 
#
#
# Author: Igor Tellez
# Contact: igor.tellez@gmail.com
# 
# Contributor:
#   - Karsten Tabelow (karsten.tabelow@wias-berlin.de) - AWS & AWSOM developer
#
# Created: 2023-12-02
# Last Modified: 2024-05-23
# --------------------------------------------------------------------------------------------------------------------

# Define the Bash variables
copePath=/pathToCopeFile/
residualsPath=/pathToResidualsFile/
varianceCopePath=/pathToVarCopeFile/
maskPath=/pathToMaskFile/
numRegressors=numberOfRegressors # indicates the number of regressors used in the analysis as integer
outputPath=/pathToOutputFolder/

# Run the R script with the Bash variables as arguments
Rscript awsom.R "$copePath" "$residualsPath" "$varianceCopePath" "$maskPath" "$numRegressors" "$outputPath"

# Switching p-values map to z-score map
fslmaths ${outputPath}/AWSOM_pVal.nii -ptoz ${outputPath}/AWSOM_zstat.nii.gz

# End of script
