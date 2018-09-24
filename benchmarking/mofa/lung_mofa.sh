#!/bin/bash

# run MOFA using the lung SNF datasets

###################
## START EDITING ##
###################

# Input files as plain text format
inFolder="/Users/asingh/Dropbox/Manuscript/diablo_datasets/benchmarking/mofa"
inFiles=( "$inFolder/lung_datasets/lung_mrna.txt" "$inFolder/lung_datasets/lung_mirna.txt" "$inFolder/lung_datasets/lung_cpg.txt")

# Options for the input files
delimiter="\t" # delimiter, such as "\t", "" or " "
header_rows=1 # set to 1 if the files contain row names
header_cols=1 # set to 1 if the files contain column names

# Output file path, please use the .hdf5 extension
outFile=( "$inFolder/lungMOFA.hdf5" )

# Data options
center_features=0   # center the features to zero-mean? (not necessary as long as learnMean=1)
scale_views=0 	    # scale the views to unit variance (not necessary as long as there no massive differences in scale)

# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian gaussian )

# Define view names
views=( mrna mirna cpg )

# Define file with covariates (not implemented yet, please ignore)
# covariatesFile="/tmp/covariates.txt"

# Maximum number of iterations
iter=1000 # we recommend to set this to a large enough value (>1000)

# Convergence criterion
# Recommendation: a 'tolerance' of 0.01 is quite strict and can take a bit of time, for initial testing we recommend increasing it to 0.1
tolerance=0.01 # training will stop when the change in the evidence lower bound (deltaELBO) is smaller than 0.01
nostop=0       # if nostop=1 the training will complete all iterations even if the convergence criterion is met

# Define the initial number of factors and how inactive factors are dropped during training.
# The model automatically removes inactive factors during training if they explain a fraction of variance smaller than 'dropR2'
# Recommendation: 
# (1) If you remove inactive factors (dropR2>0), then the initial number of factors should be large enough
# (2) If you want to get the most strong drivers of variation then we recommend dropR2 to be at least 0.05 (5%), but if you want to capture more subtle sources of variation you should decrease it to 0.01 (1%) or 0.03 (3%)
factors=2   # initial number of facotrs
startDrop=1  # initial iteration to start shutting down factors
freqDrop=1 	 # frequency of checking for shutting down factors 
dropR2=0.00  # threshold on fractionof variance explained

# Define hyperparameters for the feature-wise spike-and-slab sparsity prior 
learnTheta=( 1 1 1 ) 	# 1 means that sparsity is active whereas 0 means the sparsity is inactivated; each element of the vector corresponds to a view
initTheta=( 1 1 1 ) 	# initial value of sparsity levels (1 corresponds to a dense model, 0.5 corresponds to factors ); each element of the vector corresponds to a view
startSparsity=250 		# initial iteration to activate the spike and slab, we recommend this to be significantly larger than 1.

# Learn an intercept term (feature-wise means)?
# Recommendation: always leave it active. If all your views are gaussian you can set this to 0 and center the features, it does not matter. 
# But for non-gaussian views we noticed that this is very useful, so set it to 1
learnIntercept=1

# Random seed 
seed=0 # if 0, the seed is automatically generated using the current time


####################
## FINISH EDITING ##
####################

# Prepare command
cmd='mofa
	--delimiter "$delimiter"
	--inFiles ${inFiles[@]}
	--outFile $outFile
	--likelihoods ${likelihoods[@]}
	--views ${views[@]}
	--iter $iter
	--tolerance $tolerance
	--learnTheta ${learnTheta[@]}
	--initTheta ${initTheta[@]}
	--startSparsity ${startSparsity[@]}
	--factors $factors
	--startDrop $startDrop
	--freqDrop $freqDrop
	--dropR2 $dropR2
	--seed $seed
'

if [[ $header_rows -eq 1 ]]; then cmd="$cmd --header_rows"; fi
if [[ $header_cols -eq 1 ]]; then cmd="$cmd --header_cols"; fi
# if [ -n "$covariatesFile" ]; then cmd="$cmd --covariatesFile $covariatesFile"; fi
if [[ $center_features -eq 1 ]]; then cmd="$cmd --center_features"; fi
if [[ $scale_views -eq 1 ]]; then cmd="$cmd --scale_views"; fi
if [[ $nostop -eq 1 ]]; then cmd="$cmd --nostop"; fi
if [[ $learnIntercept -eq 1 ]]; then cmd="$cmd --learnIntercept"; fi

eval $cmd
