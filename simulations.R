## stIHC 
source("stIHC.R")
library(fda)
library(fdaPDE)
library(mclust)
library(clusterCrit)
library(cluster) 

##################################################################################################################################################################
# Simulation from Section 3.1: Equally Sized Modules
data = read.csv("Section 3.1 Equally Sized Modules.csv")
id = as.matrix(read.csv("groundtruth1.csv", header = F))

stihc = stIHC(data)
ARI = adjustedRandIndex(stihc$label,id) #Compute Adjusted Rand Index


##################################################################################################################################################################
# Simulation from Section 3.2: Imbalanced Modules
data = read.csv("Section 3.2 Imbalanced Modules.csv")
id = as.matrix(read.csv("groundtruth23.csv", header = F))

stihc = stIHC(data)
ARI = adjustedRandIndex(stihc$label,id) #Compute Adjusted Rand Index

##################################################################################################################################################################
# Simulation from Section 3.3: Sparse and Imbalanced Modules
data = read.csv("Section 3.3 Sparse and Imbalanced Modules.csv")
id = as.matrix(read.csv("groundtruth23.csv", header = F))

stihc = stIHC(data)
ARI = adjustedRandIndex(stihc$label,id) #Compute Adjusted Rand Index

