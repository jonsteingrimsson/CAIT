# Goal: Load required libraries
#       Could add more libraries in a similar fashion

if (!require(MASS)) {install.packages("MASS"); library(MASS)} 
if (!require(rpart)) {install.packages("rpart"); library(rpart)} 
if (!require(randomForestSRC)) {install.packages("randomForestSRC"); library(randomForestSRC)} 
if (!require(mgcv)) {install.packages("mgcv"); library(mgcv)} 
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)} 
if (!require(plm)) {install.packages("plm"); library(plm)} 
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)} 
if (!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)} 
# if (!require(libcoin)) {install.packages("libcoin"); library(libcoin)} 
# if (!require(inum)) {install.packages("inum"); library(inum)} 
if (!require(partykit)) {install.packages("partykit"); library(partykit)} 
if (!require(randomForest)) {install.packages("randomForest"); library(randomForest)} 
