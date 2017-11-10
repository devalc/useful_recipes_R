############################################################################################################

#Created on Mon Sep 25 14:53:29 2017
#Object: Check if list of packages is already installed and if not, install all the missing ones	#
#@author: Chinmay Deval												#

############################################################################################################

#clear environment
rm(list = ls())

#clear console
cat("\014")

# function to check if package exists and if not install it.
inst_pkg <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "gridExtra", "cowplot", "trend", "plyr", "XLConnect", "RColorBrewer", "scales", "grid")
inst_pkg(packages)
