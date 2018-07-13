# INSTALLING LIBRARIES
#install.packages("devtools")
#install.packages("roxygen2")

# LIBRARIES
library(devtools)
library(roxygen2)

# PROCESS DOCUMENTATION
setwd("mrRevBayes")
document()

# INSTALL MY OWN PACKAGE
setwd("..")
install("mrRevBayes")
