scCFGL

Description:
scCFGL is a method designed to estimate gene coexpression networks from single-cell RNA-sequencing data. Building upon scLink and CFGL, the method is based on the Gaussian Graphical Model (GGM), where genes are represented as nodes, and their coexpression relationships are captured as edges in the estimated networks. This method is specifically tailored for single-cell gene expression data and can jointly estimate networks while accounting for heterogeneity across biological conditions. The resulting networks allow for meaningful comparisons of gene relationships across different biological states, enabling enhanced biological insights from single-cell genomics data.

Installation:
To install and use scCFGL, you need to have R installed along with the necessary dependencies.

Install from GitHub:
To install scCFGL from GitHub, run the following commands in your R console:

# Install devtools if you don't have it already
install.packages("devtools")

# Install scCFGL
devtools::install_github("esb5324/scCFGL")

Dependencies:
matrixcalc
glmnet
parallel
Make sure these packages are installed before using the package:

install.packages(c("matrixcalc", "glmnet", "parallel"))

To load the package, use:

# Load the package
library(scCFGL)

Refer to the package documentation and R scripts for more details on available functions, including the sccfgl() function.

Authors:
Elle Tang (Author and Creator)
For more details, feel free to contact me at elletang22@gmail.com

Acknowledgements:
scLink (https://github.com/Vivianstats/scLink) and CFGL (https://github.com/Yafei611/CFGL)
Advisor: Dr. Qunhua Li
