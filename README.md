You can install this R package using the following two commands in R:

library(devtools)

install_github("rezadrikvandi/enrichedInference")

and then load it using library(enrichedInference)

Details about this R package can be found in the DESCRIPTION file and the main R file. The package implements the enriched inference for high dimensional data. The main function is called "enrichedInference" which applies both the enriched estimation and inference. Also, a test data set is included together with the real data set "riboflavin" analysed in the paper.

To apply the enriched method to a test data set in the package called testdata (with n=100, p=500 and only the first three covariates being truly significant each with coefficient 1 and with no true intercept), try this command:
enrichedInference(X=testdata[,-1], Y=testdata[,1], intercept=FALSE, Xcorr="marginal")

Also, to apply the enriched method to the real data (riboflavin data with n=71 and p=4088) in the package, try the following commands both with some default settings for a quicker run (one with partial correlations and another with marginal corelations):

enrichedInference(X=riboflavin$x, Y=riboflavin$y, k=1, Xcorr="partial", Bonf.correction=TRUE)
enrichedInference(X=riboflavin$x, Y=riboflavin$y, k=1, Xcorr="marginal", Bonf.correction=TRUE)

