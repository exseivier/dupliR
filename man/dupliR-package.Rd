\name{dupliR-package}
\alias{dupliR-package}
\alias{dupliR}
\docType{package}
\title{
\packageTitle{dupliR}
}
\description{
\packageDescription{dupliR}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{dupliR}
\packageIndices{dupliR}
This package calculates the differences in expression observed along the experiments (libraries)
between homoeologs or/and duplicates. It returns a plot of the Wilcoxon test p-value of comparing
the median of the diferences observed between a pair of genes of RAAC genes and No-RAAC genes for every library.
Also it calculates the Pearson correlation of the expression level between a gene pair along several
experiments/treatments (libraries).
}
\author{
\packageAuthor{dupliR}

Maintainer: \packageMaintainer{dupliR}
}
\references{
No references yet.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
DET function creates an S4 object that holds the gene name and the expression level observed in
both duplicates / homoeologs.
det <- DET(exp)

data(exp) function loads the "exp" dataset.
"exp" must be a data.frame or matrix of the reads counts normalized by gene length whose columnames
and rownames must be set with library name or treatment condition name and gene names respectively.
gene names must be taged at the end with .L or .S in order to diferentiate them and the sufix after
.L or .S must be the same name.

data(raac) function loads the raac dataset (RAAC genes). This datastructure must be a character vector

calDif() function calculates the difference in expression level between every pair of genes.
dif <- calDif(exp)

calCor() function calculates the Pearson correlation of the expression level between every pair
of genes along a set of libraries.
corr <- calCor(exp)

plotWTPV() function calculates and plots the Wilcoxon test p-values of comparing the median
between the RAAC genes and NO-RAAC genes .
plotWTPV(dif, raac)

plotCor() function plots the Perason correlation of expression for every gene pair and separates the data
between RAAC and No-RAAC genes. It also calculates the Wilcoxon test p-values comparing both data sets
(RAAC and No-RAAC genes).
}
