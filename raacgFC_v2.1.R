#!/usr/bin/env Rscript

#	LOADING DEPENDENCIES
library(edgeR)
library(Rsubread)
library(lib.loc="/Users/javier/libs/R-3.3.0/", package="dupliR")

#	ARGUMENTS FROM BASH
args <- commandArgs(trailingOnly=T)

RDATA <- args[3]		# Argument position 3: fc RData object
FC.OBJ <- args[4]		# Argument position 4: featureCounts counts object name
WD <- args[1]			# Argument position 1: Working directory
SRAFILE <- args[2]		# Argument position 2: Name of the file containing the SRR files name
SN2N <- "xlae.sn2n.txt"	# X. laevis systematic gene name to common name
INTER.GENES <- "inter.genes.txt"
DEBUG <- args[5]

#	LOADING .RData
load(RDATA)				# Loading featureCounts RData object

#	EVALUATING CHR. STR. TO ENV. R OBJECT
fc.obj <- eval(parse(text=FC.OBJ))		# Allocating FC object name
counts <- fc.obj$counts	# Allocating counts table to the counts variable

#	DEBUG!!!
if (DEBUG == "DEBUG") {
	counts <- counts[,1:8]
}

#	PARSING DATA
gene_names <- rownames(counts)
colnames <- colnames(counts)
counts <- counts/(fc.obj$annotation$Length / 1000)

head(counts)
if (dim(counts)[1] != length(gene_names)) {
	print("Row dimension does not fit with Gene names number!")
	print("Execution halted")
	q("n")
}
if(dim(counts)[2] != length(colnames)) {
	print("Column dimension does not fit with Experiments names number!")
	print("Execution halted!")
	q("n")
}
dim(counts)
print("Table dimensions OK!")

#//~~ END -> PROCESSING DATA

#################################################################################

#//~~ edgeR ANALYSIS -> BEGINS

# Experiments names
sraData <- read.delim(paste(WD, SRAFILE, sep="/"), header=F, sep="\t")
# groups names
print("edgeR analysis...")
t <- 5
pb <- txtProgressBar(min=0, max=t, style=3)
DGE.counts <- DGEList(counts=counts)
keep <- c()
for (i in 1:length(rownames(DGE.counts))) {
	meet <- c()
	for (j in 1:length(DGE.counts$counts[i,])){
		if (DGE.counts$counts[i,j] > 30){
			meet <- c(meet, TRUE)
		}
		else {
			meet <- c(meet, FALSE)
		}
	}
	if (sum(meet) > 2) {
		keep <- c(keep, TRUE)
	}
	else {
		keep <- c(keep, FALSE)
	}
}

#keep <- DGE.counts$counts[DGE.counts$counts > 1,]
DGE.counts <- DGE.counts[keep, , keep.lib.sizes=FALSE]
setTxtProgressBar(pb, 1)
#	CHECKPOINT... !!!!!
#head(DGE.counts$counts)
group <- c()
for (name in rownames(DGE.counts$samples)) {
	group <- c(group, as.character(sraData[sraData$V1 == name, 2]))
}
DGE.counts$samples$group <- group
#DGE.counts$samples

design <- model.matrix(~0+group)
#design
rownames(design) <- rownames(DGE.counts$samples)
setTxtProgressBar(pb, 2)
#design

DGE.counts <- calcNormFactors(DGE.counts)
setTxtProgressBar(pb, 3)
DGE.counts <- estimateDisp(DGE.counts, design)
setTxtProgressBar(pb, 4)
DGE.fit <- glmQLFit(DGE.counts, design)
setTxtProgressBar(pb, 5)

print("")
print("Summary of fitted and normalized data")
summary(DGE.fit$fitted.values)

#//~~ END -> edgeR ANALYSIS

###############################################################################

#//~~ dupliR ANALYSIS -> BEGINS
fitted <- DGE.fit$fitted.values
gene_names <- rownames(fitted)
# Sort fitted by rownames
fitted <- fitted[order(rownames(fitted)),]
# Regulated As A Couple genes
raac <- readLines(INTER.GENES)
# Check point!!!
head(fitted)
# Creating DET S4 object
det <- DET(fitted)
# Calculating differences between homoeologs
dif <- calDif(det, raac)
# Calculating Pearson correlation between homoeologs
corr <- calCor(det, raac)
# Plotting data.
plotWTPV(dif)
plotCor(corr)

#//~~ END -> dupliR ANALYSIS
