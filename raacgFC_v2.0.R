#!/usr/bin/env Rscript

library(edgeR)
library(Rsubread)
#library(lib.loc="/Users/javier/libs/R-3.3.0/", package="dupliR")
library(lib.loc="/Users/javier/libs/R-3.3.0/", package="dupliR")
#source(file="/Users/javier/Documents/clase/dupliR/R/dupliR.R")

args <- commandArgs(trailingOnly=T)


# Loading FeatureCounts R Object
RDATA <- args[3]		# Argument position 3: fc RData object
FC.OBJ <- args[4]		# Argument position 4: featureCounts counts object name
WD <- args[1]			# Argument position 1: Working directory
SRAFILE <- args[2]		# Argument position 2: Name of the file containing the SRR files name
SN2N <- "xlae.sn2n.txt"	# X. laevis systematic gene name to common name
INTER.GENES <- "inter.genes.txt"
load(RDATA)				# Loading featureCounts RData object
fc.obj <- eval(parse(text=FC.OBJ))		# Allocating FC object name
counts <- fc.obj$counts	# Allocating counts table to the counts variable

# DEBUG!!!
#counts <- counts[,1:8]

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

#//~~ END PROCESSING COUNTS

#//~~ edgeR analysis

# Experiments names
sraData <- read.delim(paste(WD, SRAFILE, sep="/"), header=F, sep="\t")
# groups names
print("edgeR analysis...")
t <- 5
pb <- txtProgressBar(min=0, max=t, style=3)
DGE.counts <- DGEList(counts=counts)
#	CHECKPOINT... !!!!!
#rownames(DGE.counts)
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

#//~~ END EDGER ANALYSIS

#//~~ Z-SCORE CALCULATING
fitted <- DGE.fit$fitted.values
gene_names <- rownames(fitted)
#-----------------------------------------------------------------\\
# Silenced because it is not needed in raacgFC.R				   \\

# Modifying rownames of fitted table
#head(fitted)
#rownames(fitted) <- gsub("gnl\\|", "", rownames(fitted))
#head(fitted)
#rownames(fitted) <- gsub("\\|PAC4GC.*\\|", "", rownames(fitted))
#head(fitted)
#str(fitted)

#sn2n <- read.delim(SN2N, header=F, sep="|")
#sn2n$V1 <- as.character(sn2n$V1)
#sn2n$V2 <- as.character(sn2n$V2)
#str(sn2n)

#gene_names <- c()
#counter <- 1
#t <- length(rownames(fitted))
#pb <- txtProgressBar(min=1, max=t, style=3)

#print("Formating data...")
#for (name in rownames(fitted)) {
#	if(name %in% sn2n$V1) {
#		gene_names <- c(gene_names, sn2n[sn2n$V1 == name, "V2"])
#	}
#	else {
#		gene_names <- c(gene_names, name)
#	}
#	if(counter > pivot){
#		print(counter)
#		pivot <- pivot + 5000
#	}
#	setTxtProgressBar(pb, counter)
#	counter <- counter + 1
#}
#close(pb)
#rm(pb)

#length(rownames(fitted)) == length(gene_names)
#rownames(fitted) <- gene_names									   //
#-----------------------------------------------------------------//

# Sort fitted by rownames
fitted <- fitted[order(rownames(fitted)),]
raac <- readLines(INTER.GENES)
head(fitted)
det <- DET(fitted)
dif <- calDif(det, raac)
corr <- calCor(det, raac)
plotWTPV(dif)
plotCor(corr)


######	SILENCED CODE	##########

.silenced <- function() {
# Filtering fitted by only L and S members
fitted <- fitted[grepl("\\.[LS]$", rownames(fitted)),]
head(fitted)
# Filtering gene_names: only L and S members
gene_names <- gsub("\\.[LS]$", "", gene_names)
gene_names <- sort(gene_names, decreasing=F)
# Selecting duplicates
gene_dup <- gene_names[duplicated(gene_names)]
#gene_sin <- gene_names[!duplicated(gene_names)]

#length(gene_dup)

# Data structures

fitdup <- c()
counter <- 1
diffTable <- c()
corrTable <- c()
row_gene_names <- c()
# TODO. Prepare a Data structure that holds the gene name (a character slot) and expression data of these genes (a data frame slot)
# Data frame slots should be stored at a list slot called listData that is inside of the main structure.


# Progress Bar
t <- length(gene_dup)
pb <- txtProgressBar(min=1, max=t, style=3)
print("Calculating differences and correlation coefficients...")
for (gene in gene_dup) {
#	print(gene)
#	print(fitted[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(fitted)), ])
	tmp_dups <- fitted[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(fitted)), ]
	#print(tmp_dups)
	if (is.null(dim(tmp_dups)) || dim(tmp_dups)[1] == 0 || dim(tmp_dups)[2] == 0) {
		next
	}
	fitdup <- rbind(fitdup, tmp_dups)
	row_gene_names <- c(row_gene_names, gene)
	# Calculating differences of all pairs for all experiments.
	diffs <- c()
	for(i in 1:length(tmp_dups[1,])){
		diffs <- c(diffs, abs(tmp_dups[2,i]-tmp_dups[1,i]))
	}
	diffTable <- rbind(diffTable, diffs)

	
	# Calculating correlation coefficients for every paralogue pair along all experiments.
	corr <- c()
	corr <- cor(tmp_dups[2,], tmp_dups[1,])
	corrTable <- rbind(corrTable, corr)

#	if(counter >= pivot) {
#		print(counter)
#		pivot <- pivot + 200
#	}
	setTxtProgressBar(pb, counter)
	counter <- counter + 1
}
close(pb)
rm(pb)
#print(counter)
rownames(diffTable) <- row_gene_names
colnames(diffTable) <- colnames(fitted)
rownames(corrTable) <- row_gene_names
colnames(corrTable) <- "Pearson_cor"
head(diffTable)
head(corrTable)
head(fitdup)

# Calculating variance, standard deviance and miu from the entire population.

total_sum <- sum(rowSums(diffTable))
miu <- total_sum/((length(colnames(diffTable)) * length(rownames(diffTable))) - 1)
stdev <- sqrt(var(as.vector(diffTable)))
sterr <- sqrt(var(as.vector(diffTable))) / sqrt((length(colnames(diffTable)) * length(rownames(diffTable))))

#print(total_sum)
#print(miu)
#print(stdev)
#print(sterr)

# Calculating Z-scores for every element of diffTable
print("Calculating Z-scores for every element in diffTable")
Z_scoreTable <- c()

# Progress Bar
t <- length(diffTable[,1])
pb <- txtProgressBar(min = 1, max=t, style=3)
for (i in 1:length(diffTable[,1])) {
	Z_scores <- c()
	for (j in 1:length(diffTable[1,])) {
		Z_score <- (diffTable[i,j] - miu) / stdev
		Z_scores <- c(Z_scores, Z_score)
	}
	setTxtProgressBar(pb,i)
	Z_scoreTable <- rbind(Z_scoreTable, Z_scores)
}
close(pb)
rm(pb)

rownames(Z_scoreTable) <- rownames(diffTable)
colnames(Z_scoreTable) <- colnames(diffTable)

# Filtering Regulated-as-a-couple duplicated genes and the other
inter.genes <- readLines(INTER.GENES)
testDup <- c()

# Progress Bar
cat("Filtering Duplicates and Singles from data...\n")
t <- length(rownames(Z_scoreTable))
pb <- txtProgressBar(min=1, max=t, style=3)
counter <- 1
for(gene in rownames(Z_scoreTable)){
	testDup <- c(testDup, ifelse(gene %in% inter.genes, 1, 0))
	setTxtProgressBar(pb, counter)
	counter <- counter + 1
}

close(pb)
rm(pb)
Z_scoreTable <- cbind(Z_scoreTable, testDup)
diffTable <- cbind(diffTable, testDup)
corrTable <- cbind(corrTable, testDup)

Z_scoreTable <- Z_scoreTable[order(Z_scoreTable[,"testDup"]),]
diffTable <- diffTable[order(diffTable[,"testDup"]),]
corrTable <- corrTable[order(corrTable[,"testDup"]),]

cat("Z_scoreTable...\n")
head(Z_scoreTable)
cat("diffTable...\n")
head(diffTable)
cat("corrTable...\n")
head(corrTable)

write.table(x=fitdup, file="fitdup.txt", sep="\t")
write.table(x=diffTable, file="diffTable.txt", sep="\t")
write.table(x=corrTable, file="corrTable.txt", sep="\t")
write.table(x=Z_scoreTable, file="Z_scoreTable.txt", sep="\t")


# Preparing data to plot expression diferences between homeologues

pdf("Express_diference..pdf", width=17, height=7)
par(mfrow=c(1,2))
length(diffTable[diffTable[,"testDup"] == 1, 1])
length(diffTable[diffTable[,"testDup"] == 0, 1])
boxplot(diffTable[diffTable[,"testDup"] == 1, ], main="Differences", ylim=c(0, 10000))
boxplot(diffTable[diffTable[,"testDup"] == 0, ], main="Differences", ylim=c(0, 10000))

dev.off()
pdf("corrTable.pdf", width=4, height=11)
boxplot(corrTable[corrTable[,"testDup"] == 1, 1], corrTable[corrTable[,"testDup"] == 0, 1], main="Pearson correlation")
dev.off()


# Preparing data to plot Wilcox log10(p.values)

pdf("W_pvalues.pdf", width=15, height=7)
pv.greaters <- c()
pv.lesses <- c()
colnames <- colnames(diffTable)
colnames <- colnames[-length(colnames)]
for (i in 1:length(colnames)) {
	raacg <- diffTable[diffTable[,"testDup"] == 1, i]
	dups <- diffTable[diffTable[,"testDup"] == 0, i]
	greater <- wilcox.test(raacg, dups, alternative="greater", conf.level=0.99)
	less <- wilcox.test(raacg, dups, alternative="less", conf.level=0.99)
	pv.greaters <- c(pv.greaters, greater$p.value)
	pv.lesses <- c(pv.lesses, less$p.value)
}
print(pv.greaters)
pv.greaters <- log10(pv.greaters) * -1
pv.greaters[pv.greaters == Inf] <- 0

pv.lesses <- log10(pv.lesses)
pv.lesses[pv.lesses == Inf] <- 0

p.values <- rbind(pv.greaters, pv.lesses)
print(p.values)
colnames(p.values) <- colnames
rownames(p.values) <- c("greater", "less")

print("Printing p.values")
print(p.values)

barplot(p.values, names.arg=colnames, beside=TRUE, col=c("green","red"))


dev.off()

print("Everything is OK!")
}
