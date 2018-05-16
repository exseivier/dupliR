# OBJECTS
setClass("pair", representation(geneName="character", df="data.frame"),
validity = function(object){
	# Validating that colnames length is equal to geneName length
	errors <- character()
	len_geneName <- length(rownames(object@df))
	if (len_geneName > 2) {
		errors <- c(errors, "gene names length is higher than 2")
	}
	if (length(errors) == 0) TRUE else errors
})

setClass("ExpMet", representation(ExpName="character"))
setClass("pairList", representation(ExpName="ExpMet", Data="list"),
validity = function(object){
	# Validating ExpName length is equal to columns length of every item in Data
	errors <- character()
	len_ExpName <- length(object@ExpName@ExpName)
	for (item in object@Data) {
		len_cols_df <- length(colnames(item@df))
		if (len_ExpName != len_cols_df) {
			errors <- c(errors, paste("Expermient names length does not match rownames in data frame in ", item@geneName, sep=""))
		}
	}
	if (length(errors) == 0) TRUE else errors
})


# METHODS

# DET: Object creator
setGeneric("DET", function(table) standardGeneric("DET"))
setMethod("DET", signature("data.frame"),
	function(table) {
		if (is.null(colnames(table)) || is.null(rownames(table))) {
			return("No rownames or/and colnames in dataset")
		}
		experiment_names <- new("ExpMet", ExpName=colnames(table))
		table <- table[grepl("\\.[LS]$", rownames(table)),]
		table <- table[order(rownames(table)),]
		genes <- rownames(table)
		genes <- gsub("\\.[LS]$", "", genes)
		genes <- sort(genes)
		genes <- genes[duplicated(genes)]
		genes <- unique(genes)
		table <- table[gsub("\\.[LS]$","", rownames(table)) %in% genes,]
		dataList <- list()
		min <- 0
		max <- length(genes)
		cat("Creating DET object\n")
		pb <- txtProgressBar(min=min, max=max, style=3)
		genes_out <- c()
		for (gene in genes) {
	#		print(table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)),])
			if (length( table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), 1]) > 2) {
				min <- min + 1
				setTxtProgressBar(pb, min)
				next
			}
			else {
				data <- new("pair", geneName=gene, df=table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), ])
				table <- table[-grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), ]
				dataList <- c(dataList, data)
				min <- min + 1
				setTxtProgressBar(pb, min)
				genes_out <- c(genes_out, gene)
			}
		}
		object <- new("pairList", ExpName=experiment_names, Data=dataList)
		names(object@Data) <- genes_out
		cat("\nDET object created.\n")
		object
})

setMethod("DET", signature("matrix"),
	function(table) {
		if (is.null(colnames(table)) || is.null(rownames(table))) {
			return("No rownames or/and colnames in dataset")
		}
		table <- as.data.frame(table)
		experiment_names <- new("ExpMet", ExpName=colnames(table))
		table <- table[grepl("\\.[LS]$", rownames(table)),]
		table <- table[order(rownames(table)),]
		genes <- rownames(table)
		genes <- gsub("\\.[LS]$", "", genes)
		genes <- sort(genes)
		genes <- genes[duplicated(genes)]
		genes <- unique(genes)
		table <- table[gsub("\\.[LS]$","", rownames(table)) %in% genes,]
		min <- 0
		max <- length(genes)
		cat("Creating DET object\n")
		pb <- txtProgressBar(min=min, max=max, style=3)
		dataList <- list()
		genes_out <- c()
		for (gene in genes) {
	#		print(table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)),])
			if (length( table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), 1]) > 2) {
				min <- min + 1
				setTxtProgressBar(pb, min)
				next
			}
			else {
				data <- new("pair", geneName=gene, df=table[grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), ])
				table <- table[-grepl(paste(gene, "\\.[LS]$", sep=""), rownames(table)), ]
				dataList <- c(dataList, data)
				min <- min + 1
				setTxtProgressBar(pb, min)
				genes_out <- c(genes_out, gene)
			}
		}
		object <- new("pairList", ExpName=experiment_names, Data=dataList)
		names(object@Data) <- genes_out
		cat("\nDET object created.\n")
		object
})

# calDif: Calculates differences between homoeologs/duplicates
setGeneric("calDif", function(object, raac) standardGeneric("calDif"))
setMethod("calDif", signature("pairList", "character"),
	function(object, raac) {
		result <- c()
		genes <- c()
		cat("\nCalculating absolute differences between homoeologs/duplicates.\n")
		min <- 0
		max <- length(object@Data)
		pb <- txtProgressBar(min=min, max=max, style=3)
		for (item in object@Data){
			d <- calDif(item, raac)
			result <- rbind(result, d)
			genes <- c(genes, item@geneName)
			min <- min + 1
			setTxtProgressBar(pb, min)
		}
		rownames(result) <- genes
		colnames(result) <- c(object@ExpName@ExpName, "testDup")
		cat("\n")
		result
})

setMethod("calDif", signature("pair", "character"),
	function(object, raac) {
		result <- c()
		for (i in 1:length(object@df[1,])) {
			d <- abs(object@df[1,i] - object@df[2,i])
			result <- c(result, d)
		}
		if (object@geneName %in% raac) {
			result <- c(result, 1)
		}
		else {
			result <- c(result, 0)
		}
		result
})


# calCor: Calculates the Pearsson correlation for every homoeolog / duplicate along the entire set of experiments
setGeneric("calCor", function(object, raac) standardGeneric("calCor"))
setMethod("calCor", signature("pairList", "character"),
	function(object, raac) {
		result <- c()
		genes <- c()
		cat("\nCalculating pearson correlation between homoeologs/duplicates along experiments.\n")
		min <- 0
		max <- length(object@Data)
		pb <- txtProgressBar(min=min, max=max, style=3)
		for (item in object@Data) {
			c <- calCor(item, raac)
			result <- rbind(result, c)
			genes <- c(genes, item@geneName)
			min <- min + 1
			setTxtProgressBar(pb, min)
		}
		rownames(result) <- genes
		colnames(result) <- c("Pearson R", "testDup")
		cat("\n")
		result
})

setMethod("calCor", signature("pair", "character"),
	function(object, raac){
		result <- cor(as.numeric(object@df[1,]), as.numeric(object@df[2,]))
		if (object@geneName %in% raac) {
			result <- c(result, 1)
		}
		else {
			result <- c(result, 0)
		}
		result
})


# plotWTPV: plots the wilcoxon test p-value
setGeneric("plotWTPV", function(table) standardGeneric("plotWTPV"))
setMethod("plotWTPV", signature("matrix"),
	function(table) {
	if (is.null(colnames(table))) {
		return("Data set has no colnames")
	}
	pvals_l <- c()
	pvals_g <- c()
	raacg <- table[table[,"testDup"] == 1,]
	dups <- table[table[,"testDup"] == 0,]
	for(i in 1:(length(colnames(table))-1)) {
		pval_l <- wilcox.test(as.numeric(raacg[,i]), as.numeric(dups[,i]), exact=F, alternative="less")$p.value
		pval_g <- wilcox.test(as.numeric(raacg[,i]), as.numeric(dups[,i]), exact=F, alternative="greater")$p.value
		pvals_l <- c(pvals_l, pval_l)
		pvals_g <- c(pvals_g, pval_g)
	}
	pvals_l[pvals_l == Inf] <- 0
	pvals_g[pvals_g == Inf] <- 0
	pvals_l <- log10(pvals_l)
	pvals_g <- -log10(pvals_g)
	pvals <- rbind(pvals_g, pvals_l)
	colnames(pvals) <- colnames(table)[-length(colnames(table))]
	rownames(pvals) <- c("Greater", "Less")
	cat("\nPlotting Wilcoxon test p values\n")
	print(pvals)
	pdf("WTPV.pdf", width=17, height=7)
	barplot(pvals, xlab="", ylab="log10(Wilcoxon test pvalue)", ylim=c(min(pvals_l), max(pvals_g)), beside=TRUE, col=c("green", "red"))
	dev.off()
})

# plotDif: plots differences between homoeologs raac and not raac along experiments
setGeneric("plotDif", function(table) standardGeneric("plotDif"))
setMethod("plotDif", signature("matrix"),
	function(table) {
		cat("\nAre you honorable enough to plot those results [y|n]? ")
		answer <- readLines()
		cat("\nI do not think you are honorable!\nSee you next time!\n")
})

# plotCor: plots Pearson correlation between raac and no-raac genes along experiments
setGeneric("plotCor", function(table) standardGeneric("plotCor"))
setMethod("plotCor", signature("matrix"),
	function(table) {
		pdf("RAAC_pearson.pdf", width=7, height=7)
		raacg <- abs(as.numeric(as.character(table[table[,"testDup"] == 1,1])))
		dups <- abs(as.numeric(as.character(table[table[,"testDup"] == 0,1])))
		w <- wilcox.test(raacg, dups, alternative="greater")
		boxplot(raacg, dups, names=c("RAAC genes", "No-RAAC genes"), ylab="abs(Pearson correlation Coefficient)", xlab="")
		text(1.5, 0.2, paste("W: ", round(w$p.value, 4), sep=""))
		dev.off()
})

setGeneric("dist.logfc", function(DGE) standardGeneric("dist.logfc"))
setMethod("dist.logfc", signature("data.frame"),
	function(DGE) {
		H <- list()
		S <- list()
		for (name in rownames(DGE)) {
			name <- gsub("\\.[LSPX]$", "", name)
			tmp.data <- DGE[gsub("\\.[LSPX]$", "", rownames(DGE)) %in% name,1]
			if (length(tmp.data) == 2) {
				H[[name]] <- as.numeric(as.character(dist(tmp.data)))
			}
			else {
				S[[name]] <- tmp.data
			}
		}
		SimpleList(homeologs=unlist(H), singletons=unlist(S))
})

setGeneric("test.wilcoxon", function(dist_obj, targets, cutoff) standardGeneric("test.wilcoxon"))
setMethod("test.wilcoxon", signature("SimpleList", "SimpleList", "numeric"),
	function(dist_obj, targets, cutoff) {
		WT <- SimpleList()
		for (name in names(targets)) {
			if (length(dist_obj[[1]][names(dist_obj[[1]]) %in% targets[[name]]]) > 0) {
				wt <- wilcox.test(dist_obj[[1]][names(dist_obj[[1]]) %in% targets[[name]]], dist_obj[[1]][!names(dist_obj[[1]]) %in% targets[[name]]], alternative="less")$p.value
				if (wt <= cutoff) {
					WT[[name]] <- wt
				}
			}
		}
		WT <- WT[order(unlist(WT@listData))]
		WT
})

setGeneric("create.targetGene.matrix", function(path, total_genes, targs) standardGeneric("create.targetGene.matrix"))
setMethod("create.targetGene.matrix", signature("character", "character", "SimpleList"),
	function(path, total_genes, targs) {
		colnames <- names(targs)
		total_genes <- gsub("\\|.*$", "", total_genes)
		mat <- matrix(0, ncol=length(colnames), nrow=length(total_genes))
		rownames(mat) <- total_genes
		colnames(mat) <- colnames
		pb <- txtProgressBar(min=0, max=length(names(targs)), style=3, label="Progress in colnames", char="*")
		counter <- 0
		for (name in names(targs)) {
			counter <- counter + 1
			setTxtProgressBar(pb, counter, label="Progress in colnames")
			for(gname in targs[[name]]) {
				mat[gname, name] <- 1
			}
		}
		mat
})
