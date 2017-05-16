##Test with ctrl
##invoked by VISITs_gene
###################################################
##Chunk number 1: load the library
###################################################
#cd /cclab_nas/jiyu/mypipeline/test
#load("temp_Asun/VISITs_gene.RData") 

Install_And_Load <- function(Required_Packages)
{
    Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
    source("http://bioconductor.org/biocLite.R")
    if(length(Remaining_Packages)) 
    {
        biocLite(Remaining_Packages)
    }
    for(package_name in Required_Packages)
    {
        suppressMessages(suppressWarnings(library(package_name,character.only=TRUE,quietly=TRUE)));
    }
}
Required_Packages=c("genefilter", "gridExtra", "metap", "knitrBootstrap","rmarkdown", "vsn", "edgeR", "GenomeGraphs", "DESeq2",
		    "ggplot2", "easyGgplot2", "DSS", "biomaRt", "BiocParallel");
Install_And_Load(Required_Packages);

###################################################
### code chunk number 2: read parameter
###################################################
#args <- c("temp_Jae_m_d0/exp_info.txt", "1", "temp_Jae_m_d0", "ture", "20","20", "","2500","analyze")
#args <- c("testdata/output/exp_info.txt", "1", "testdata/output", "ture", "20","20","","2500","analyze","DESeq")
args <- commandArgs(trailingOnly = TRUE)
INPUT=args[1]
READS=as.numeric(args[2])
OUTPUT=args[3]
VERBOSE=args[4]
THREADS=as.numeric(args[5])
NUMBER=as.numeric(args[6])
INFO=args[7]
INCLUSION=as.numeric(args[8])
MODE=args[9]
NORM=args[10]
SCR=ifelse(args[11]=="Pos", "greater", "lesser")
if (VERBOSE) {
	#print (paste("ChromInfo:",args[5]))
	message (paste("Normalization Method:" , NORM))
	message (paste("RUNNING MODE:", MODE))
	message (paste("SCREENING:", SCR))
	message (paste("INPUT:", INPUT))
	message (paste("READS:", READS))
	message (paste("Output Dir:", OUTPUT))
	message (paste("PWD:",getwd()))
	message (paste("VERBOSE:", VERBOSE))
	message (paste("SELECTED FILE NAME:", INPUT))
	message (paste("THREADS:", THREADS))
	message (paste("NUMBER of GENES:", NUMBER))
	message (paste("Running Info:", INFO))
	message (paste("Including:", INCLUSION))
}

###################################################
### code chunk number 10: Generate Report
###################################################
if (MODE == "report") {
	rmarkdown::render(input=paste(OUTPUT, "report.rmd", sep="/"), 'knitrBootstrap::bootstrap_document',  clean=T, params = list(number = NUMBER, info = INFO, inclusion = INCLUSION, threads = THREADS))
	q()
}


###################################################
### code chunk number 3: read file
##################################################
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Reading File..."$NORMAL"')
name <- c("chr","start","end","id","exon_num","strand","count")
exp_info <- read.table(INPUT, header=T, comment.char = "#")
if (!is.na(VERBOSE)) {
	                system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Experiment Info Loaded Successfully"$NORMAL"')
}
count <- bplapply(as.character(exp_info[,1]), read.table, header = F, BPPARAM = MulticoreParam(workers = THREADS))
if (!is.na(VERBOSE)) {
	        system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Data Loaded Successfully"$NORMAL"')
}
count <- lapply(1:length(count), function(x) { count[[x]][order(count[[x]]$V4),]})
count.dm <- do.call(cbind.data.frame, count)
colnames(count.dm) <- paste(rep(exp_info[,3], each=7), rep(name, nrow(exp_info)), sep="_")
count.dm[,4] <- sapply(strsplit(as.character(count.dm[,4]), split="\\."), "[[", 1)

#ENV
ID2Symbol_Mus <- read.table("annotation/ID2Symbol_Mus.txt", header=T, sep="\t")  ##GeneID2Symbol Mapping for Mouse
ID2Symbol_Hm  <- read.table("annotation/ID2Symbol_Hm.txt", header=T, sep="\t") #GeneID2Symbol Mapping for Human 
TransID2Symbol_Mus <- read.table("annotation/TransID2Symbol_Mus.txt", header=T, sep="\t") #TranscriptID2Symbol Mapping for Mouse
TransID2Symbol_Hm <- read.table("annotation/TransID2Symbol_Hm.txt", header=T, sep="\t") #TranscriptID2Symbol Mapping for Human
if (!is.na(VERBOSE)) {
	system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Annotation Loaded Successfully"$NORMAL"')
}

###################################################
### code chunk number 4: Get Count
##################################################
####Get Read Count for all (Effective + nonEffective)
all <- count.dm[,c(1,2,3,4,6)]; colnames(all) <- name[c(1,2,3,4,6)];
for (i in 1:(nrow(exp_info)/3)) {
        all[, 5 + i] <- count.dm[, i*21 - 14] + count.dm[, i*21 - 7] + count.dm[, i*21]                             
}
colnames(all)[6:ncol(all)] <- as.character(unique(exp_info[,3]))
####Get Read Count for Effective Insertion (Exon+SenseIntron)
Eff <- count.dm[,c(1,2,3,4,6)]; colnames(Eff) <- name[c(1,2,3,4,6)];
for (i in 1 : (nrow(exp_info)/3)) {
	Eff[, 5 + i] <- count.dm[, i*21 - 7] + count.dm[, i*21] 
}
colnames(Eff)[6:ncol(Eff)] <- as.character(unique(exp_info[,3]))
####Get Read Count for non-Effective Insertion (only AntiIntron)
nonEff <- count.dm[,c(1,2,3,4,6)]; colnames(nonEff) <- name[c(1,2,3,4,6)];
for (i in 1 : (nrow(exp_info)/3)) {
        nonEff[, 5 + i] <- count.dm[, i*21 - 14]
}
colnames(nonEff)[6:ncol(nonEff)] <- as.character(unique(exp_info[,3]))
if (!is.na(VERBOSE)) {
	        system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Count Table Generated"$NORMAL"')
}
#ID2Length_Mus <- read.table("annotation/ID2Length_Mus.txt", header=T, sep="\t", comment.char="#")
#ID2Length_Hm <- read.table("annotation/ID2Length_Hm.txt", header=T, sep="\t", comment.char="#")
#ID2Length_Hm[,1] <- sapply(strsplit(as.character(ID2Length_Hm[,1]), split="\\."), "[[", 1)
#ID2Length_Hm <- merge(ID2Length_Hm, Eff, all.y=T, by.x="Geneid", by.y="id")
#ID2Length_Hm$EffLength <- as.numeric(ID2Length_Hm[,2]) + as.numeric(ID2Length_Hm[,5]) - as.numeric(ID2Length_Hm[,4])
#
#Eff.rpkm <- rpkm(Eff[,6:ncol(Eff)], gene.length=ID2Length_Hm[,9]); Eff.rpkm <- cbind(Eff[,1:5], Eff.rpkm)
#Eff.tpm <- sapply(6:ncol(Eff.rpkm), function(x) {Eff.rpkm[,x]/sum(Eff.rpkm[,x])*10^6}); Eff.tpm <- cbind(Eff[,1:5], Eff.tpm)
#all.rpkm <- rpkm(all[, 6:ncol(all)], gene.length=ID2Length_Hm[,5] - ID2Length_Hm[,4]); all.rpkm <- cbind(all[,1:5], all.rpkm)
#all.tpm <- sapply(6:ncol(all.rpkm), function(x) {all.rpkm[,x]/sum(all.rpkm[,x])*10^6}); all.tpm <- cbind(all[,1:5], all.tpm)

###################################################
### code chunk number 5: Normalization Methods
###################################################
Eff.data <- Eff[,6:ncol(Eff)]
nonEff.data <- nonEff[,6:ncol(nonEff)]
###Filtering
if (nrow(exp_info) <= 6) {
		isexpr.sense <-  Eff[,6] > READS & nonEff[, 7] > READS
} else {
		isexpr.sense <-  rowSums(Eff[, seq(6, ncol(Eff)-1,2)] > 1) >= READS & rowSums(nonEff[,seq(7, ncol(nonEff),2)] > 1) >= READS
}
filter <- rowSums(Eff.data > READS) >= READS
###Different Methods
if(NORM == "TC") {
#norm.fact <- colSums(all[filter,6:ncol(all)])/10^6 
norm.fact <- colSums(Eff[filter,6:ncol(Eff)])/mean(colSums(Eff[filter,6:ncol(Eff)]))
norm.fact.sense <- colSums(all[isexpr.sense, 6:ncol(all)])/mean(colSums(all[isexpr.sense,6:ncol(all)]))
} else if(NORM == "DESeq") {
norm.fact <- estimateSizeFactorsForMatrix(Eff[filter,6:ncol(Eff)])
norm.fact.sense <- estimateSizeFactorsForMatrix(all[isexpr.sense, 6:ncol(all)])
} else if(NORM == "TMM") {
norm.fact <- calcNormFactors(Eff[filter,6:ncol(Eff)]) * colSums(Eff[filter,6:ncol(Eff)])/mean(colSums(Eff[filter,6:ncol(Eff)]))
norm.fact.sense <- calcNormFactors(all[isexpr.sense,6:ncol(all)]) * colSums(all[isexpr.sense,6:ncol(all)])/mean(colSums(all[isexpr.sense, 6:ncol(all)])) 
} else if(NORM == "upper") {
	norm.fact <- calcNormFactors(Eff[filter,6:ncol(Eff)], method="upperquartile")
	norm.fact.sense <- calcNormFactors(all[isexpr.sense,6:ncol(all)], method="upperquartile")
} else if (NORM == "ctrl") {
	controlGenes <- rowSums(sapply(6:ncol(Eff), function(x) {
			Eff.tpm[filter, x] <= quantile(Eff.tpm[filter,x], probs=0.95) & Eff.tpm[filter, x] >= quantile(Eff.tpm[filter,x], probs=0.5)})) == length(6:ncol(Eff))
	norm.fact <- estimateSizeFactorsForMatrix(Eff[filter,6:ncol(Eff)], controlGenes = controlGenes)	
        controlGenes <- rowSums(sapply(6:ncol(all), function(x) {
			all.tpm[isexpr.sense, x] <= quantile(all.tpm[isexpr.sense,x], probs=0.95) & all.tpm[isexpr.sense, x] >= quantile(all.tpm[isexpr.sense,x], probs=0.5)})) == length(6:ncol(all))
	norm.fact.sense <- estimateSizeFactorsForMatrix(all[isexpr.sense, 6:ncol(all)], controlGenes = controlGenes)

} else if (NORM == "CGAT") {
	norm.fact <- sapply(6:ncol(Eff), function(x) { sum(Eff[filter,][Eff.tpm[filter,x] <= quantile(Eff.tpm[filter,x], probs=0.95),x])
	})
	norm.fact <- norm.fact / mean(norm.fact)
	norm.fact.sense <- sapply(6:ncol(all), function(x) { sum(all[isexpr.sense,][all.tpm[isexpr.sense,x] <= quantile(all.tpm[isexpr.sense, x], probs=0.95),x])
				         })
	norm.fact.sense <- norm.fact.sense / mean(norm.fact.sense)
} else {
norm.fact <- rep(1, nrow(exp_info)/3)
norm.fact.sense <- rep(1, nrow(exp_info)/3)
}
###Calculate Normalized Value
Eff.data.norm  <- sapply(1:ncol(Eff.data), function(x) {Eff.data[filter,x] / norm.fact[x]})
colnames(Eff.data.norm) <- colnames(Eff[,6:ncol(Eff)])
nonEff.data.norm <- sapply(1:ncol(nonEff.data), function(x) {nonEff.data[filter,x] / norm.fact[x]})
colnames(nonEff.data.norm) <- colnames(Eff.data)
Eff.prop.norm <- sapply(1:ncol(Eff.data), function(x) {Eff.data[filter,x] / ( Eff.data[filter,x] + nonEff.data[filter,x])})
Eff.prop.norm <- replace(Eff.prop.norm, is.na(Eff.prop.norm), 0)
colnames(Eff.prop.norm) <- colnames(Eff.data)
if (!is.na(VERBOSE)) {
	        system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Normalization Done"$NORMAL"')
}

#############################################################
### code chunk number 6: Statistical Test without replicates
############################################################
if (MODE == "analyze") {
if (nrow(exp_info) <= 6) {
####Fisher Test
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Performing test for count enrichment without replicates..."$NORMAL"')
#temp <- data.filtered[,c("selected.InactCount","ctrl.InactCount","selected.rest","ctrl.rest")]
if (grep("Selected", exp_info$SampleType, ignore.case = T)[1] == 1 ) {##Selected Sample ranked first
	temp <- cbind(Eff.data.norm[,1], Eff.data.norm[,2], sum(Eff.data.norm[,1]) - Eff.data.norm[,1], sum(Eff.data.norm[,2]) - Eff.data.norm[,2]) 
}
if (grep("Selected", exp_info$SampleType, ignore.case = T)[1] == 4 ) {##Ctrl Sample ranked first
        temp <- cbind(Eff.data.norm[,2], Eff.data.norm[,1], sum(Eff.data.norm[,2]) - Eff.data.norm[,2], sum(Eff.data.norm[,1]) - Eff.data.norm[,1])
}
temp <- sapply(1:ncol(temp), function(x) as.integer(temp[,x]))
fisher.res <- bplapply(1:nrow(temp), function(x) fisher.test(matrix(as.numeric(temp[x,]), ncol=2),alternative = "greater"), 
			BPPARAM = MulticoreParam(workers = THREADS))
data.filtered <- Eff[filter,1:5]
data.filtered$p <- unlist(sapply(fisher.res, "[[", "p.value"))
data.filtered$fc <- sapply(fisher.res, "[[", "estimate")
data.filtered$fdr <- p.adjust(data.filtered$p,method="BH")
###Independent Filtering
#Ind_Filter = filtered_R(alpha=0.01, filter=rowMeans(all[filter,6:7]), test=data.filtered$p, theta=seq(from=0, to=1, by=0.01), method="BH")
#Filter_cut = as.numeric(gsub("%", "", names(which(Ind_Filter == max(Ind_Filter))[1]))) / 100
#data.filtered <- data.filtered[rowMeans(all[filter, 6:7]) > quantile(rowMeans(all[filter, 6:7]), Filter_cut), ]
#data.filtered$fdr <- p.adjust(data.filtered$p, method="BH")

###Binomial Test
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Performing test for sense enrichment without replicates..."$NORMAL"')
if (grep("Selected", exp_info$SampleType, ignore.case=T)[1] < grep("Ctrl", exp_info$SampleType, ignore.case=T)[1]) {##Selected Sample ranked first {
	bino.res <- bplapply(1:nrow(Eff.data[isexpr.sense,]), function(x) {
		binom.test(as.integer(Eff.data[isexpr.sense,][x, 1]/norm.fact.sense[1]), n=as.integer(all[isexpr.sense,][x,6]/norm.fact.sense[1]), p=Eff.data[isexpr.sense,][x,2]/all[isexpr.sense,][x,7], alternative=SCR)
		}, 
	BPPARAM = MulticoreParam(workers = THREADS))
}
if (grep("Selected", exp_info$SampleType, ignore.case=T)[1] > grep("Ctrl", exp_info$SampleType, ignore.case=T)[1]) {##Ctrl Sample ranked first {
        bino.res <- bplapply(1:nrow(Eff.data[isexpr.sense,]), function(x) {
                binom.test(as.integer(Eff.data[isexpr.sense,][x, 2]/norm.fact.sense[2]), n=as.integer(all[isexpr.sense,][x,7]/norm.fact.sense[2]), p=Eff.data[isexpr.sense,][x,1]/all[isexpr.sense,][x,6], alternative=SCR)
                },
        BPPARAM = MulticoreParam(workers = THREADS))
}
for (i in 1:length(bino.res)) { ##Remove null value caused by 0 read after normalization
	        if (is.null(bino.res[[i]]$estimate)) bino.res[[i]]$estimate <- 0
        if (is.null(bino.res[[i]]$null.value)) bino.res[[i]]$null.value <- 0
	        if (is.null(bino.res[[i]]$p.value)) bino.res[[i]]$p.value <- 1
}
p2 <- sapply(bino.res, "[[", "p.value")
#enrich <- ifelse((sapply(bino.res, "[[", "estimate") - sapply(bino.res, "[[", "null.value") >= 0), 1, -1)
enrich <- sapply(bino.res, "[[", "estimate") / sapply(bino.res, "[[", "null.value")
fdr2 <- p.adjust(p2, method="BH")
bino.res.dm <- data.frame(ID=Eff[isexpr.sense,]$id, p2=p2, enrich=enrich, fdr2=fdr2)
###Independent Filtering
#Ind_Filter = filtered_R(alpha=0.01, filter=rowMeans(all[isexpr.sense,6:7]), test=bino.res.dm$p2, theta=seq(from=0, to=1, by=0.01), method="BH")
#Filter_cut = as.numeric(gsub("%", "", names(which(Ind_Filter == max(Ind_Filter))[1]))) / 100
#bino.res.dm <- bino.res.dm[rowMeans(all[isexpr.sense, 6:7]) > quantile(rowMeans(all[isexpr.sense, 6:7]), Filter_cut), ]
#bino.res.dm$fdr2 <- p.adjust(bino.res.dm$p2, method="BH")

########Merge Results
data.filtered <- merge(data.filtered, bino.res.dm, by.x="id", by.y="ID")  
if (length(grep("ENSMUSG", data.filtered$id[1])))  {
	data.filtered <- merge(data.filtered, ID2Symbol_Mus, all.x=T, by.x="id", by.y=1)
} else if (length(grep("ENSG", data.filtered$id[1]))) {
	data.filtered <- merge(data.filtered, ID2Symbol_Hm, all.x=T, by.x="id", by.y=1)
} else if (length(grep("ENSMUST", data.filtered$id[1]))) {
        data.filtered <- merge(data.filtered, TransID2Symbol_Hm, all.x=T, by.x="id", by.y=1)
} else if (length(grep("ENST", data.filtered$id[1]))) {
        data.filtered <- merge(data.filtered, TransID2Symbol_Hm, all.x=T, by.x="id", by.y=1)
}
data.filtered.cpm <- cpm(cbind(Eff[,6:7], nonEff[,6:7]), normalized.lib.sizes=F, log=F); row.names(data.filtered.cpm) <- Eff$id
data.filtered <- merge(data.filtered, data.filtered.cpm, all.x=T, by.x="id", by.y=0)
colnames(data.filtered)[13:14] <- paste(colnames(data.filtered)[13:14], ".Eff", sep="")
colnames(data.filtered)[15:16] <- paste(colnames(data.filtered)[15:16], ".nonEff", sep="")
colnames(data.filtered)[12] <- "Gene_Symbol"
for (i in 1:nrow(data.filtered)) {
	if (is.na(data.filtered[i,]$p2)) { ###if sense P is NA, use count P instead
		data.filtered$sumz_p[i] <- data.filtered$p[i]
	} else if (data.filtered[i,]$p ==1 | data.filtered[i,]$p2 ==1 ) { ##P in metaP has to be 0-1
	data.filtered$sumz_p[i] <- 1
	} else if (data.filtered[i,]$p < 1e-30 & data.filtered[i,]$p2 < 1e-30) { ##both p very small
	data.filtered$sumz_p[i] <-0
	} else {	
	data.filtered$sumz_p[i] <- sumz(
		c(data.filtered[i,]$p+1e-30, data.filtered[i,]$p2+1e-30), c(rowMeans(data.filtered[i,13:14]), rowMeans(data.filtered[i,13:16]))
		)$p
	}
}
data.filtered$sumz_fdr <- p.adjust(data.filtered$sumz_p, method='BH')

system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Test done..."$NORMAL"')
}
#############################################################
### code chunk number 7: Statistical Test with replicates
############################################################
if (nrow(exp_info) > 6) {
############edgeR+DSS
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Performing test for count enrichment with replicates..."$NORMAL"')
SampleType=factor(exp_info$SampleType)
colData <- unique(exp_info[,4:ncol(exp_info)])
colData <- data.frame(sapply(1:ncol(colData), function(x) {factor(colData[,x])}))
colnames(colData) <- colnames(exp_info)[4:ncol(exp_info)]
Eff.filtered <- Eff[filter,]
row.names(colData) <- colnames(Eff.filtered[,6:ncol(Eff.filtered)])
row.names(Eff.filtered) <- Eff.filtered[,4]
#design <- as.formula(paste("~", paste(colnames(colData)[ncol(colData):1],collapse="+"))) 
#design.matrix <- model.matrix(design, data=colData)
#seqData=newSeqCountSet(as.matrix(Eff.filtered[,6:ncol(Eff.filtered)]), as.data.frame(design.matrix))
#normalizationFactor(seqData) <- norm.fact
#seqData=estDispersion(seqData)
#edgeR.DSS.fit <- glmQLFit(as.matrix(Eff.filtered[,6:ncol(Eff.filtered)]), design.matrix,robust=T,
#			  dispersion=dispersion(seqData), lib.size=norm.fact) # * mean(colSums(Eff.filtered[,6:ncol(Eff.filtered)])))
#edgeR.DSS <- glmQLFTest(glmfit=edgeR.DSS.fit)
#edgeR.DSS$table$FDR <- p.adjust(edgeR.DSS$table$PValue, method='BH')
#edgeR.DSS <- edgeR.DSS$table
#for (i in 1:nrow(edgeR.DSS)) {
#	  if (SCR=="greater" & edgeR.DSS$logFC[i] > 0) {edgeR.DSS$PValue[i] <- edgeR.DSS$PValue[i]/2
#	  	} else if (SCR=="lesser" & edgeR.DSS$logFC[i] < 0 )  {edgeR.DSS$PValue[i] <- edgeR.DSS$PValue[i]/2
#		} else {edgeR.DSS$PValue[i] <- 1 - edgeR.DSS$PValue[i]}
#}
#edgeR.DSS$FDR <- p.adjust(edgeR.DSS$PValue, method="BH")
#Eff.res <- edgeR.DSS
design = as.formula(paste("~", paste(colnames(colData)[ncol(colData):1],collapse="+")))
counts = as.matrix(Eff.filtered[,6:ncol(Eff.filtered)])
DESeq_DSS.data <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = design)
seqData=newSeqCountSet(counts, as.data.frame(model.matrix(design, data = colData)))
#nf <- colSums(counts)/mean(colSums(counts))
nf <- matrix(data=rep(norm.fact, nrow(counts)), byrow=T, ncol=ncol(counts))
normalizationFactors(DESeq_DSS.data) <- nf
seqData=estDispersion(seqData)
DESeq_DSS.data <- estimateDispersions(DESeq_DSS.data, fitType="local") 
dispersions(DESeq_DSS.data) <-  dispersion(seqData)
DESeq_DSS <- nbinomWaldTest(DESeq_DSS.data, betaPrior=F)
Eff.res <- results(DESeq_DSS, lfcThreshold=0,  cooksCutoff=FALSE, altHypothesis=SCR,  independentFiltering=F) 
Eff.res <- Eff.res[order(Eff.res$padj),]

if (length(grep("ENSMUSG", all$id[1]))) { 
	Eff.res <- merge(as.data.frame(Eff.res), ID2Symbol_Mus, all.x=T, by.x=0, by.y=1)
} else if (length(grep("ENSG", all$id[1]))) {
	Eff.res <- merge(as.data.frame(Eff.res), ID2Symbol_Hm, all.x=T, by.x=0, by.y=1)
} else if (length(grep("ENSMUST", all$id[1]))) {
	Eff.res <- merge(as.data.frame(Eff.res), TransID2Symbol_Mus, all.x=T, by.x=0, by.y=1)
} else if (length(grep("ENST", all$id[1]))) {
        Eff.res <- merge(as.data.frame(Eff.res), TransID2Symbol_Hm, all.x=T, by.x=0, by.y=1)
}
Eff.res <- merge(Eff.res, Eff[,1:5], by.x=1, by.y="id", all.x=T)
############DSS
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Performing test for sense enrichment with replicates..."$NORMAL"')
#isexpr.sense <-  rowSums(Eff[,6:(5+nrow(exp_info)/6)] > 1) >=2  & rowSums(nonEff[,(ncol(nonEff)-nrow(exp_info)/6+1):ncol(nonEff)] > 1) >=2 ##reads in both Eff-Sel and nonEff-Ctrl
DSS <- lapply(6:ncol(Eff),function(i) {
	data.frame(chr=all[isexpr.sense,]$chr, pos=as.numeric(gsub("ENST", "", gsub("ENSMUST", "", gsub("ENSG", "", gsub("ENSMUSG","", all[isexpr.sense,]$id))))), 
		   N = as.integer(all[isexpr.sense,i]/norm.fact.sense[i-5]), X = as.integer(Eff[isexpr.sense,i]/norm.fact.sense[i-5])) #chr, pos, coverage, Eff
})
all_DSS <- makeBSseqData(DSS, colnames(all)[6:ncol(all)])
temp <- sapply(data.frame(assays(all_DSS)$Cov), function(x){replace(x, x == 0,1)})
colnames(temp) <- colnames(assays(all_DSS)$Cov)
assays(all_DSS)$Cov <- temp
all_DSS.fit <- DMLfit.multiFactor(all_DSS, design=colData, formula = as.formula(paste("~", paste(colnames(colData)[1:ncol(colData)],collapse="+"))))             
all_DSS.res = DMLtest.multiFactor(all_DSS.fit, coef=2)
all_DSS.res <- all_DSS.res[order(all_DSS.res$pos),]
all_DSS.res$id <- all[isexpr.sense,]$id 
for (i in 1:nrow(all_DSS.res)) {
	if (SCR=="greater") {all_DSS.res$pvals[i] <- pnorm(-all_DSS.res$stat[i])}
        else {all_DSS.res$pvals[i] <- 1-pnorm(-all_DSS.res$stat[i])}
}

############Independent Filtering, not used
#Ind_Filter = filtered_R(alpha=0.01, filter=rowMeans(cpm(assays(all_DSS)$Cov,log=T)), test=all_DSS.res$pvals, theta=seq(from=0, to=1, by=0.01), method="BH")
#Filter_cut = as.numeric(gsub("%", "", names(which(Ind_Filter == max(Ind_Filter))[1]))) / 100
#all_DSS.res <- all_DSS.res[rowMeans(assays(all_DSS)$Cov) > quantile(rowMeans(cpm(assays(all_DSS)$Cov, log=T)), Filter_cut), ] 
all_DSS.res$fdrs <- p.adjust(all_DSS.res$pvals, method="BH")
#nrow(subset(all_DSS.res, fdrs < 0.05 & stat > 1)) #47

############Merge DESeq2 with DSS and data in each sample
#Eff.res <- merge(Eff.res[,c(1,2,6,7,8,9,10,11)], all_DSS.res[,3:6], all.x=T, by.x=1, by.y="id") #merge edgeR with DSS
data.filtered.cpm <- cpm(cbind(Eff[filter,6:ncol(Eff)], nonEff[filter, 6:ncol(nonEff)]),  normalized.lib.sizes=F, log=F); 
row.names(data.filtered.cpm) <- Eff[filter,]$id
Eff.res <- merge(Eff.res[,c(1,3,6,7,8,9,10,11,12)], all_DSS.res[,3:6], all.x=T, by.x=1, by.y="id")
Eff.res <- merge(Eff.res, data.filtered.cpm, by.x=1, by.y=0)
colnames(Eff.res)[13:(12+nrow(exp_info)/3)] <- paste(colnames(Eff.res)[13:(12+nrow(exp_info)/3)], ".Eff", sep="") #edgeR: 13->11
colnames(Eff.res)[(13+nrow(exp_info)/3):ncol(Eff.res)] <- paste(colnames(Eff.res)[(13+nrow(exp_info)/3):ncol(Eff.res)], ".nonEff", sep="")#edgeR: 11->12
for (i in 1:nrow(Eff.res)) {
	if (is.na(Eff.res[i,]$pvals)) {
		Eff.res$sumz_p[i] <- Eff.res$pvalue[i]
	} else if (Eff.res[i,]$pvalue ==1 | Eff.res[i,]$pvals ==1 ) {
		Eff.res$sumz_p[i] <- 1
	} else if (Eff.res[i,]$pvalue < 1e-30 & Eff.res[i,]$pvals < 1e-30) {
		Eff.res$sumz_p[i] <- 0
	} else {
		Eff.res$sumz_p[i] <- sumz(
					c(Eff.res[i,]$pvalue+1e-30, Eff.res[i,]$pvals+1e-30), c(mean(data.filtered.cpm[i,1:(ncol(data.filtered.cpm)/2)]), mean(data.filtered.cpm[i,])))$p
		}
}

Eff.res$sumz_fdr <- p.adjust(Eff.res$sumz_p, method='BH')

system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Test done..."$NORMAL"')
colnames(Eff.res)[1:12] <- c("ID","Log2FC","p.count","fdr.count","Gene_Symbol","chr","start","end","strand","DSS_stat","p.sense","fdr.sense")
}
###################################################
### code chunk number 8: Output
###################################################
system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Output..."$NORMAL"')
#if (exists("bedfile")) {rm(bedfile)} ##Clean big data
#RData=paste(paste(READS, "withCtrl", sep="_"),"RData",sep=".")
if (nrow(exp_info) <= 6) {
	write.table(data.filtered[,c(-2:-5)], paste(OUTPUT,"results_gene.txt",sep="/"),sep="\t",row.names=F,quote=F)
} else {
	write.table(Eff.res, paste(OUTPUT,"results_gene.txt",sep="/"),sep="\t",row.names=F,quote=F)
}
###Bed File for Visualization
bedfile_Name <- unique(paste(OUTPUT, paste("bed/", paste(exp_info$SampleName, ".bed", sep=""), sep=""),sep="/"))
bedfile <- bplapply(bedfile_Name, read.table, header = F, BPPARAM = MulticoreParam(workers = THREADS))
names(bedfile) <- unique(exp_info$SampleName)
save.image(paste(OUTPUT, "gene.RData", sep="/"))
}
###################################################
### code chunk number 9: Diagnose Plot
###################################################
if (MODE == "diagnose") {
message("Generating Diagnose Plots")
########################MA Plot 
###Read in Proliferation Genes
Prof.hs <- read.table("annotation/Hs_Proliferate.txt", sep="\t", header=T)
Prof.mm <- read.table("annotation/Mm_Proliferate.txt", sep="\t", header=T)
row.names(Eff.data.norm) <- Eff[filter,]$id
###Get M-A value for Proliferation Genes
if (length(grep("ENSMUSG", Eff$id[1]))) {
			Prof.points <- Eff.data.norm[match(Prof.mm[,1], row.names(Eff.data.norm)),]
			Prof.points <- Prof.points[!is.na(Prof.points[,1]),]
} else if (length(grep("ENSG", Eff$id[1]))) {
			Prof.points <- Eff.data.norm[match(Prof.hs[,1], row.names(Eff.data.norm)),]
			Prof.points <- Prof.points[!is.na(Prof.points[,1]),]
} else if (length(grep("ENSMUST", Eff$id[1]))) {
	                #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Mus, by.x=0, by.y=1)
} else if (length(grep("ENST", Eff$id[1]))) {
	                #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Hm, by.x=0, by.y=1)
}
row.names(Eff.prop.norm) <- Eff[filter,]$id
if (length(grep("ENSMUSG", Eff$id[1]))) {
	                Prof.sense.points <- Eff.prop.norm[match(Prof.mm[,1], row.names(Eff.prop.norm)),]
                        Prof.sense.points <- Prof.sense.points[!is.na(Prof.sense.points[,1]),]
} else if (length(grep("ENSG", Eff$id[1]))) {
	                Prof.sense.points <- Eff.prop.norm[match(Prof.hs[,1], row.names(Eff.prop.norm)),]
                        Prof.sense.points <- Prof.sense.points[!is.na(Prof.sense.points[,1]),]
} else if (length(grep("ENSMUST", Eff$id[1]))) {
	                #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Mus, by.x=0, by.y=1)
} else if (length(grep("ENST", Eff$id[1]))) {
	                #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Hm, by.x=0, by.y=1)
}
###Plot
pdf(paste(OUTPUT,"diagnose/MA-Plot.pdf", sep="/"), height=10, width=16)
par(mfrow=c(1,2))
if( nrow(exp_info) > 6) {
	    scatter.smooth(rowMeans(log(Eff.data.norm+1, 2)),
			   rowMeans(log(Eff.data.norm[,1:(ncol(Eff.data.norm)/2)]+1, 2) - log(Eff.data.norm[,(ncol(Eff.data.norm)/2+1):ncol(Eff.data.norm)]+1, 2)),
				    main="M-A Plot for Count Enrichment", xlab="A",ylab="M", lpars=list(col = "blue", lwd = 4)) ##M-A plot
            points(rowMeans(log(Prof.points+1, 2)),
		           rowMeans(log(Prof.points[,1:(ncol(Eff.data)/2)]+1,2) - log(Prof.points[,(ncol(Eff.data)/2+1):ncol(Eff.data)]+1,2)), col=3, lwd=3)
} else {
	    scatter.smooth(rowMeans(log(Eff.data.norm+1,2)), log(Eff.data.norm[,1]+1,2) - log(Eff.data.norm[,2]+1, 2),
			   main="M-A Plot for Count Enrichment", xlab="A",ylab="M", lpars=list(col = "blue", lwd = 4))
            points(rowMeans(log(Prof.points+1,2)), log(Prof.points[,1]+1,2) - log(Prof.points[,2]+1,2), col=3, lwd=3)
}
legend("bottomright", pch=1, bty="n", col=3, "Proliferation")
abline(h=0, lty=3, col=2, lwd=2)
if( nrow(exp_info) > 6) {
	    scatter.smooth(rowMeans(Eff.prop.norm), rowMeans(Eff.prop.norm[,1:(ncol(Eff.prop.norm)/2)] - Eff.prop.norm[,(ncol(Eff.prop.norm)/2+1):ncol(Eff.prop.norm)]),                           main="M-A Plot for Sense Enrichment", xlab="A",ylab="M", lpars=list(col = "blue", lwd = 4)) ##M-A plot
            points(rowMeans(Prof.sense.points),
		   rowMeans(Prof.sense.points[,1:(ncol(Eff.data)/2)] - Prof.sense.points[,(ncol(Eff.data)/2+1):ncol(Eff.data)]),col=3, lwd=3)
} else {
	    scatter.smooth(rowMeans(Eff.prop.norm), Eff.prop.norm[,1] - Eff.prop.norm[,2],
			   main="M-A Plot for Sense Enrichment", xlab="A",ylab="M", lpars=list(col = "blue", lwd = 4))
            points(rowMeans(Prof.sense.points), Prof.sense.points[,1] - Prof.sense.points[,2], col=3, lwd=3)
}
legend("bottomright", pch=1, bty="n", col=3, "Proliferation")
abline(h=0, lty=3, col=2, lwd=3)
dev.off()
#############################Dispersion Plot
d1 <- meanSdPlot(log(Eff.data.norm+1, 2), rank=F)
d2 <- meanSdPlot(Eff.prop.norm,rank=F)
pdf(paste(OUTPUT,"diagnose/Mean-Variance.pdf", sep="/"),width=16)
grid.arrange(d1$gg,d2$gg, ncol=2)
dev.off()
##############################GC
Hs.ID2Length <- read.table("annotation/Hs_ID_EffLength_GC.txt", sep="\t", header=T)
Mm.ID2Length <- read.table("annotation/Mm_ID_EffLength_GC.txt", sep="\t", header=T)
if (length(grep("ENSMUSG", Eff$id[1]))) {
	        Eff.cqn <- merge(Eff[filter,], Mm.ID2Length, by.x="id", by.y="Geneid")
} else if (length(grep("ENSG", Eff$id[1]))) {
	        Eff.cqn <- merge(Eff[filter,], Hs.ID2Length, by.x="id", by.y="Geneid")
} else if (length(grep("ENSMUST", Eff$id[1]))) {
	        #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Mus, by.x=0, by.y=1)
} else if (length(grep("ENST", Eff$id[1]))) {
	        #Eff.cqn <- merge(Eff[filter,], TransID2Symbol_Hm, by.x=0, by.y=1)
}
#cqn.subset <- cqn(Eff.cqn[, 6:(ncol(Eff.cqn)-2)], lengths = Eff.cqn$Eff_Size,
#		  x = Eff.cqn$percentage_gc_content, sizeFactors = NULL, verbose = TRUE)
#cqnRes <- data.frame(cqn.subset$y + cqn.subset$offset)
q <- list()
for (i in 1:length(6:(ncol(Eff.cqn)-2))) {
	name=colnames(Eff.cqn)[i+5]
	pdf(paste(paste(OUTPUT, "diagnose", paste(name,"_GC.pdf",sep=""), sep="/")), width=12, height=8) ##OUTPUT, "diagnose",sep="/"
	print(qplot(Eff.cqn$percentage_gc_content, cpm(Eff.cqn[,i+5], log=T), 
		    xlab="GC%", ylab="Log2-Normalized Count", geom = c("point", "smooth"), main=name))
	dev.off()
}
#pdf(paste(OUTPUT, "diagnose", "GC.pdf", sep="/"), width=12, height=2*length(q))
#print(q[[name]])
#print(marrangeGrob(q, nrow=length(q)/2, ncol=2))
#dev.off()
##########################BCV between ctrl and sel, only for dataset with replicates
###Count Enrichment using BCV; Sense Enrichment using SD
if( nrow(exp_info) > 6) {
getBCV <- function(counts, norm.fact) {
	edgeR.data <- DGEList(counts = counts, genes= row.names(counts))       
	edgeR.data <- calcNormFactors(edgeR.data) 
	edgeR.data$samples$norm.factors <- norm.fact/edgeR.data$samples$lib.size
	edgeR.data <- estimateCommonDisp(edgeR.data)
	edgeR.data <- estimateTrendedDisp(edgeR.data)
	edgeR.data <- estimateTagwiseDisp(edgeR.data) 
	return(sqrt(edgeR.data$tagwise.dispersion))
}
BCV <- data.frame(matrix(data=NA, ncol=4, nrow=table(filter)[2]))
colnames(BCV) <- c("Count_sel","Count_ctrl","Sense_sel","Sense_ctrl")
counts=as.matrix(Eff.data[filter,])
BCV[,1] <- getBCV(counts[, 1:(ncol(counts)/2)], norm.fact[1:(ncol(counts)/2)])
BCV[,2] <- getBCV(counts[, (ncol(counts)/2+1):ncol(counts)], norm.fact[(ncol(counts)/2+1):ncol(counts)])
BCV[,3] <- apply(Eff.prop.norm[, 1:(ncol(counts)/2)], 1, sd)
BCV[,4] <- apply(Eff.prop.norm[, (ncol(counts)/2+1):ncol(counts)], 1, sd)
BCV.long <- data.frame(IntraVariance = c(BCV[,1], BCV[,2], BCV[,3], BCV[,4]))
BCV.long$Type <- sapply(strsplit(rep(colnames(BCV)[1:4], each=nrow(BCV)), split="_"), "[[", 1)
BCV.long$Sample <- sapply(strsplit(rep(colnames(BCV)[1:4], each=nrow(BCV)), split="_"), "[[", 2)

pdf(paste(OUTPUT,"diagnose/BCV.pdf", sep="/"))
print(ggplot2.boxplot(data=BCV.long, xName="Type", yName='IntraVariance', groupName='Sample', 
		position=position_dodge(0.8), #interval between boxplot of the same group
		backgroundColor="white", brewerPalette="Paired",#groupColors=c('#999999','#E69F00','#56B4E9')
		addDot=F)) #,legendPosition="right" 
dev.off() 
}
#pdf(paste(OUTPUT,"diagnose/SENSE_Box-Plot.pdf", sep="/"),width=12)
#boxplot(Eff.prop.norm)
#dev.off()
###Boxplot
#pdf(paste(OUTPUT,"diagnose/COUNT_Box-Plot.pdf", sep="/"),width=12)
#boxplot(log(Eff.data.norm+1,2))
#dev.off()
##############P-value Histogram

}

system('RED="\\033[1;31m";NORMAL="\\033[0;39m";echo "$RED"Done"$NORMAL"')
