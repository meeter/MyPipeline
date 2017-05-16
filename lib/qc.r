###################################################
##Chunk number 1: load the library
###################################################
#load("pca.RData")
Install_And_Load <- function(Required_Packages)
{
    Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
    source("http://bioconductor.org/biocLite.R")
    if(length(Remaining_Packages))
    {
        biocLite(Remaining_Packages);
    }
    Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
    if(length(Remaining_Packages))
    {
        install.packages(Remaining_Packages);
    }
    for(package_name in Required_Packages)
    {
        suppressMessages(library(package_name,character.only=TRUE,quietly=TRUE));
    }
}
Required_Packages=c("ggplot2", "DESeq2", "subSeq", "edgeR", "reshape2")
Install_And_Load(Required_Packages)


###################################################
### code chunk number 2: read parameter
###################################################
args <- commandArgs(trailingOnly = TRUE)

input=args[1]
output=args[2]
VERBOSE=args[3]

if (!is.na(VERBOSE)) {
        print (paste("Input file:",args[1]))
        print (paste("Output:",args[2])) 
}
exp_info <- read.table(input, comment.char = "#",header=T)
#print("Experiment Info Table Loaded")

###################################################
### code chunk number 3: read file
###################################################
ExonCount <- list()
IntronCount <- list()
AntiIntronCount <- list()
for (i in 1:nrow(exp_info)) 
{
	 ExonCountFileName <- paste(output, paste("count", gsub("bam", "ExonCount.bed",basename(as.character(exp_info[i,1]))), sep="/"), sep="/")   
	 ExonCount[[i]] <- read.table(ExonCountFileName, sep="\t", header=F)
	 IntronCountFileName <- paste(output, paste("count", gsub("bam", "IntronCount.bed",basename(as.character(exp_info[i,1]))), sep="/"), sep="/")
         IntronCount[[i]] <- read.table(IntronCountFileName, sep="\t", header=F)
	 AntiIntronCountFileName <- paste(output, paste("count", gsub("bam", "AntiIntronCount.bed",basename(as.character(exp_info[i,1]))), sep="/"), sep="/")
         AntiIntronCount[[i]] <- read.table(AntiIntronCountFileName, sep="\t", header=F)
}

###################################################
### code chunk number 4: Summarize reads
###################################################
Count <- ExonCount
for (i in 1:nrow(exp_info))
{
	Count[[i]][,7] <- ExonCount[[i]][,7] + IntronCount[[i]][,7] + AntiIntronCount[[i]][,7]
}
Count <- do.call(cbind.data.frame, Count)
Count <- Count[,c(1:4,6,seq(7,7*nrow(exp_info), by=7))]
Label <- as.character(exp_info[,3])
print("Reading Count Files Finished")
###################################################
### code chunk number 5: PCA Plot
###################################################
print ("Generating PCA Plot")
cell_type <- factor(exp_info[,4])
colData <- data.frame(cell_type=cell_type) 
row.names(colData) <- Label 
batch.dds <- DESeqDataSetFromMatrix(countData = as.matrix(Count[,6:ncol(Count)]), colData = colData, design = ~ cell_type)
batch.vsd <- varianceStabilizingTransformation(batch.dds)
batch.vsd.pca <- plotPCA(batch.vsd, intgroup=c("cell_type"),returnData=TRUE)

pdf(paste(output,"pca.pdf", sep="/"),height=10, width=15)
percentVar <- round(100 * attr(batch.vsd.pca, "percentVar"))
ggplot(batch.vsd.pca, aes(PC1, PC2, color=cell_type)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	geom_text(aes(label = Label),hjust=0.25, vjust=-0.5, size=5, show_guide = F) + 
	ggtitle("PCA Plot") +
	theme(plot.title=element_text(size=20), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
	title=element_text(size=20, face="bold"), legend.text=element_text(size=10), axis.line = element_line(colour = "black"))
dev.off()

###################################################
### code chunk number 6: Saturation Plot
###################################################
print ("Generating Saturation Plot")
proportion=seq(0.01,1,0.03)
nonhit <- lapply(proportion, function(x) {
	temp <- generateSubsampledMatrix(as.matrix(Count[, 6:ncol(Count)]), x, seed=123)
	return(colSums(temp!=0)/nrow(temp)*100)
	}
)
nonhit <- do.call(rbind.data.frame, nonhit)
nonhit <- cbind(proportion, nonhit)
colnames(nonhit) <- c("Proportion", Label) 
nonhit_long <- melt(nonhit, id="Proportion")  # convert to long format
nonhit_long <- merge(nonhit_long, exp_info[,3:4], all.x=T, by.x="variable", by.y="SampleName")
pdf(paste(output, "saturation.pdf", sep="/"), height=10, width=15)
ggplot(data=nonhit_long, aes(x=Proportion, y=value, colour = variable, linetype=SampleType)) +
	geom_line(size=2) +
	xlab("Proportion of Reads") + 
	ylab("Hit Genes %") +
	ggtitle("Saturation Plot") +  
	theme(plot.title=element_text(size=20), axis.title.y=element_text(size=20), axis.title.x=element_text(size=20),
	title=element_text(size=20, face="bold"), legend.text=element_text(size=20), axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
	 axis.line = element_line(colour = "black"))
dev.off()

##################################################
## code chunk number 7: Manhatton Plot
##################################################
print("Generating Manhattan Plot")
Count.rpkm <- rpkm(Count[,6:ncol(Count)], gene.length=abs(Count[,2] - Count[,3]), normalized.lib.sizes=F, log=T)
#Count.tpm <- Count.rpkm
#for (i in 1:ncol(Count.tpm) ) {
#       Count.tpm[,i] <- Count.rpkm[,i]/apply(Count.rpkm, 2, sum)[i]*10^6
#}
lapply(1:ncol(Count.rpkm), function(i) {
	man.rpkm <- cbind(Count[,4], Count[,1:2], Count.rpkm[,i])
	colnames(man.rpkm) <- c("SNP","CHR","BP","P")
	#man.rpkm$P <- 1/(10^man.rpkm$P)
	#man.rpkm$P <- log(man.rpkm$P + 1, 2)
	man.rpkm$CHR <- gsub("chr", "", man.rpkm$CHR)
	pdf(paste(output, paste(Label[i],".pdf",sep=""), sep="/"), heigh=12, width=16)
	print(ggplot(man.rpkm, aes(x=CHR,y=P,col=CHR)) + 
	geom_point(position = position_jitter(w = 0.4, h = 0.1)) + 
	geom_abline(intercept = median(man.rpkm$P), slope = 0) + 
	ggtitle("") +
	ylim(-0.1,max(man.rpkm$P)) + 
	ylab("Log2-RPKM") + 
	theme(plot.title=element_text(size=20), axis.title.y=element_text(size=20), axis.title.x=element_text(size=20),legend.position="none",
	      axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),		
              title=element_text(size=20, face="bold"), axis.line = element_line(colour = "black")))
	dev.off();
	}
)


##########################################################
### code chunk number 7: Gini Plot (not used) Lorenz Curve
##########################################################
#print ("Generating Lorenz Curve")
#Count.rpkm <- rpkm(Count[,6:ncol(Count)], gene.length=abs(Count[,2] - Count[,3]), normalized.lib.sizes=T, log=F)
#Count.tpm <- Count.rpkm
#for (i in 1:ncol(Count.tpm) ) {
#	Count.tpm[,i] <- Count.rpkm[,i]/apply(Count.rpkm, 2, sum)[i]*10^6
#}
#ineq <- 1:4
#for (i in 1:ncol(Count.tpm) ) {
#	ineq[i] <- ineq(Count.tpm[,i])
#}
#pdf(paste(output, "gini.pdf", sep="/"), height=10, width=15)
#plot(1:length(ineq), ineq)
#text(1:length(ineq), ineq+0.1,labels=Label)
#dev.off()
#Count.Lorenz <- Count.tpm
#for (i in 1:ncol(Count.tpm)) {
#	Count.tpm <- Count.tpm[order(Count.tpm[,i], decreasing=T),]
#	Count.Lorenz[,i] <- cumsum(Count.tpm[,i]/sum(Count.tpm[,i])) * 100
#}
#Count.Lorenz.p <- (1:nrow(Count.tpm))/nrow(Count.tpm) * 100
#
#lorenz <- data.frame(cbind(Count.Lorenz.p, Count.Lorenz))
#colnames(lorenz) <- c("Proportion", Label)
#lorenz_long <- melt(lorenz, id="Proportion")  # convert to long format
#lorenz_long <- merge(lorenz_long, exp_info[,3:4], all.x=T, by.x="variable", by.y="SampleName")
#pdf(paste(output, "lorenz.pdf", sep="/"), height=10, width=15)
#ggplot(data=lorenz_long, aes(x=Proportion, y=value, colour = variable, linetype=SampleType)) +
#        geom_line(size=2) +
#	geom_line(aes(x = Proportion, y = Proportion), color=1) +
#        xlab("Proportion of Genes %") +
#        ylab("Proportion of Reads %") +
#        ggtitle("Lorenz Curve") +
#        theme(plot.title=element_text(size=20), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
#        title=element_text(size=20, face="bold"), legend.text=element_text(size=10),axis.line = element_line(colour = "black"))
#dev.off()

save.image(paste(output, "qc.RData", sep="/"))
q("no")
