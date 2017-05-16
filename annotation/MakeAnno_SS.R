### Use GenomicFeatures to retrieve Ensembl/RefSeq
### Generate synthetic Exon, interval and reduced
### intron, output bed file
###################################################

###################################################
### code chunk number 1: load the library
###################################################
GI <- require(GenomicFeatures,warn.conflicts = F, quietly = T)
IR <- require(IRanges,warn.conflicts = F, quietly = T)
BR <- require(biomaRt,warn.conflicts = F, quietly = T)
BP <- require(BiocParallel,warn.conflicts = F, quietly = T)
if (GI & IR & BR & BP) {print("SUCCESS IN LOADING PACKAGES")}else{stop("ERROR IN LOADING PACKAGES")}
						 
#tx <- loadDb("hg38.sqlite")
###################################################
### code chunk number 2: Get Gene Annotation
###################################################
args <- commandArgs(trailingOnly = TRUE)
#args <- c("hg18","ensGene")
tx <- makeTxDbFromUCSC(genome=args[1],tablename=args[2]) #########Change Gene Build Version Here (hg18/hg19/mm9/mm10,ensGene/refGene)
output <- paste(args[1],args[2],sep="_")
dir.create(output)
NCORE=as.numeric(args[3])
#saveDb(tx, file="hg38.sqlite") 
#columns(tx)
#keytypes(tx)
txname <- elementMetadata((transcripts(tx)))[,2]  

###################################################
### code chunk number 3: Label refGene with Multi-Loci
### If yes, add "_num" to the end of its Gene ID
### Use TXID to distinguish isoform and gene copy
###################################################
if (args[2] == "refGene")
	{refseq2gene <- select(tx, keys = txname, columns=c("TXID","TXCHROM","TXSTRAND","TXSTART","TXEND","TXNAME","GENEID"),keytype="TXNAME") 						
	 refseq2gene <- refseq2gene[order(refseq2gene$GENEID,refseq2gene$TXNAME),];row.names(refseq2gene) <- 1:dim(refseq2gene)[1] 
	 print ("Get RefGene <-> Transcripts Mapping")
	 #####Get genes with multi-loci, these genes have same refseq ID, but in different locations.
	 refseq2gene.List <- split(row.names(refseq2gene),refseq2gene$TXNAME)
	 dup.ID <- as.numeric(unlist(refseq2gene.List[elementLengths(refseq2gene.List) > 1]))
	 dup.GENEID <- refseq2gene[dup.ID,"GENEID"]
	 dup.anno <- refseq2gene[!is.na(match(refseq2gene$GENEID,dup.GENEID)),]
	 dup.anno <- dup.anno[order(dup.anno$GENEID,dup.anno$TXID),]; row.names(dup.anno) <- 1:dim(dup.anno)[1]		
	 #####Add "_dup" to these multi-loci genes
	 count=1
	 for (i in 1:(dim(dup.anno)[1]-1)) {
	 	if (dup.anno$GENEID[i] == dup.anno$GENEID[i+1]) 
	 	{
	 		#temp <-	paste(dup.anno$GENEID[i],"_",count,sep="")
	 		if ((dup.anno$TXID[i+1]-dup.anno$TXID[i]) == 1) {	dup.anno$NEWID[i] <- paste(dup.anno$GENEID[i],"_",count,sep="")}
	 		if ((dup.anno$TXID[i+1]-dup.anno$TXID[i]) > 1) 
	 			 {dup.anno$NEWID[i] <- paste(dup.anno$GENEID[i],"_",count,sep=""); count=count+1; dup.anno$NEWID[i+1] <- paste(dup.anno$GENEID[i+1],"_",count,sep="")}
	 	}
	 	if (dup.anno$GENEID[i] != dup.anno$GENEID[i+1])		
	 	{ dup.anno$NEWID[i] <- paste(dup.anno$GENEID[i],"_",count,sep=""); count=1 }
	 }
}	
	 
###################################################
### code chunk number 4: Get Exon Info to do 
### collapsing, run seperatedly for ensGene and refGene
###################################################
#get all exon information
print ("Generating Synthetic Exon")
anno.refseq <- select(tx, keys = txname, columns=c("TXCHROM","TXSTRAND","EXONSTART","EXONEND","TXNAME","GENEID","TXID"),keytype="TXNAME")
############For ensGene############
if (args[2] == "ensGene") {
	 #generate rangeList object to do reduce (collapse exon)	
	 rngList.exon <- split(IRanges(start=anno.refseq$EXONSTART,end=anno.refseq$EXONEND),factor(anno.refseq$GENEID, levels=unique(anno.refseq$GENEID)))
	 print ("Exon Number Before Collapsing")
	 print(summary(elementLengths(rngList.exon)))
	 #Reduce
	 rngList.exon<- reduce(rngList.exon)
	 print ("Exon Number After Collapsing")
	 print(summary(elementLengths(rngList.exon)))
	 ####preparing output bed file
	 syntheticGeneModel <- anno.refseq[rep(
	   match(names(rngList.exon),anno.refseq$GENEID), 
	   elementLengths(rngList.exon)),] #Parent
	 #update the coordinates
	 syntheticGeneModel$EXONSTART<- unlist(start(rngList.exon))
	 syntheticGeneModel$EXONEND<- unlist(end(rngList.exon))
	 #options(scipen=10)
	 #print ("Writing Collapsed Exon Bed")
	 #syntheticGeneModel.output <- paste(output,"syntheticGeneModel.bed",sep="/")
	 #write.table(syntheticGeneModel[,c(6,2,3,4,5,7)],syntheticGeneModel.output,sep="\t",row.names=F,col.names=F,quote=F)
	 ######Add splicing sites for each exon (for the first exon,3'end, for the last exon, only 5'end)
	 #get exon order and number
	 ExonID <- lapply(elementLengths(rngList.exon),":",1)
	 sel <- syntheticGeneModel$TXSTRAND[cumsum(elementLengths(rngList.exon))] == "+"
	 ExonID[sel] <- sapply(ExonID[sel],rev)
	 syntheticGeneModel$ExonID <- unlist(ExonID)
	 syntheticGeneModel$ExonNum <- rep(elementLengths(rngList.exon), times=elementLengths(rngList.exon))
	 #Extend 2bp for genes > 1 exons
	 sel.FirstPlus <- syntheticGeneModel$ExonID ==1 & syntheticGeneModel$TXSTRAND == "+" & syntheticGeneModel$ExonNum > 1 #First exon in plus, end + 2bp
	 syntheticGeneModel[sel.FirstPlus, ]$EXONEND = syntheticGeneModel[sel.FirstPlus, ]$EXONEND + 2
	 sel.FirstMinus <- syntheticGeneModel$ExonID ==1 & syntheticGeneModel$TXSTRAND == "-" & syntheticGeneModel$ExonNum > 1 #First exon in minus, start - 2bp
	 syntheticGeneModel[sel.FirstMinus, ]$EXONSTART= syntheticGeneModel[sel.FirstMinus, ]$EXONSTART - 2
	 sel.LastPlus <- syntheticGeneModel$ExonID == syntheticGeneModel$ExonNum & syntheticGeneModel$TXSTRAND == "+" & syntheticGeneModel$ExonNum > 1 #Last exon in plus, start - 2bp
	 syntheticGeneModel[sel.LastPlus, ]$EXONSTART = syntheticGeneModel[sel.LastPlus, ]$EXONSTART - 2
	 sel.LastMinus <- syntheticGeneModel$ExonID == syntheticGeneModel$ExonNum & syntheticGeneModel$TXSTRAND == "-" & syntheticGeneModel$ExonNum > 1 #Last exon in minus, end + 2bp
	 syntheticGeneModel[sel.LastMinus, ]$EXONEND = syntheticGeneModel[sel.LastMinus, ]$EXONEND + 2
	 rest <- syntheticGeneModel$ExonNum > 1 & !sel.FirstPlus & !sel.FirstMinus & !sel.LastPlus & !sel.LastMinus #For rest exons, start - 2bp, end + 2bp
	 syntheticGeneModel[rest, ]$EXONSTART <- syntheticGeneModel[rest, ]$EXONSTART -2; syntheticGeneModel[rest, ]$EXONEND <- syntheticGeneModel[rest, ]$EXONEND + 2
	 options(scipen=10)
	 print ("Writing Collapsed Exon Bed with splicing sites")
	 syntheticGeneModel.output.ss <- paste(output,"syntheticGeneModel.bed",sep="/")
	 write.table(syntheticGeneModel[,c(6,2,3,4,5,7)],syntheticGeneModel.output.ss,sep="\t",row.names=F,col.names=F,quote=F)
}
############For refGene############
if (args[2] == "refGene") {
	 ###First replace geneID with geneID+"_num" (NEWID)
	 anno.refseq$NEWID <- anno.refseq$GENEID
	 temp <- match(anno.refseq$TXID,dup.anno[,"TXID"])
	 anno.refseq[!is.na(temp),"NEWID"] <- dup.anno[temp[!is.na(temp)],"NEWID"]
	 #generate rangeList object to do reduce (collapse exon)	
	 rngList.exon <- split(IRanges(start=anno.refseq$EXONSTART,end=anno.refseq$EXONEND),factor(anno.refseq$NEWID, levels=unique(anno.refseq$NEWID)))
	 print ("Exon Number Before Collapsing")
	 print(summary(elementLengths(rngList.exon)))
	 #####Reduce
	 rngList.exon<- reduce(rngList.exon)
	 print ("Exon Number After Collapsing")
	 print (summary(elementLengths(rngList.exon)))
	 ######Preparing output bed file
	 syntheticGeneModel <- anno.refseq[rep(
	 match(names(rngList.exon),anno.refseq$NEWID), 
	 elementLengths(rngList.exon)),] #Parent
	 #update the coordinates
	 syntheticGeneModel$EXONSTART <- unlist(start(rngList.exon))
	 syntheticGeneModel$EXONEND <- unlist(end(rngList.exon))
	 #options(scipen=10)
	 #print ("Writing Collapsed Exon Bed")
	 #syntheticGeneModel.output <- paste(output,"syntheticGeneModel.bed",sep="/")
	 #write.table(syntheticGeneModel[,c(6,2,3,8,5,7)],syntheticGeneModel.output,sep="\t",row.names=F,col.names=F,quote=F)
	 ######Add splicing sites for each exon (for the first exon,3'end, for the last exon, only 5'end)
	 #get exon order and number
	 ExonID <- lapply(elementLengths(rngList.exon),":",1)
	 sel <- syntheticGeneModel$TXSTRAND[cumsum(elementLengths(rngList.exon))] == "+"
	 ExonID[sel] <- sapply(ExonID[sel],rev)
	 syntheticGeneModel$ExonID <- unlist(ExonID)
	 syntheticGeneModel$ExonNum <- rep(elementLengths(rngList.exon), times=elementLengths(rngList.exon))
	 #Extend 2bp for genes > 1 exons
	 sel.FirstPlus <- syntheticGeneModel$ExonID ==1 & syntheticGeneModel$TXSTRAND == "+" & syntheticGeneModel$ExonNum > 1 #First exon in plus, end + 2bp
	 syntheticGeneModel[sel.FirstPlus, ]$EXONEND = syntheticGeneModel[sel.FirstPlus, ]$EXONEND + 2
	 sel.FirstMinus <- syntheticGeneModel$ExonID ==1 & syntheticGeneModel$TXSTRAND == "-" & syntheticGeneModel$ExonNum > 1 #First exon in minus, start - 2bp
	 syntheticGeneModel[sel.FirstMinus, ]$EXONSTART= syntheticGeneModel[sel.FirstMinus, ]$EXONSTART - 2
	 sel.LastPlus <- syntheticGeneModel$ExonID == syntheticGeneModel$ExonNum & syntheticGeneModel$TXSTRAND == "+" & syntheticGeneModel$ExonNum > 1 #Last exon in plus, start - 2bp
	 syntheticGeneModel[sel.LastPlus, ]$EXONSTART = syntheticGeneModel[sel.LastPlus, ]$EXONSTART - 2
	 sel.LastMinus <- syntheticGeneModel$ExonID == syntheticGeneModel$ExonNum & syntheticGeneModel$TXSTRAND == "-" & syntheticGeneModel$ExonNum > 1 #Last exon in minus, end + 2bp
	 syntheticGeneModel[sel.LastMinus, ]$EXONEND = syntheticGeneModel[sel.LastMinus, ]$EXONEND + 2
	 rest <- syntheticGeneModel$ExonNum > 1 & !sel.FirstPlus & !sel.FirstMinus & !sel.LastPlus & !sel.LastMinus #For rest exons, start - 2bp, end + 2bp
	 syntheticGeneModel[rest, ]$EXONSTART <- syntheticGeneModel[rest, ]$EXONSTART -2; syntheticGeneModel[rest, ]$EXONEND <- syntheticGeneModel[rest, ]$EXONEND + 2
	 options(scipen=10)
	 print ("Writing Collapsed Exon Bed with splicing sites")
	 syntheticGeneModel.output.ss <- paste(output,"syntheticGeneModel.bed",sep="/")
	 write.table(syntheticGeneModel[,c(6,2,3,8,5,7)],syntheticGeneModel.output.ss,sep="\t",row.names=F,col.names=F,quote=F)
}


###################################################
### code chunk number 5:get geneInterval bed could 
### also use shell to generate
### but here by R in order to get the reduced exon
###################################################
#Eg in Shell: cat /cclab_nas/Reference/hg38/gff/hg38.RefSeq.gff |awk '{if ($3=="transcript"||$3=="mRNA") print $0}' >hg38.Refseq_GeneInterval.gff
print ("Create Gene Interval Bed...")
############For ensGene############
if (args[2] == "ensGene") {
	 GeneIntervalModel <- syntheticGeneModel[,c(6,2,3,4,5,7)]; 
	 GeneIntervalModel$EXONSTART <- 0.0000; GeneIntervalModel$EXONEND <- 0.0000;GeneIntervalModel$TXID <- 0.0000
	 GeneIntervalModel <- unique(GeneIntervalModel)
	 GeneIntervalModel$EXONSTART <- unlist(bplapply(rngList.exon,function(x){min(unlist(start(x)))}, BPPARAM= MulticoreParam(workers = NCORE)))
	 GeneIntervalModel$EXONEND <- unlist(bplapply(rngList.exon,function(x){max(unlist(end(x)))}, BPPARAM= MulticoreParam(workers = NCORE)))
	 print ("Writing Gene Interval Bed")
	 GeneIntervalModel.output <- paste(output,"GeneIntervalModel.bed",sep="/")
	 options(scipen=10)
	 GeneIntervalModel <- GeneIntervalModel[order(GeneIntervalModel[,1], GeneIntervalModel[,2]),] #From version 2.22 of bedtools, input must be sorted
	 write.table(GeneIntervalModel,GeneIntervalModel.output,sep="\t",row.names=F,col.names=F,quote=F)
}
############For refGene############
if (args[2] == "refGene") {
	 GeneIntervalModel <- syntheticGeneModel[,c(6,2,3,8,5,7)]; 
	 GeneIntervalModel$EXONSTART <- 0.0000; GeneIntervalModel$EXONEND <- 0.0000;GeneIntervalModel$TXID <- 0.0000
	 GeneIntervalModel <- unique(GeneIntervalModel)
	 GeneIntervalModel$EXONSTART <- unlist(bplapply(rngList.exon,function(x){min(unlist(start(x)))}, BPPARAM= MulticoreParam(workers = NCORE)))
	 GeneIntervalModel$EXONEND <- unlist(bplapply(rngList.exon,function(x){max(unlist(end(x)))}, BPPARAM= MulticoreParam(workers = NCORE)))
	 print ("Writing Gene Interval Bed")
	 GeneIntervalModel.output <- paste(output,"GeneIntervalModel.bed",sep="/")
 	 GeneIntervalModel <- GeneIntervalModel[order(GeneIntervalModel[,1], GeneIntervalModel[,2]),] #From version 2.22 of bedtools, input must be sorted 
	 options(scipen=10)
	 write.table(GeneIntervalModel,GeneIntervalModel.output,sep="\t",row.names=F,col.names=F,quote=F)     
}	 
	 
###################################################
###### Code chunk number 6:get reduced intron
###################################################
#Remove transcripts with only 1 exon
print ("Create Intron Bed, this may take a while...")
############For ensGene############
if (args[2] == "ensGene") {
	 rngList.GeneInterval <- split(IRanges(start=GeneIntervalModel$EXONSTART,end=GeneIntervalModel$EXONEND),factor(GeneIntervalModel$GENEID,levels=unique(GeneIntervalModel$GENEID)),drop=T)
	 #this step takes a lot of time, should be improved															
	 rngList.intron <- bpmapply(setdiff,rngList.GeneInterval,rngList.exon, BPPARAM= MulticoreParam(workers = NCORE)) # get intronic coordinates
	 ######get annotation (strand,chr, ...)
	 IntronGeneModel <- GeneIntervalModel[rep(match(names(rngList.intron),GeneIntervalModel$GENEID),elementLengths(rngList.intron)),] #Parent                   															                                
	 ######update the coordinates
	 IntronGeneModel[,2]<- as.numeric(unlist(lapply(rngList.intron,start)))
	 IntronGeneModel[,3]<- as.numeric(unlist(lapply(rngList.intron,end)))
	 print ("Writing Intron Bed")
	 IntronGeneModel.output <- paste(output,"IntronGeneModel.bed",sep="/")
	 options(scipen=10)  
	 write.table(IntronGeneModel,IntronGeneModel.output,sep="\t",row.names=F,col.names=F,quote=F)
}
############For refGene############	 		
if (args[2] == "refGene") {        	 
	 rngList.GeneInterval <- split(IRanges(start=GeneIntervalModel$EXONSTART,end=GeneIntervalModel$EXONEND),factor(GeneIntervalModel$NEWID,levels=unique(GeneIntervalModel$NEWID)),drop=T)
	 #this step takes a lot of time, should be improved															
	 rngList.intron <- bpmapply(setdiff,rngList.GeneInterval,rngList.exon, BPPARAM= MulticoreParam(workers = NCORE)) # get intronic coordinates
	 ######get annotation (strand,chr, ...)
	 IntronGeneModel <- GeneIntervalModel[rep(match(names(rngList.intron),GeneIntervalModel$NEWID),elementLengths(rngList.intron)),] #Parent                   															                                
	 ######update the coordinates
	 IntronGeneModel[,2]<- as.numeric(unlist(lapply(rngList.intron,start)))
	 IntronGeneModel[,3]<- as.numeric(unlist(lapply(rngList.intron,end)))
	 print ("Writing Intron Bed")
	 IntronGeneModel.output <- paste(output,"IntronGeneModel.bed",sep="/")
	 options(scipen=10)  
	 write.table(IntronGeneModel,IntronGeneModel.output,sep="\t",row.names=F,col.names=F,quote=F)	
}
