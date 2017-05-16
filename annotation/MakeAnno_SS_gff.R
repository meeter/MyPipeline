####Using gff from transcriptomic reconstruction
####R CMD BATCH combineGFF.R input_gff output_folder Num_Threads


###################################################
### code chunk number 1: load the library
###################################################
GI <- require(genomeIntervals,quietly = T)
IR <- require(IRanges, quietly = T)
BP <- require(BiocParallel,warn.conflicts = F, quietly = T)
if (GI & IR & BP) {print("SUCCESS IN LOADING PACKAGES")
	      }else{stop("ERROR IN LOADING PACKAGES")}
     
#load("GRCh37.75/SyntheticExon.gff.RData")
#save.image("GRCh37.75/SyntheticExon.gff.RData")

###################################################
### code chunk number 2: read gff3
###################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3){stop ("Too Many Arguments, please refer to readme")}    
if (length(args) < 2){stop ("Too Few Arguments, please refer to readme")}
if (args[3] == '') {args[3] <- 1}
print ("Loading Gff file")
gff <- readGff3(args[1])
#gff <- readGff3("merged_j.gff")
#output file name
SyntheticExon.output <- paste(args[2], 'syntheticGeneModel.bed', sep="/")
GeneInterval.output <-  paste(args[2], 'GeneIntervalModel.bed', sep="/")
IntronInterval.output <- paste(args[2], 'IntronGeneModel.bed',sep="/")

#get transcript <-> Gene Mapping
transcriptGeneMapping <- data.frame(
  getGffAttribute(gff[gff$type=="transcript"|gff$type=="mRNA" ,],"ID"), #
  getGffAttribute(gff[gff$type=="transcript"|gff$type=="mRNA" ,],"gene_id")) #;Parent
#head(transcriptGeneMapping)


###################################################
### code chunk number 3: ranges List
###################################################
print ("Generating Synthetic Exon")
sel <- gff$type=="exon"
rngList <- split(IRanges(start=gff[sel,1],end=gff[sel,2]),
                transcriptGeneMapping[match(
                  sapply(strsplit(getGffAttribute(gff[sel,],"Parent"),","),"[",1),
                  transcriptGeneMapping$ID),"gene_id"]) #Parent                 
                  
                
#rngList
print ("Exon Number and Length Before Collapsing")
summary(elementLengths(rngList))
summary(unlist(width(rngList)))

###################################################
### code chunk number 4: reduce
###################################################
rngList<- IRanges::reduce(rngList)
#rngList
print ("Exon Number and Length After Collapsing")
summary(elementLengths(rngList))
summary(unlist(width(rngList)))


###################################################
### code chunk number 5: create the gff for output
###################################################
## get the gff information; here we simply duplicate the
## first exon of every gene by the number of synthetic 
## exons per gene. The content will be updated afterwards.
exons <- gff[sel,]
syntheticGene.gff3<- exons[rep(
  match(names(rngList),
  transcriptGeneMapping[
    match(sapply(strsplit(getGffAttribute(exons,"Parent"),","),"[",1),
    transcriptGeneMapping$ID),"gene_id"]),  #Parent
  elementLengths(rngList)),]

## update the coordinates
syntheticGene.gff3[,1]<- unlist(start(rngList))
syntheticGene.gff3[,2]<- unlist(end(rngList))

###### change the source
#levels(syntheticGene.gff3$source)<- "inhouse"

###### get the exon number for the minus strand
exonNumber <- lapply(elementLengths(rngList),":",1)

###### reverse them on the plus strand
sel <- strand(syntheticGene.gff3)[cumsum(elementLengths(rngList))] == "+"
exonNumber[sel] <- sapply(exonNumber[sel],rev)

###### update the attributes
syntheticGene.gff3$gffAttributes<- paste("ID=",
  rep(names(rngList),elementLengths(rngList)),
  ":",unlist(exonNumber),";Parent=",
  rep(names(rngList),elementLengths(rngList)),sep="") #".0",

######Add splicing sites for each exon (for the first exon,3'end, for the last exon, only 5'end)
	 #get exon order and number
	 #ExonID <- lapply(elementLengths(rngList),":",1)
	 #sel <- syntheticGene.gff3$TXSTRAND[cumsum(elementLengths(rngList))] == "+"
	 #ExonID[sel] <- sapply(ExonID[sel],rev)
	 syntheticGene.gff3$score <- unlist(exonNumber) #score <- ExonID
	 syntheticGene.gff3$phase <- rep(elementLengths(rngList), times=elementLengths(rngList)) #phase <- ExonNumber
	 #Extend 2bp for genes > 1 exons
	 sel.FirstPlus <- syntheticGene.gff3$score ==1 & syntheticGene.gff3$strand == "+" & syntheticGene.gff3$phase > 1 #First exon in plus, end + 2bp
	 syntheticGene.gff3[sel.FirstPlus, 2] = syntheticGene.gff3[sel.FirstPlus, 2] + 2
	 sel.FirstMinus <- syntheticGene.gff3$score ==1 & syntheticGene.gff3$strand == "-" & syntheticGene.gff3$phase > 1 #First exon in minus, start - 2bp
	 syntheticGene.gff3[sel.FirstMinus, 1]= syntheticGene.gff3[sel.FirstMinus, 1] - 2
	 sel.LastPlus <- syntheticGene.gff3$score == syntheticGene.gff3$phase & syntheticGene.gff3$strand == "+" & syntheticGene.gff3$phase > 1 #Last exon in plus, start - 2bp
	 syntheticGene.gff3[sel.LastPlus, 1] = syntheticGene.gff3[sel.LastPlus, 1] - 2
	 sel.LastMinus <- syntheticGene.gff3$score == syntheticGene.gff3$phase & syntheticGene.gff3$strand == "-" & syntheticGene.gff3$phase > 1 #Last exon in minus, end + 2bp
	 syntheticGene.gff3[sel.LastMinus, 2] = syntheticGene.gff3[sel.LastMinus, 2] + 2
	 rest <- syntheticGene.gff3$phase > 1 & !sel.FirstPlus & !sel.FirstMinus & !sel.LastPlus & !sel.LastMinus #For rest exons, start - 2bp, end + 2bp
	 syntheticGene.gff3[rest, 1] <- syntheticGene.gff3[rest, 1] -2; syntheticGene.gff3[rest, 2] <- syntheticGene.gff3[rest, 2] + 2


## write the file
print ("Writing Collapsed Exon Gff")
#writeGff3(syntheticGene.gff3,file=SyntheticExon.output)
#syntheticGene.gff3$seq_name <- paste("chr",syntheticGene.gff3$seq_name,sep="")
syntheticGeneModel <- cbind(as.character(syntheticGene.gff3$seq_name), syntheticGene.gff3[,1]-1,syntheticGene.gff3[,2]+1, 
                            getGffAttribute(syntheticGene.gff3, "Parent"), syntheticGene.gff3$phase, as.character(syntheticGene.gff3$strand))
syntheticGeneModel[,4] <- sapply(strsplit(syntheticGeneModel[,4], "\\."), "[",1)

colnames(syntheticGeneModel) <- rep("",6)
syntheticGeneModel <- data.frame(syntheticGeneModel)
syntheticGeneModel[,2] <- gsub(" ","", format(as.numeric(as.character(syntheticGeneModel[,2])), scientific = FALSE))
syntheticGeneModel[,3] <- gsub(" ","", format(as.numeric(as.character(syntheticGeneModel[,3])), scientific = FALSE))
write.table(syntheticGeneModel, SyntheticExon.output, sep="\t",row.names=F,col.names=F,quote=F)


###################################################
### code chunk number 6: create gff for GeneInterval
###################################################
print ("Create Gene Interval Gff, this may take a while...")
GeneInterval <- bplapply(rngList, function(x) {
return(IRanges(start = min(start(x)), end = max(end(x))))
}, BPPARAM= MulticoreParam(workers = as.numeric(args[3]))
)

                                                                                     
GeneInterval.gff3 <- syntheticGene.gff3[match(names(rngList),getGffAttribute(syntheticGene.gff3,"Parent")),]  
GeneInterval.gff3[,1]<- as.numeric(unlist(lapply(GeneInterval,start)))                                             
GeneInterval.gff3[,2]<- as.numeric(unlist(lapply(GeneInterval,end)))                                               
GeneInterval.gff3$type = 'GeneInterval'
print ("Writing Gene Interval Gff")                                                                 
#writeGff3(GeneInterval.gff3,file=GeneInterval.output)      
GeneIntervalModel <- cbind(as.character(GeneInterval.gff3$seq_name), GeneInterval.gff3[,1]-1, GeneInterval.gff3[,2]+1, 
                            getGffAttribute(GeneInterval.gff3, "Parent"), GeneInterval.gff3$phase, as.character(GeneInterval.gff3$strand))
GeneIntervalModel[,4] <- sapply(strsplit(GeneIntervalModel[,4], "\\."), "[",1)
colnames(GeneIntervalModel) <- rep("",6)
GeneIntervalModel <- data.frame(GeneIntervalModel)
GeneIntervalModel[,2] <- gsub(" ", "", format(as.numeric(as.character(GeneIntervalModel[,2])), scientific = FALSE))
GeneIntervalModel[,3] <- gsub(" ", "", format(as.numeric(as.character(GeneIntervalModel[,3])), scientific = FALSE))
write.table(GeneIntervalModel, GeneInterval.output, sep="\t",row.names=F,col.names=F,quote=F)                            


###################################################
### code chunk number 7: create gff for Small Intron
###################################################

#Remove transcripts with only 1 exon
print ("Create Intron Gff, this may take a while...")   
rngList.exon <-  split(IRanges(start=syntheticGene.gff3[,1],end=syntheticGene.gff3[,2]),
                  #transcriptGeneMapping[match(
                  #sapply(strsplit(getGffAttribute(syntheticGene.gff3,"Parent"),","),"[",1),
                  getGffAttribute(syntheticGene.gff3,"Parent"))
                  #transcriptGeneMapping$ID),"gene_name"]) #Parent          
                                                                                        
#IntronInterval <- mapply(setdiff,GeneInterval,rngList.exon) #39179     
IntronInterval <- bpmapply(setdiff,GeneInterval,rngList.exon, BPPARAM= MulticoreParam(workers = as.numeric(args[3]))) # get intronic coordinates  
IfMultiExon <- sapply(IntronInterval,function(x) {ifelse(length(x) == 0,FALSE,TRUE)} )                                
IfMultiExon <- which(IfMultiExon == TRUE)                                                                    
IntronInterval <- IntronInterval[match(names(IfMultiExon),names(IntronInterval))] #25255                                                                     
#get annotation                                                                                              
IntronInterval.gff3 <- syntheticGene.gff3[rep(match(names(IntronInterval),getGffAttribute(syntheticGene.gff3,"Parent")),  
                         elementLengths(IntronInterval)),]                                                                                                                               
#update the coordinates                                                                                    
IntronInterval.gff3[,1]<- as.numeric(unlist(lapply(IntronInterval,start)))                                                
IntronInterval.gff3[,2]<- as.numeric(unlist(lapply(IntronInterval,end)))                                                  
#get the intron number for the minus strand                                                                
intronNumber<- lapply(elementLengths(IntronInterval),":",1)                                                           
#reverse them on the plus strand                                                                           
sel<- strand(IntronInterval.gff3)[cumsum(elementLengths(IntronInterval))] == "+"                                          
intronNumber[sel]<- sapply(intronNumber[sel],rev)                                                            
#update the attributes                                                                                     
IntronInterval.gff3$gffAttributes<- paste("ID=",                                                                 
  rep(names(IntronInterval),elementLengths(IntronInterval)),                                                                   
  ":",unlist(intronNumber),";Parent=",                                                                       
  rep(names(IntronInterval),elementLengths(IntronInterval)),sep="") #".0",                                                                                                               
#write the file                                                                                            
IntronInterval.gff3$type = 'intron'                                                                              
print ("Writing Intron Gff")
#writeGff3(IntronInterval.gff3,file=IntronInterval.output)                                          
IntronIntervalModel <- cbind(as.character(IntronInterval.gff3$seq_name), IntronInterval.gff3[,1]-1, IntronInterval.gff3[,2]+1, 
                            getGffAttribute(IntronInterval.gff3, "Parent"), IntronInterval.gff3$phase, as.character(IntronInterval.gff3$strand))
IntronIntervalModel[,4] <- sapply(strsplit(IntronIntervalModel[,4], "\\."), "[",1)
colnames(IntronIntervalModel) <- rep("",6)
IntronIntervalModel <- data.frame(IntronIntervalModel)
IntronIntervalModel[,2] <- gsub(" ", "", format(as.numeric(as.character(IntronIntervalModel[,2])), scientific = FALSE))
IntronIntervalModel[,3] <- gsub(" ", "", format(as.numeric(as.character(IntronIntervalModel[,3])), scientific = FALSE))
write.table(IntronIntervalModel, IntronInterval.output, sep="\t",row.names=F,col.names=F,quote=F)                            
 
