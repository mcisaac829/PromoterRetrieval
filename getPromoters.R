#source("https://bioconductor.org/biocLite.R")
#biocLite()
library(Biostrings)
library(BSgenome)
require(IRanges)

setwd('/Applications/Bioinformatics/REDUCE_Suite/examples/dgh_2018_0221/ComplexGenes')

MyGenome <- "SacCer3.fasta"
GeneExpressionData <- "data_forFig2.txt"
toMatch <- c("CUT","XUT","MUT","SUT","AS")
BedFile <- "complex_annotation_step3.bed"

#output promoters of genes containing strings shown in "toMatch" to promoters.txt
#merged promoters that overlap and output to regions.txt

#load genome 
genome <- read.delim(file=MyGenome)
PromoterLength = 300

#load gene expression data, to get names of non-canonical transcripts --> complexGenes
data <- read.delim(file = GeneExpressionData,row.names=1,header=TRUE)
complexGenes <- row.names(data[grep(paste(toMatch,collapse="|"), row.names(data)),])
numComplexGenes <- length(complexGenes)

# get chromosome #, orientation, start, and finish 
bedFile <- as.data.frame(read.table(BedFile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
complexGenesBed <- subset(bedFile, V4 %in% complexGenes)[,c(1,2,3,4,6)]
colnames(complexGenesBed) = c("Chromosome","Beginning","End","Gene","Strand")


# use positions to get promoter sequences
YeastGenome <- readDNAStringSet("SacCer3.fasta")

#subseq(YeastGenome$mitochondrion,start=1,end=5)

df <- data.frame(PromoterStart = integer(),
                 PromoterEnd = integer(),
                 Orientation = character(),
                 Chromosome = character(),
                 Gene = character(),
                 stringsAsFactors=FALSE)

#####get promoters and write to promoters.txt
sink("promoters.txt")
for(entry in 1:numComplexGenes){
#get the chromosome sequence
getChromosome = YeastGenome[names(YeastGenome) == complexGenesBed[entry,"Chromosome"]]
#get promoter beginning and end 
if(complexGenesBed[entry,"Strand"]=="+"){
  ORFstart = complexGenesBed[entry,"Beginning"]
  PromStart = ORFstart - PromoterLength
  PromEnd = ORFstart - 1
  PromoterSequence <- paste0(toString(subseq(getChromosome,start=PromStart,end=PromEnd)),"\n")
} else {
  ORFStart = complexGenesBed[entry,"End"]
  PromStart = ORFStart + 1
  PromEnd = ORFStart + PromoterLength
  PromoterSequence <- paste0(toString(reverseComplement(subseq(getChromosome,start=PromStart,end=min(width(getChromosome),PromEnd)))),"\n")
}

#extract promoter 
PromoterName <- paste0(">",complexGenesBed[entry,"Gene"]," ",
                       complexGenesBed[entry,"Chromosome"]," ",
                       complexGenesBed[entry,"Strand"]," ",
                       ", promoter sequence, start = ",PromStart,", end = ",PromEnd,"\n")

df[entry,1] = PromStart
df[entry,2] = PromEnd
df[entry,3] = complexGenesBed[entry,"Strand"]
df[entry,4] = complexGenesBed[entry,"Chromosome"]
df[entry,5] = complexGenesBed[entry,"Gene"]

cat(PromoterName)
cat(PromoterSequence)
}
sink()






rm(list = setdiff(ls(),"df")) #only keep df in memory
MyGenome <- "SacCer3.fasta"

#
mylist.names = unique(df$Chromosome)
mylist <- vector("list", length(mylist.names)) # to contain all of the intervals
names(mylist) <- mylist.names

for(k in 1:length(names(mylist))){
  Chr <- names(mylist)[k]
  mylist[[Chr]] <- df[df$Chromosome==Chr,]
}


mylist2 <- lapply(mylist,function(x){split(x,f=x$Orientation)}) #to contain merged intervals
for(k in 1:length(mylist2)){
  J <-mylist2[[k]]
  J$`-` <- as.data.frame(IRanges::reduce(IRanges(J$`-`$PromoterStart, J$`-`$PromoterEnd)))
  J$`+` <- as.data.frame(IRanges::reduce(IRanges(J$`+`$PromoterStart, J$`+`$PromoterEnd)))
  mylist2[[k]] <- J
}

#get data from genome and print to FASTA
YeastGenome <- readDNAStringSet(MyGenome)

sink("regions.txt")
for(k in 1:length(mylist2)){
  
  Chr <- names(mylist2[k])
  
  getChromosome <- YeastGenome[Chr]
  
  df <- mylist2[[Chr]] #slice  region list
  
  for(j in 1:dim(df$`+`)[1]){
    RegionName <- paste0(">",Chr,
                         ", (+ strand) start = ",df$`+`[j,1],", end = ",df$`+`[j,2],"\n")
    RegionSequence <- paste0(toString(subseq(getChromosome,start=df$`+`[j,1],end=df$`+`[j,2])),"\n")
    cat(RegionName)
    cat(RegionSequence)
  }
  
  for(i in 1:dim(df$`-`)[1]){
    RegionName <- paste0(">",Chr,
                         ", (- strand) start = ",df$`-`[i,1],", end = ",df$`-`[i,2],"\n")
    RegionSequence <- paste0(toString(reverseComplement(subseq(getChromosome,start=df$`-`[i,1],end=min(width(getChromosome),df$`-`[i,2])))),"\n")
    cat(RegionName)
    cat(RegionSequence)
  }
  
  
}
sink()
sink()




#J$`-`$group2 <- subjectHits(IRanges::findOverlaps(ir,IRanges::reduce(ir)))

# #####create distance information 
# 
# #create list that will have distance matrices
# mylist.names = unique(df$Chromosome)
# mylist <- vector("list", length(mylist.names)) # to contain all of the distance matrices
# names(mylist) <- mylist.names
# 
# for(k in 1:length(names(mylist))){
#   Chr <- names(mylist)[k]
#   d <- df[df$Chromosome==Chr,]
#   m <- matrix(nrow = dim(d)[1],ncol=dim(d)[1])
#   for(i in 1:dim(d)[1]){
#     for(j in i:dim(d)[1]){
#     
#       if(d[i,2]==d[j,2]){
#        m[i,j] <- d[i,1] - d[j,1]} else {
#           m[i,j] <- 1000000
#        }
#     
#     }
#     
#   }
#   mylist[[Chr]] <- m
#   print(k)
#   print(d)
# }
# 
# 

df_total <-data.frame()
for(k in 1:length(names(mylist))){
  df = sapply(mylist2[[k]],NROW) 
  df_total <- rbind(df_total,df)
  }
