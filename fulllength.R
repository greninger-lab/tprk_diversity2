# Creates full length ORFs that can be fed into tree.

list.of.packages <- c("treeio","ggtree","stringr","Biostrings","phylobase","pegas","tidyverse","lubridate","ape",
                      "plyr","phangorn","RColorBrewer","dplyr","optparse","data.table","tidyr", "BiocGenerics")
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))

option_list <- list(make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-r", "--relativefreq"), type="character", default=NULL, help="relative frequency cutoff", metavar="character"),
                    make_option(c("-c", "--count"), type="character", default=NULL, help="count cutoff", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/new_tprk/demuxed/20200820_pacbio/chimeric_stuff/redo-with-actual-truncated/')

#path <- opt$directory
path <- "./"

if (opt$relativefreq) {
  rf_cutoff = as.numeric(opt$relativefreq)
} else {
  rf_cutoff = 0.2
}

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are.

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files_redo"
#####
# Grabs sample names and PacBio files from the metadata file
metadata <- read.table(opt$metadata, sep=',', header=TRUE)
metadata <- read.table("metadata.csv", sep=',', header=TRUE)

PacBio_fns <- RAD_files_fix
PacBio_fns <- c(as.character(metadata$PacBio))
sample_names <- c(as.character(metadata$SampleName))

## Function to convert fasta to dataframe
BString2df=function(BString){
  names = names(BString)
  sequences = paste(BString)
  df <- data.frame(names,sequences)
  return(df)
}

## Function to convert fasta to dataframe.  Assumes the two columns are [names,sequences].
df2BString=function(df){
  BString <- BStringSet(df$sequences)
  names(BString) = paste0(df$names)
  return(BString)
}


## Function to remove translated ORFs with stop codons
removeTruncatedORF=function(AAFile){
  AAFile[,2]$sequences<-as.character(AAFile[,2]$sequences)
  
  for (num in c(1:nrow(AAFile))){
    if (!(regexpr("\\*", AAFile[num,2]) == nchar(AAFile[num,2]))){
      AAFile[num,] <- NA}
    print("deleted")
  }
  return(AAFile)
}

## Function to remove translated ORFs with stop codons
removeFrameShift=function(AAFile){
  conserved <- c("HGFKTTTDFKIVFPIVAKKD",
                 "RTREDGVQEYIKVELTGNST",
                 "VGAKVSMKLWGLCALAATDV",
                 "ADALLTLGYRWFSAGGYFAS",
                 "LETKGSDPDTSFLEGLDLGV",
                 "YFPVYGKVWGSYRHDMGEYG",
                 "WEQGKLQENSNVVIEKNVTE")
  AAFile[,2]$sequences<-as.character(AAFile[,2]$sequences)
  for (num in c(1:nrow(AAFile))){
    frameShift <- FALSE
    for (i in c(1:length(conserved))){
      if (agrepl(conserved[i],AAFile[num,2],max.distance=3,ignore.case=TRUE) == 0){
        frameShift <- TRUE
      }
    }
    if (frameShift){
      AAFile[num,] <- NA
    }
  }
  return(AAFile)
}

fasta_files <- list()
df_list <- list()
aa_list <- list()
df_aa_list <- list()
df_aa_filtered_list <- list()
allAA <- list()
allAAfilt <- list()

for (i in 1:length(PacBio_fns)) {
  #fastafile_name <- paste((substr(PacBio_fns[i],1,nchar(PacBio_fns[i])-6)),".noprimers.filtered.RAD.nolines.fix.fasta",sep="")
  fastafile_name <- PacBio_fns[i]
  fastafile <- reverseComplement(readDNAStringSet(fastafile_name))
  names(fastafile) <- paste(sample_names[i],"_",names(fastafile),sep="")
  fasta_files <- c(fasta_files,fastafile)
  df_list <- c(df_list,BString2df(fastafile))
  amino_acids <- Biostrings::translate(fastafile, getGeneticCode("11"), no.init.codon=FALSE,
                                       if.fuzzy.codon="error")
  aa_list <- c(aa_list, amino_acids)
  df_aa <- BString2df(amino_acids)
  
  df_aa <- mutate(df_aa,sample=sapply(strsplit(as.character(names),"_"),"[",1),
                  count=as.numeric(sapply(strsplit(as.character(names),"_"),"[",3)))
  
  ##Combine same TprKs
  df_aa <- data.table(df_aa)
  df_aa <- df_aa[order(-df_aa$count),]
  df_aa <- df_aa[,list(names,sample,count=sum(count)),by='sequences']
  df_aa <- df_aa[!duplicated(df_aa$sequences),]
  df_aa$names <- paste(df_aa$sample,sapply(strsplit(as.character(df_aa$names), split="_"), "[", 2),df_aa$count,sep="_")
  df_aa <- df_aa[,c(2,1,3,4)]
  df_aa <- mutate(df_aa,percentage=round(count / sum(count)*100,3))
  df_aa_filt <- filter(df_aa,percentage>= rf_cutoff)
  df_aa_list[[i]] <- df_aa
  df_aa_filtered_list[[i]] <- df_aa_filt
  if (i == 1) {
    allAA <- df_aa
    allAAfilt <- df_aa_filt
  } else {
    allAA <- base::rbind(allAA, df_aa)
    allAAfilt <- base::rbind(allAAfilt, df_aa_filt)
  }
}

## Remove TprKs that are not complete ORFs
allAA_fullORFs_df <- drop_na(removeTruncatedORF(allAA))
allAAfilt_fullORFs_df <- drop_na(removeTruncatedORF(allAAfilt))

## Remove TprKs that have two frameshifts that put them 
allAA_fullORFs_df <- drop_na(removeFrameShift(allAA_fullORFs_df))
allAAfilt_fullORFs_df <- drop_na(removeFrameShift(allAAfilt_fullORFs_df))

allAAfilt_fullORFs_df[which(duplicated(allAAfilt_fullORFs_df$sequences) == TRUE),]

write.table(allAAfilt_fullORFs_df,"Table_allAAfilt_fullORFs.tsv",sep='\t',row.names=FALSE)
allAA_fullORFs_BString <- df2BString(allAA_fullORFs_df)
allAAfilt_fullORFs_BString <- df2BString(allAAfilt_fullORFs_df)

AAoutfile <- "Isolates_aa_fullORFs.fasta"
AAoutfilefilt <- "Isolates_aa_filt_fullORFs.fasta"

writeXStringSet(allAA_fullORFs_BString, AAoutfile, append=FALSE,format="fasta")
writeXStringSet(allAAfilt_fullORFs_BString, paste(AAoutfilefilt,sep=""), append=FALSE,format="fasta")

## Need to include our Illumina check right here!  ###

AAaln_outfile <- paste0(substr(AAoutfile,1,nchar(AAoutfile)-5),"aln.fasta")
AAaln_outfile_filt <- paste0(substr(AAoutfilefilt,1,nchar(AAoutfilefilt)-5),"aln.fasta")

#mafft_command <- paste0("mafft --auto ",AAoutfile," > ",AAaln_outfile)
#system(mafft_command)

mafft_command <- paste0("/usr/local/miniconda/bin/mafft --auto ",AAoutfilefilt," > ",AAaln_outfile_filt)
system(mafft_command)