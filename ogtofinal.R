## Load packages
list.of.packages <- c("JuliaCall", "reticulate", "devtools", "optparse", "dada2", "ggplot2", "ShortRead",
                      "reshape2", "optparse")
lapply(list.of.packages,require,character.only = TRUE)


## Points to Julia install in docker "quay.io/greninger-lab/tprk"
julia <- julia_setup(JULIA_HOME = "/Applications/Julia-0.6.app/Contents/Resources/julia/bin/")

##Specifying Illumina vs. PacBio files, and what the sample name is.
option_list <- list(make_option(c("-s", "--script_path"), type="character", default=NULL, help="Directory where scripts are located.", 
                                metavar="character"),
                    make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"),
                    make_option(c("-p", "--pacbio"), type="character", default=FALSE, help="Specify if these files are only PacBio.", 
                                metavar="character", action="store_true"),
                    make_option(c("-i", "--illumina"), type="character", default=FALSE, help="Specify if these files are only Illumina.", 
                                metavar="character", action="store_true"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$directory
script.dir <- '/Users/uwvirongs/Documents/tprk-master/'
script.dir <- opt$script_path

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are. (This is the -directory option).
## script.dir refers to the folder where all the script files are located. This should point to where you saved the cloned GitHub.

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files_redo3"
script.dir <- "/Users/uwvirongs/Documents/tprK-master/"

## This script can also be run from the command line.
## Usage: rscript \path\to\og_files_to_all_reads.R -s [script_path] -d [directory]

#####
path <- "./"
getwd()
setwd("/Users/uwvirongs/Documents/Michelle/tprk_files/illumina_check/rad_filtee1")
metadata <- read.table("metadata.csv", sep=',', header=TRUE)
metadata <- read.table(opt$metadata, sep=',', header=TRUE)
PacBio_fns <- c(as.character(metadata$PacBio))
Illumina_fns <- c(as.character(metadata$Illumina))
sample_names <- c(as.character(metadata$SampleName))
syph_path <- paste(script.dir,"/syph_r.py",sep='')

## Identify primers to go from ATG to stop in tprK
tprKF <- "GGAAAGAAAAGAACCATACATCC"
tprKR <- "CGCAGTTCCGGATTCTGA"
rc <- dada2:::rc
noprimer_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.fastq",sep ='')
nop <- file.path(noprimer_filenames)

if(opt$illumina == FALSE) {
  ## Remove primers
  for (count in c(1:length(nop))) {
    if(file.exists(nop[count])) {
      print(paste(noprimer_filenames[count], " already exists. Skipping removing primers step...", sep=""))
    } else {
      print("Removing primers from PacBio...")
      prim <- removePrimers(PacBio_fns, nop, primer.fwd=tprKF, primer.rev=rc(tprKR), orient=TRUE, verbose=TRUE)
    }
  }
  
  print("Filtering PacBio reads...")
  ## Setting up file names to filter.
  filter_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.fastq",sep ='')
  filterEE1_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.EE1.fastq",sep ='')
  filt <- file.path(filter_filenames)
  filtEE1 <- file.path(filterEE1_filenames)

    ## Filter reads for tprK length and do not worry about expected errors.
  for (count in c(1:length(filt))) {
    if (file.exists(filt[count])) {
      print(paste(filter_filenames[count]," already exists. Skipping filtering step..."), sep="")
    } else {
      print(paste("Filtering ",nop[count],"...",sep=""))
      track <- fastqFilter(nop[count], filt[count], minLen=1400,maxLen=1800,
                           maxN=0,
                           compress=FALSE, multithread=TRUE)
    }
  }
  
  ##Consider: Filter reads for tprK length and allow only 1 expected error for the entire read.
  for (count in c(1:length(filtEE1))) {
     track <- fastqFilter(nop[count], filtEE1[count], minLen=1400,maxLen=1800,
                          maxN=0, maxEE=1,
                          compress=FALSE, multithread=TRUE)
  }
  
  RAD_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.RAD.fasta",sep ='')
  RAD_files <- file.path(RAD_filenames)
  
  print("Constructing RAD files...")
  ## Build RAD files for each PacBio sample. This step takes forever!!!
  for (count in c(1:length(filt))) {
    to_rad_name <- paste(RAD_filenames[count])
    # Skips RAD step if files already exist, because it takes forever.
    if(file.exists(to_rad_name)) {
      print(paste(to_rad_name, " already exists. Skipping RAD step...", sep=""))
    } else{
      print("Setting up Julia...")
      #julia_command("using Pkg")
      #julia_command("Pkg.add(\"IJulia\")")
      #julia_command("Pkg.build(\"SpecialFunctions\")")
      #julia_command("Pkg.init(); Pkg.update(); Pkg.clone(\"https://github.com/MurrellGroup/NextGenSeqUtils.jl\"); using NextGenSeqUtils")
      #julia_command("Pkg.clone(\"https://github.com/MurrellGroup/DPMeansClustering.jl.git\")")
      #julia_command("Pkg.clone(\"https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git\"); using RobustAmpliconDenoising")
      
      #julia_command("Pkg.add(PackageSpec(name=\"NextGenSeqUtils\", rev= \"1.0\", url = \"https://github.com/MurrellGroup/NextGenSeqUtils.jl.git\"))")
      #julia_command("Pkg.add(PackageSpec(name=\"DPMeansClustering\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/DPMeansClustering.jl.git\"))")
      #julia_command("Pkg.add(PackageSpec(name=\"RobustAmpliconDenoising\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git\"))")
      
      #julia_command("using NextGenSeqUtils")
      #julia_command("using RobustAmpliconDenoising")
      
      
      #julia_readfastq <- paste("seqs, QVs, seq_names = read_fastq(\"",filt[count],'")',sep="")
      julia_readfastq <- paste("seqs, QVs, seq_names = read_fastq(\"",filtEE1[count],'")',sep="")
      
      julia_command(julia_readfastq)
      julia_command("templates,template_sizes,template_indices = denoise(seqs)")
      julia_writefasta <- paste("write_fasta(\"",RAD_files[count],'",templates,names = ["seqs$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])',sep="")
      julia_command(julia_writefasta)
    }
  }
  
  #RAD_files<- c(as.character(metadata$PacBio))
  ## RAD denoised files are written.  Let's get some frequencies of different variable regions
  RAD_files_nolines <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fasta",sep ='')
  RAD_files_fix <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fix.fasta",sep ='')
  #RAD_files_fix <-c(as.character(metadata$PacBio))
  # Fixes up the fastas so they wrap and don't have awkward new lines.
  # TODO: fix this section so it works. For some reason the pipeline currently runs without it? But probably should fix this anyway.
  awk_command <- paste("awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' < ",RAD_files," > ",RAD_files_nolines," ;")
  fix_firstline <- paste("tail -n+2 ",RAD_files_nolines," > ",RAD_files_fix)
  for (count in c(1:length(awk_command))) {
    system(awk_command[count])
    system(fix_firstline[count])
  }
  
  print("Making PacBio frequency files...")
  # Calls syph_r to make the final_data.csv for each PacBio sample.
  # TODO: right now this doesn't work if it's a repeat run because of previous .fastq files. Fix so you can rerun safely
  # without having to delete stuff when rerunning.
  for (num in c(1:length(RAD_files_fix))) {
    syphrPacBio_command <- paste("python3 ",syph_path," -i fasta -pacbio -d ./",sep='')
    system(syphrPacBio_command)
  }
  
  ## Make PacBio frequency comparison file
  PacBio_freq_files = list.files("./",pattern="*_final_data.csv",full.names=T)
  PacBio_freq_files_fullpath <- PacBio_freq_files
  PacBio_freq_files_fullpath = list.files(".",pattern="*_final_data.csv",full.names=T)
  compare_PacBio_df <- data.frame(Region=character(),Read=character())
  sample_names <- sort(sample_names)
  
  
  i=1
  
  
  # Renames the columns attaching PB and sample name to relative frequency and count so we can 
  # mash everything together into a big dataframe.
  for (i in 1:length(PacBio_freq_files)) {
    PacBioFreqtitle <- paste("PB_",sample_names[i],"_RelativeFreq",sep='')
    PacBioCounttitle <- paste("PB_",sample_names[i],"_Count",sep='')
    pacbiodf <- read.csv(PacBio_freq_files_fullpath[i],col.names = c("Region","Read",PacBioFreqtitle,PacBioCounttitle),check.names = FALSE)
    pacbiodf <- pacbiodf[order(pacbiodf$Region,-pacbiodf[[PacBioFreqtitle]]),][]
    compare_PacBio_df <- merge(compare_PacBio_df,pacbiodf,all=TRUE)
    # assign(PacBio_freq_files[i], read.csv(PacBio_freq_files_fullpath[i],col.names = c("Region","Read",PacBioFreqtitle,PacBioCounttitle)))
    # compare_PacBio_df <- merge(compare_PacBio_df,get(PacBio_freq_files[i]),all=TRUE)
  } 
} else {
  print("Illumina option specified. Skipping making PacBio frequency files...")
}

#setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_v6/chinese_data/trimmed')
if (opt$pacbio == FALSE) {
  print("Making Illumina frequency files...")
  Illumina_freq_path <- "Illumina/"
  Illumina_freq_path <- "./"

  # Creates frequency files for Illumina (final_data.csvs).
  syphrIllumina_command <- paste("python3 ",syph_path," -i fastq -illumina -d ",Illumina_freq_path,"/ ;",sep='')
  system(syphrIllumina_command)
  
  # Makes Illumina frequency comparison file
  Illumina_freq_files = list.files(Illumina_freq_path,pattern="*_final_data.csv")
  Illumina_freq_files_fullpath = list.files(Illumina_freq_path,pattern="*_final_data.csv",full.names=T)
  Illumina_freq_files_fullpath <- Illumina_freq_files
  compare_Illumina_df <- data.frame(Region=character(),Read=character())
  
  #svgIllumina_command <- paste("/usr/bin/python /Users/uwvirongs/Documents/Michelle/syphilis/syph_visualizer_single_svg.py ",Illumina_freq_files_fullpath,sep='')
  #for (num in length(svgIllumina_command)) system(svgIllumina_command[num])
  
  # Renames columns to attach Ill_ and sample name onto relative frequency and count columns so we can mash everything
  # together into a big dataframe.
  for (i in 1:length(Illumina_freq_files)) {
    IlluminaFreqtitle <- paste("Ill_",sample_names[i],"_RelativeFreq",sep='')
    IlluminaCounttitle <- paste("Ill_",sample_names[i],"_Count",sep='')
    Illuminadf <- read.csv(Illumina_freq_files_fullpath[i],col.names = c("Region","Read",IlluminaFreqtitle,IlluminaCounttitle),check.names = FALSE)
    Illuminadf <- Illuminadf[order(Illuminadf$Region,-Illuminadf[[IlluminaFreqtitle]]),][]
    compare_Illumina_df <- merge(compare_Illumina_df,Illuminadf,all=TRUE)
  } 
} else {
  print("Pacbio option specified. Skipping making Illumina frequency comparison files.")
}

#setwd("/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_five_billion_maybe_final/chinese_data")
setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_five_billion_maybe_final')
old_table <- read.table("allreads_illumina.csv", sep=',', header=TRUE)
allreads <- read.table("allreads_chinesetrimmed_filtered.csv",sep=',',header=TRUE)
#dunno <- cbind(old_table, allreads)
#colnames(old_table)
newtable <- full_join(old_table,allreads,by=c("Read","Region"))

old_table <- read.table("unchecked_Table_allAAfilt_fullORFs.tsv",header=T)
new_table <- read.table("Table_allAAfilt_fullORFs.tsv",sep="\t",header=T)
kms <- full_join(old_table,new_table,by=c("sample","sequences"))

#newtable2 <- newtable[!with(newtable,is.na(PB_AS12old_RelativeFreq)& is.na(Ill_AS12_RelativeFreq)& 
#                              is.na(PB_AS8old_RelativeFreq)& is.na(Ill_AS8_RelativeFreq)& 
#                              is.na(PB_AS12_RelativeFreq)& is.na(PB_AS8_RelativeFreq)),]
#newtable3 <- select(newtable2, -Region.y)
#colorder <- c("Region.x","Read","PB_AS12old_RelativeFreq","PB_AS12_RelativeFreq","Ill_AS12_RelativeFreq",
#              "PB_AS8old_RelativeFreq","PB_AS8_RelativeFreq","Ill_AS8_RelativeFreq")
#newtable3 <- newtable3[, colorder]


if (opt$pacbio) {
  allreads <- compare_PacBio_df
} else if (opt$illumina) {
  allreads <- compare_Illumina_df
} else {
  # Merges the PacBio and Illumina comparison tables into the all-important allreads.csv!
  print("Merging to allreads.csv...")
  allreads <- merge(compare_PacBio_df,compare_Illumina_df,all=T)
}
path<-"."
allreads_out <- paste(path,"/allreads.csv",sep='')
write.csv(newtable,file="allreads_combined.csv",row.names=FALSE,quote=FALSE)

getwd()
file.size("compare_pacbio_df.csv")
# allreads$`Illumina
# allreads$`IlluminaRelativeFreq_148B-tprK` <- as.numeric(allreads$`IlluminaRelativeFreq_148B-tprK`)
# allreads$`IlluminaRelativeFreq_148B2-tprK` <- as.numeric(allreads$`IlluminaRelativeFreq_148B2-tprK`)
# ggplot(allreads,aes(x=X148B_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + xlim(0,10) + ylim(0,10)
# 
# ggplot(allreads,aes(x=X148B_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + scale_y_log10(limits = c(0.05,100)) + scale_x_log10(limits = c(0.05,100)) #+ xlim(0.01,100) + ylim(0.01,100)
# 
# ggplot(allreads,aes(x=allreads$X148B2.tprK_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B2-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + xlim(0,10) + ylim(0,10)
