# Verifies that reads are contained in both Illumina preps at a mean rf of >0.1%.
# Also makes Illumina_squared files for each individual sample.
# Then validates PacBio files against Illumina files.

list.of.packages <- c("optparse", "dplyr", "tidyr", "tibble", "foreach","iterators","doParallel")
lapply(list.of.packages,library,character.only=T)

option_list <- list(make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"));

metadata <- read.table("metadata.csv", sep=',', header=TRUE)
metadata <- read.table(opt$metadata, sep=',', header=TRUE)
sample_names = as.character(metadata$SampleName)
pacbio_table = as.character(metadata$PacBio)
illumina1_table = as.character(metadata$Illumina1)
illumina2 = as.character(metadata$Illumina2)

setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_five_billion_maybe_final/bunnies')
setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_no_denoising')

num_cores <- 15
registerDoParallel(cores=num_cores)

# Changed from 0.2 to 0.1
rf_cutoff = 0.1

#foreach(i=1:length(illumina2)) %dopar% {
#  make_table_command <- paste("/Users/uwvirongs/miniconda3/bin/python3 /Users/uwvirongs/Documents/Michelle/tprk_files/illumina_check/syph_r_illumina_check.py -i fastq -d ./ -s ",illumina2[i],sep="")
#  system(make_table_command)
#}

illumina1_table <- paste(substr((illumina1_table),1,nchar(illumina1_table)-6),"_over5count_final_dna_data.csv",sep='')
illumina2_table <- paste(substr((illumina2),1,nchar(illumina2)-6),"_over5count_final_dna_data.csv",sep='')

## Comparing Illuminas ##
Illumina_freq_files = list.files('./',pattern="*_final_dna_data.csv")
Illumina_freq_files_fullpath <- Illumina_freq_files
compare_Illumina_df <- data.frame(Region=character(),Read=character())

#setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_five_billion_maybe_final/Illumina')
i=1
for (i in 1:length(sample_names)) {
  IlluminaFreqtitle <- paste("Ill_",sample_names[i],"_RelativeFreq",sep='')
  IlluminaCounttitle <- paste("Ill_",sample_names[i],"_Count",sep='')
  IlluminaFreq2title <- paste("Ill2_",sample_names[i],"_RelativeFreq",sep='')
  IlluminaCount2title <- paste("Ill2_",sample_names[i],"_Count",sep='')
  Illuminadf <- read.csv(illumina1_table[i],col.names = c("Region","Read",IlluminaFreqtitle,IlluminaCounttitle),check.names = FALSE)
  Illuminadf <- Illuminadf[order(Illuminadf$Region,-Illuminadf[[IlluminaFreqtitle]]),][]
  
  Illumina2df <- read.csv(illumina2_table[i],col.names = c("Region","Read",IlluminaFreq2title,IlluminaCount2title),check.names = FALSE)
  Illumina2df <- Illumina2df[order(Illumina2df$Region,-Illumina2df[[IlluminaFreq2title]]),][]
  Illumina_vs_Illumina_df <- merge(Illuminadf,Illumina2df,all=TRUE)
  
  write.csv(Illumina_vs_Illumina_df,file=paste(sample_names[i],"_prefiltered.csv",sep=""),row.names=FALSE,quote=FALSE)
  
  # Only take reads that are in both Illumina dataframes
  commondfIllumina <- filter(Illumina_vs_Illumina_df,!(is.na(Illumina_vs_Illumina_df[[IlluminaFreq2title]]) | is.na((Illumina_vs_Illumina_df[[IlluminaFreqtitle]]))))
  Illumina_vs_Illumina_name <- paste(sample_names[i],"_Illumina_squared.csv",sep="")
  write.csv(commondfIllumina, file=Illumina_vs_Illumina_name,row.names=FALSE,quote=FALSE)
  
  system(paste("Rscript /Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_v6_code/recalculate_frequency_amin.R -m ../metadata.csv -f ",
               Illumina_vs_Illumina_name))
  
  commondfIllumina <- read.csv(Illumina_vs_Illumina_name,sep=",")

  #common_mean_Illumina <- commondfIllumina %>% transmute(Read,Region,
  #                                                       rf = rowMeans(select(.,IlluminaFreqtitle,IlluminaFreq2title)),
  #                                                       count = rowMeans(select(.,IlluminaCounttitle,IlluminaCount2title)))
  
  # Use our second Illumina set of sample for counts and mean the rfs
  common_mean_Illumina <- commondfIllumina %>% transmute(Read,Region,
                                                         rf = rowMeans(select(.,IlluminaFreqtitle,IlluminaFreq2title)),
                                                         count = rowMeans(select(.,IlluminaCounttitle,IlluminaCount2title)))

  common_mean_Illumina <- common_mean_Illumina %>% filter(rf >= rf_cutoff)
  names(common_mean_Illumina)[names(common_mean_Illumina) == "rf"] <- IlluminaFreqtitle
  names(common_mean_Illumina)[names(common_mean_Illumina) == "count"] <- IlluminaCounttitle
  common_mean_Illumina <- common_mean_Illumina[c(2,1,3,4)]
  common_mean_Illumina_name <- paste(sample_names[i],"_validated_over5_count.csv",sep="")
  write.csv(common_mean_Illumina,file=common_mean_Illumina_name,row.names=FALSE,quote=FALSE)
  system(paste("Rscript /Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_v6_code/recalculate_frequency_amin.R -m ../metadata.csv -f ",
               common_mean_Illumina_name))
  
  compare_Illumina_df <- merge(compare_Illumina_df,common_mean_Illumina,all=TRUE)
}

write.csv(compare_Illumina_df,file="Illumina_validated_Illumina_allreads.csv",row.names=FALSE,quote=FALSE)
system(paste("Rscript /Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_v6_code/recalculate_frequency_amin.R -m ../metadata.csv -f Illumina_validated_Illumina_allreads.csv"))

#blah <- read.csv("AS10-amplicon-2_S48_L001_R1_001_over10count_final_dna_data.csv",sep=",",header=T)
#blah2 <- read.csv("AS10-amplicon_S31_L001_R1_001_over10count_final_dna_data.csv",sep=",",header=T)
#idk <- dplyr::full_join(blah, blah2, by = c("Region", "Read"))
#write.csv(idk, "as10_join.csv", row.names = F, quote=F)

## Comparing PacBio to Illumina^2 ##
for (i in 1:length(pacbio_table)) {
  make_table_command <- paste("/Users/uwvirongs/miniconda3/bin/python3 /Users/uwvirongs/Documents/Michelle/tprk_files/illumina_check/RAD_to_csv.py -i ",sample_names[i],"_validated_over5_count.csv -s ",pacbio_table[i],sep="")
  print(make_table_command)
  system(make_table_command)
}

#amin checking

setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/tprk_redo_five_billion_maybe_final/PacBio')

make_table_command <- paste("/Users/uwvirongs/miniconda3/bin/python3 /Users/uwvirongs/Documents/Michelle/tprk_files/illumina_check/RAD_to_csv-fuckaround.py -i AS22_Illumina_squared.csv -s AS22_m54179_190829_055030.Q20.noprimers.filtered.RAD.nolines.fix.fasta",sep="")
print(make_table_command)
system(make_table_command)
