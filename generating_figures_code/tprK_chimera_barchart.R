library(ggplot2)

setwd('/Users/uwvirongs/Documents/Michelle/tprk_files/new_tprk/demuxed/20200820_pacbio/chimera_bar_graph/')
tprk <- read.csv("tprk-chimera_rename.csv", header=TRUE, stringsAsFactors = F)

#M <- c("AS10,MI04,Chimera1,Chimera2")
#tprk$Clone <- factor(tprk$Clone, levels = M)

plot <- ggplot(tprk, aes(x=Replicate, y=Percent, fill=reorder(Clone,Percent))) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 0, vjust=5), axis.title.x = element_text(vjust = -0.2), axis.ticks.x = element_blank()) + 
  scale_fill_brewer(palette="Spectral", name = "Clone") + guides(fill = guide_legend(reverse = TRUE)) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) 
ggsave("chimera_barchart.pdf",plot=plot,width=5,height=3,units="in")
