# Load libraries
library(ggbio)
library(rtracklayer)
library(biovizBase)
library(BSgenome.Drerio.UCSC.danRer10)
library(dplyr)
library(reshape)

# Set directory
setwd("~/Documents/GitHub/Met_zebrafish_germline/Data")

############################################################
# Grand plot                                               #
############################################################

# Load data and transform variables
data_Fig4A <- read.csv2("probes_1Mb.txt",sep = "\t",stringsAsFactors = FALSE,dec=".")
data_Fig4A$norm_21d_egfp <- data_Fig4A$X21d_egfp/(sum(data_Fig4A$X21d_egfp)/nrow(data_Fig4A))
data_Fig4A$norm_25d_egfp <- data_Fig4A$X25d_egfp/(sum(data_Fig4A$X25d_egfp)/nrow(data_Fig4A))
data_Fig4A$norm_28d_egfp <- data_Fig4A$X28d_egfp/(sum(data_Fig4A$X28d_egfp)/nrow(data_Fig4A))

# Create genomic ranges table
gr.table <- transformDfToGr(data_Fig4A, seqnames = "Chromosome", start = "Start",end="End")
head(gr.table)

# Calculate seq lenghts
seqlengths(gr.table)

danRer10IdeogramCyto <- getIdeogram("danRer10")
seqlengths(gr.table) <- as.numeric(seqlengths(danRer10IdeogramCyto)[1:26])

seqlengths(gr.table)

# Plot GrandLinear GSC 21d
plotGrandLinear(gr.table, geom = "area",aes(y = log10(norm_21d_egfp+1))) + ylab("") + ylim(0,3) + theme_gray()+
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank(),legend.position="none")+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Plot GrandLinear GSC 25d
plotGrandLinear(gr.table, geom = "area",aes(y = log10(norm_25d_egfp+1))) + ylab("") + ylim(0,3) + theme_gray()+
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank(),legend.position="none")+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# Plot GrandLinear GSC 28d
plotGrandLinear(gr.table, geom = "area",aes(y = log10(norm_28d_egfp+1))) + ylab("") + ylim(0,3) + theme_gray()+
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank(),legend.position="none")+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

# 1500 x 250

############################################################
# Ideogram                                                 #
############################################################
danRer10IdeogramCyto <- getIdeogram("danRer10", cytoband = FALSE)

plotIdeogram(danRer10IdeogramCyto,"chr4",xlabel = TRUE,zoom.region = c(76310000, 76335000),
             color = "red",fill = "red",alpha = 0.5,size = 2,aspect.ratio = 1/10)

# 750 x 125

############################################################
# Coverage plot                                            #
############################################################

bg <- BSgenome.Drerio.UCSC.danRer10
gr <- GRanges(seqnames = "4", ranges = IRanges(start = 76310000, end = 76335000))

CTOT.bam <- "CTOT.bam"
autoplot(CTOT.bam, bsgenome = bg,which = gr, geom = "area") + xlim(gr) + ylim(c(0,100))

CTOB.bam <- "CTOB.bam"
autoplot(CTOB.bam, bsgenome = bg,which = gr, geom = "area") + xlim(gr) + ylim(c(0,100))

# 700 x 200

############################################################
# Wilson data fem rDNA                                     #
############################################################

# Load data and transform variables
Wilson_data <- read.csv2("Wilson_data.csv",sep = ",",stringsAsFactors = FALSE,dec=".")
str(Wilson_data)

p <- ggplot(Wilson_data, aes(Coordinates,P_value))
p + geom_point(aes(colour = factor(Strain)),size = 3) +
  geom_hline(yintercept=-log(0.01), linetype="dashed", color = "black") +
  xlim(73900000, 76625712)

# 700 x 550

############################################################
# Relation amplification - methylation                     #
############################################################

# Load data and transform variables
data_fig_5b <- read.csv2("data_fem_rDNA.csv",sep = ",",stringsAsFactors = FALSE,dec=",")

data_cleaned_fig_5b = data_fig_5b[data_fig_5b$CG_calls>10,]

data_cleaned_fig_5b$Amplification
p <- ggplot(data_cleaned_fig_5b, aes(Amplification,CG_met))

p + geom_point(aes(colour = factor(Label)),size = 6) +
  scale_color_manual(values=c("#808080","#FFB6C1","#add8e6","#FF00FF",
                              "#ff0000","#0000ff","#00AA00"))+
  geom_vline(xintercept = 1, linetype= "dashed", color = "black")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.5)) +
  scale_x_log10()

##########################################################
# Summary amplification                                    #
############################################################

# Load data and transform variables
fig_6b <- read.csv2("Data_summary_B.csv",sep = ",",stringsAsFactors = FALSE,dec=".")
str(fig_6b)

p <- ggplot(fig_6b, aes(Day,Amplification))
p + geom_point(aes(colour = factor(Sex)),size = 8) +
  scale_color_manual(values=c("#ff0000","#0000ff","#808080"))+
  xlim(0,61) + ylim(0,25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.5))
