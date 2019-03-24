# Load libraries
library(ggbio)
library(rtracklayer)
library(biovizBase)
library(BSgenome.Drerio.UCSC.GRCz11)
library(dplyr)
library(reshape)
library(ggplot2)

# Set directory
setwd("~/Documents/GitHub/Met_zebrafish_germline")

############################################################
# Grand plot                                               #
############################################################

# Load data and transform variables
data_1Mb_probes <- read.csv2("Data/Probes_1Mb/Probes_1Mb.txt",sep = "\t",stringsAsFactors = FALSE,dec=".")
data_1Mb_probes$norm_21d_egfp <- data_1Mb_probes$X21d_egfp/(sum(data_1Mb_probes$X21d_egfp)/nrow(data_1Mb_probes))
data_1Mb_probes$norm_25d_egfp <- data_1Mb_probes$X25d_egfp/(sum(data_1Mb_probes$X25d_egfp)/nrow(data_1Mb_probes))
data_1Mb_probes$norm_28d_egfp <- data_1Mb_probes$X28d_egfp/(sum(data_1Mb_probes$X28d_egfp)/nrow(data_1Mb_probes))

# Create genomic ranges table
gr.table <- transformDfToGr(data_1Mb_probes, seqnames = "Chromosome", start = "Start",end="End")
head(gr.table)

# Calculate seq lenghts
seqlengths(gr.table)

danRer11IdeogramCyto <- getIdeogram("danRer11")
seqlengths(gr.table) <- as.numeric(seqlengths(danRer11IdeogramCyto)[1:26])

seqlengths(gr.table)

# Plot GrandLinear germline 21d
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

# Plot GrandLinear germline 25d
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
danRer11IdeogramCyto <- getIdeogram("danRer11", cytoband = FALSE)

plotIdeogram(danRer11IdeogramCyto,"chr4",xlabel = TRUE,zoom.region = c(start = 77545000, end = 77570000),
             color = "red",fill = "red",alpha = 0.5,size = 2,aspect.ratio = 1/10)

# 750 x 125

############################################################
# Coverage plot                                            #
############################################################

bg <- BSgenome.Drerio.UCSC.GRCz11
gr <- GRanges(seqnames = "4", ranges = IRanges(start = 77545000, end = 77570000))

CTOT.bam <- "Data/Mapping/CTOT.bam"

# Limits
#gr <- GRanges(seqnames = "4", ranges = IRanges(start = 77549888, end = 77549895))
#gr <- GRanges(seqnames = "4", ranges = IRanges(start = 77567273, end = 77567279))

autoplot(CTOT.bam, bsgenome = bg,which = gr, geom = "area") + xlim(gr) + ylim(c(0,100))

CTOB.bam <- "Data/Mapping/CTOB.bam"

# Limits
#gr <- GRanges(seqnames = "4", ranges = IRanges(start = 77549894, end = 77549900))
#gr <- GRanges(seqnames = "4", ranges = IRanges(start = 77567273, end = 77567279))

autoplot(CTOB.bam, bsgenome = bg,which = gr, geom = "area") + xlim(gr) + ylim(c(0,100))

# 700 x 200

############################################################
# Wilson data fem rDNA                                     #
############################################################

# Load data and transform variables
Wilson_data <- read.csv2("Data/Wilson_data.csv",sep = ",",stringsAsFactors = FALSE,dec=".")

p <- ggplot(Wilson_data, aes(Coordinates,P_value))
p + geom_point(aes(colour = factor(Strain),shape = factor(Strain)),size = 3) +
  scale_shape_manual(values=c(rep(19,2),108,rep(19,3))) +
  geom_hline(yintercept=-log(0.01), linetype="dashed", color = "black") +
  xlim(75000000,78090000)

# 700 x 550

############################################################
# Relation amplification - methylation                     #
############################################################

# Load data and transform variables
data_fem_rDNA <- read.csv("Data/data_fem_rDNA.csv",sep = ",",stringsAsFactors = FALSE)

data_cleaned_fem_rDNA = data_fem_rDNA[data_fem_rDNA$n_calls_fem_rDNA>10,]

p <- ggplot(data_cleaned_fem_rDNA, aes(per_reads_fem_rDNA,met_rDNA))

p + geom_point(aes(colour = factor(Label)),size = 6) +
  scale_color_manual(values=c("#00AA00","#FF00FF","#808080","#FFB6C1",
                              "#ff0000","#add8e6","#0000ff"))+
  geom_vline(xintercept = 1, linetype= "dashed", color = "black")+
  geom_hline(yintercept = 35, linetype= "dashed", color = "black")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.5)) +
  scale_x_log10()

data_cleaned_fem_rDNA_germline = data_cleaned_fem_rDNA[data_cleaned_fem_rDNA$type=="egfp",]
str(data_cleaned_fem_rDNA_germline)

sum(data_cleaned_fem_rDNA_germline$per_reads_fem_rDNA>1)

# 1000 x 800 (937.5 x 748.75) -> 1344.086 x 1073.477

# Group 1

group_1 <- data_cleaned_fem_rDNA[data_cleaned_fem_rDNA$met_rDNA > 35,]
mean(group_1$met_rDNA)
sd(group_1$met_rDNA)

# Group 2

group_2 <- data_cleaned_fem_rDNA[(data_cleaned_fem_rDNA$met_rDNA < 35)&
                                   (data_cleaned_fem_rDNA$per_reads_fem_rDNA < 1),]
mean(group_2$met_rDNA)
sd(group_2$met_rDNA)
min(group_2$per_reads_fem_rDNA)
max(group_2$per_reads_fem_rDNA)

# Group 3

group_3 <- data_cleaned_fem_rDNA[data_cleaned_fem_rDNA$per_reads_fem_rDNA > 1,]
mean(group_3$met_rDNA)
sd(group_3$met_rDNA)

############################################################
# Summary amplification                                    #
############################################################

# Load data and transform variables

data_fem_rDNA <- read.csv("Data/data_fem_rDNA.csv",sep = ",",stringsAsFactors = FALSE)

data_fem_rDNA_germline = data_fem_rDNA[data_fem_rDNA$type == "egfp",]

p <- ggplot(data_fem_rDNA_germline, aes(Day,per_reads_fem_rDNA ))

p + geom_point(aes(colour = factor(sex)),size = 4) +
  scale_color_manual(values=c("#ff0000","#0000ff","#00AA00"))+
  xlim(0,65) + ylim(0,25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.5))

# 1200 x 600

############################################################
# Simulation Low coverage WGBS                             #
############################################################

# Figure simulation

# Data 1 - Sperm
data_sperm <- read.csv("Data/Sampling/sperm_simulation.csv",sep = ",", stringsAsFactors = FALSE)
data_sperm$Total_Cs_CpG = data_sperm$Total_CsM_CpG + data_sperm$Total_CsU_CpG
data_sperm$CsM_CpG = data_sperm$Total_CsM_CpG * 100 / data_sperm$Total_Cs_CpG

data_sperm_clean = data.frame("ID" = data_sperm$Sample, "Pred_Cs" = data_sperm$Group,
                              "Total_Cs_CpG" = data_sperm$Total_Cs_CpG,
                              "CsM_CpG" = data_sperm$CsM_CpG)

data_sperm_summary <- as.data.frame(data_sperm_clean %>% group_by(Pred_Cs) %>%
                                      summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                                quantile_99_95 = quantile(CsM_CpG, c(.9995)),
                                                quantile_00_05 = quantile(CsM_CpG, c(.0005)),
                                                dif = quantile_99_95 - quantile_00_05,
                                                ME = dif / 2
                                      ))

data_sperm_summary

# Data 2 - Muscle
data_muscle <- read.csv("Data/Sampling/muscle_simulation.csv",sep = ",", stringsAsFactors = FALSE)
data_muscle$Total_Cs_CpG = data_muscle$Total_CsM_CpG + data_muscle$Total_CsU_CpG
data_muscle$CsM_CpG = data_muscle$Total_CsM_CpG * 100 / data_muscle$Total_Cs_CpG

data_muscle_clean = data.frame("ID" = data_muscle$Sample, "Pred_Cs" = data_muscle$Group,
                               "Total_Cs_CpG" = data_muscle$Total_Cs_CpG,
                               "CsM_CpG" = data_muscle$CsM_CpG)

data_muscle_summary <- as.data.frame(data_muscle_clean %>% group_by(Pred_Cs) %>%
                                       summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                                 quantile_99_95 = quantile(CsM_CpG, c(.9995)),
                                                 quantile_00_05 = quantile(CsM_CpG, c(.0005)),
                                                 dif = quantile_99_95 - quantile_00_05,
                                                 ME = dif / 2
                                       ))

data_muscle_summary


data_summary <- rbind(data_sperm_summary,data_muscle_summary)
data_summary["Tissue"] = c(rep("Sperm",12),rep("Muscle",12))

# Analyse data
data_all = rbind(data_sperm_clean, data_muscle_clean)
data_all["Tissue"] = c(rep("Sperm",12000),rep("Muscle",12000))

x_1 = c(data_sperm_summary$Pred_Cs, rev(data_sperm_summary$Pred_Cs))
y_1 = c(data_sperm_summary$quantile_99_95, rev(data_sperm_summary$quantile_00_05))

pol_data_sperm = data.frame("x" = x_1, "y"= y_1, "group" = c(rep("low",12),rep("high",12)))

x_2 = c(data_muscle_summary$Pred_Cs, rev(data_muscle_summary$Pred_Cs))
y_2 = c(data_muscle_summary$quantile_99_95, rev(data_muscle_summary$quantile_00_05))

pol_data_muscle = data.frame("x" = x_2, "y"= y_2, "group" = c(rep("low",12),rep("high",12)))

pol_all = rbind(pol_data_sperm, pol_data_muscle)
pol_all["Tissue"] = c(rep("Sperm",24),rep("Muscle",24))

# Plot simulation
plot_A = ggplot(data_all, aes(x = Total_Cs_CpG, y = CsM_CpG, color = Tissue)) +
  geom_point(alpha = 0.5) +
  xlab("Number of Cs in CG context") + ylab("CG methylation (%)") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    axis.line = element_line(size = 0.5)
  ) + 
  scale_color_manual(values=c('#F8766D','#619CFF')) +
  geom_polygon(data = pol_all, aes(x = x, y = y, fill = Tissue, color = NA),alpha = 0.2) +
  scale_fill_manual(values=c('#F8766D','#619CFF')) +
  ylim(50,100) + 
  geom_hline(yintercept=94.34,linetype="dashed", color="blue") + 
  geom_hline(yintercept=81.12,linetype="dashed", color="red")

# 1100 x 1000
plot_A

# Construct model
x1 = data_sperm_summary$Pred_Cs
y1 = data_sperm_summary$ME
xy1 = data.frame("x" = x1, "y" = y1)

mod_new1 <- lm(y1 ~ I(1/sqrt(x1))+0, xy1)
mod_new1$coefficients/100

x2 = data_muscle_summary$Pred_Cs
y2 = data_muscle_summary$ME
xy2 = data.frame("x" = x2, "y" = y2)

mod_new2 <- lm(y2 ~ I(1/sqrt(x2))+0, xy2)
mod_new2$coefficients/100

options(digits=4)

x_mod_1 <- seq(0,30000,10)
y_mod_1 <- c()
for (i in seq(0,30000,10)){
  y_temp = 1.207/(sqrt(i))*100
  y_mod_1 = c(y_mod_1,y_temp)
}

x_mod_2 <- seq(0,30000,10)
y_mod_2 <- c()
for (i in seq(0,30000,10)){
  y_temp = 2.109/(sqrt(i))*100
  y_mod_2 = c(y_mod_2,y_temp)
}

mod_1 = data.frame("X_mod"=x_mod_1,"Y_mod"=y_mod_1)
mod_2 = data.frame("X_mod"=x_mod_2,"Y_mod"=y_mod_2)

mod_all <- rbind(mod_1,mod_2)
mod_all["Tissue"] = c(rep("Sperm",3001),rep("Muscle",3001))

plot_B = ggplot(data_summary, aes(x = Pred_Cs, y = ME, color = Tissue)) +
  geom_point(alpha = 0.9, size = 4) +
  geom_line(data = mod_all, aes(x = X_mod, y = Y_mod, color = Tissue)) + 
  xlab("Number of Cs in CG context") + ylab("Margin of error (± %) 99% CI") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    axis.line = element_line(size = 0.5)
  ) +
  geom_hline(yintercept=1,linetype="dashed", color="blue", size = 0.75) +
  geom_hline(yintercept=2.5,linetype="dashed", color="green", size = 0.75) +
  geom_hline(yintercept=5,linetype="dashed", color="red",size = 0.75) +
  geom_vline(xintercept=10000,linetype="dashed", color="black",size = 0.5) +
  ylim(0,30) +
  scale_color_manual(values=c('#F8766D','#619CFF'))

# 1100 x 1000
plot_B
