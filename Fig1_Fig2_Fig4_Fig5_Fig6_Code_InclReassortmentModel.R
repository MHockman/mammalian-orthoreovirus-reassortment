library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readr)
library(gtable)
library(ggnewscale)
library(scales)
   
#Color Presets
Predictioncolor <- "#841c1c"

T3Dcolor <- "#220e1c"
T3DcolorLite <- "#872e6c"
T3DS208Pcolor <- "#548894"

Nococolor <- "#c8e07a"

T1Lcolor <- "#009cb9"
T1LcolorLite <- "#00ccf2"
T1LP208Scolor <- "#e78e10"
T1Lvarcolor <- "#df3522"

#Set basic plot parameters
BasePlot_legend <- ggplot(data = NULL) +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 40),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 40),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = c(0.2, 0.9),
        legend.background=element_blank(),
        legend.key = element_blank(),
        plot.margin = margin(10, 10, 50, 50)) 

#Get the Data
AllSerotypes_ReassortmentData_NarrowGate <- read_csv("data/AllSerotypes_ReassortmentData_NarrowGate.csv")

T1L_Growth_Curve_CSV <- read_csv("data/T1L Growth Curve CSV.csv")
T1L_Growth <- melt(T1L_Growth_Curve_CSV, id.vars = c("Time", "Virus"))

#Plot Figure 1###########################################################################################
Fig_1A_T1LGrowthCurve <- BasePlot_legend +
  geom_smooth(data = T1L_Growth, aes(x = Time, y = log10(value), color = Virus)) +
  geom_point(data = T1L_Growth, aes(x = Time, y = log10(value), color = Virus)) +
  scale_color_manual(values = c(T1Lvarcolor, T1Lcolor)) +
  labs(x = "Time (h)",
       y = expression('Titer (log'[10]*'(PFU/mL))'))
Fig_1A_T1LGrowthCurve

ggsave("1A_T1LGrowthCurve_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)

#Plot Figure 2##########################################################################################

#Generate model predictions of reassortment
moilist <- seq(-1, 2.5, by = 0.05)
moilist <- 10^moilist 

reassortment.function <- function(mu){
  
  limit <- qpois(0.9999999999, lambda = mu)
  reassortment <- 0
  numsum <- 0
  denomsum <- 0
  totvir <- 0
  
  for(i in 0:limit)
  {
    pi <- dpois(i, lambda = 0.5*mu)
    
    for (j in 0:limit){
      
      pj <- dpois(j, lambda = 0.5*mu)
      
      totvir <- i + j
      
      if(totvir != 0){
        
        c <- k * totvir
        
        denom <- pi*pj*c
        numer <- denom*(1-((i/(i+j))^10) - ((j/(i+j))^10))
        
        numsum <- numsum + numer
        denomsum <- denomsum + denom
        
      }
      
    }
    
  }
  
  reassortment <- numsum/denomsum
  return(reassortment)
}

k <- 11
T3Dreassortment_prediction <- sapply(X = moilist, FUN = reassortment.function)
T3Dreassortment_prediction.df <- data.frame(reassortment = T3Dreassortment_prediction,
                                            MOI = moilist)

k <- 93
T1Lreassortment_prediction <- sapply(X = moilist, FUN = reassortment.function)
T1Lreassortment.df_prediction <- data.frame(reassortment = T1Lreassortment_prediction,
                                            MOI = moilist)

Figure_2A_T1LReassortment <- BasePlot_legend +
  geom_line(data = T1Lreassortment.df_prediction, size = 3, color = Predictioncolor, aes(x = MOI, y = reassortment, lty = 'Prediction')) +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% "T1L"), color = T1Lcolor, size = 6, aes(x = MOI, y = PercR)) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.title = element_blank())+
  labs(title = "T1L",
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_2A_T1LReassortment  

ggsave("2A_T1LReassortment_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)


Figure_2B_T3DReassortment <- BasePlot_legend +
  geom_line(data = T3Dreassortment_prediction.df, size = 3, color = Predictioncolor, aes(x = MOI, y = reassortment, lty = 'Prediction')) +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% "T3D"), color = T3Dcolor, size = 6, aes(x = MOI, y = PercR)) +
  theme(legend.title = element_blank()) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "T3D", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_2B_T3DReassortment

ggsave("2B_T3DReassortment_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)  

#Plot Figure 4#####################################################################################

Figure_4C_T1LReassortment <- BasePlot_legend +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% c("T1L", "T1LP208S")), size = 6, aes(x = MOI, y = PercR, color = Virus)) +
  scale_color_manual(values = c(T1Lcolor, T1LP208Scolor)) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "T1L", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_4C_T1LReassortment

ggsave("4C_T1LReassortmentMorphology_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)  

Figure_4D_T3DReassortment <- BasePlot_legend +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% c("T3D", "T3DS208P")), size = 6, aes(x = MOI, y = PercR, color = Virus))+
  scale_color_manual(values = c(T3Dcolor, T3DS208Pcolor)) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "T3D", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_4D_T3DReassortment

ggsave("4D_T3DReassortmentMorphology_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)  

#Data for 4E
AllSerotypes_PercI_Replicates_Narrow <- read_csv("data/AllSerotypes_PercI_Replicates_Narrow.csv")
T3DPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T3D")
T3DS208PPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T3DS208P")
T1LPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T1L")
T1LP208SPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T1LP208S")


Ratios.T3D.df <- data.frame(Serotype = "T3D",
                            MOI = T3DPerci$MOI,
                            Ratio = T3DS208PPerci$PercI/T3DPerci$PercI)
Ratios.T1L.df <- data.frame(Serotype = "T1L",
                            MOI = T1LPerci$MOI,
                            Ratio = T1LP208SPerci$PercI/T1LPerci$PercI)

Ratios <- bind_rows(Ratios.T3D.df, Ratios.T1L.df)
Ratios_average <- aggregate(Ratio ~ Serotype, Ratios, mean)

Figure_4E_InfectionRatios <- BasePlot_legend +
  geom_col(data = Ratios_average, aes(Serotype, Ratio, fill = Serotype)) +
  geom_point(data = Ratios, size = 3, aes(x = Serotype, y = Ratio, color = Serotype)) +
  scale_color_manual(values = c(T3DcolorLite, T1LcolorLite)) +
  scale_fill_manual(values = c(T3Dcolor, T1Lcolor)) +
  geom_point(color = "black", pch = 21, size = 3.3, data = Ratios, aes(x = Serotype, y = Ratio)) +
  theme(legend.position = "None") +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(x = "\nSerotype",
       y = "Fold change in \nproportion infected cells\n")
Figure_4E_InfectionRatios

ggsave("4E_InfectionRatios_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)  

Figure_4F_AllSerotypeReassortment <- BasePlot_legend +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% c("T1L", "T1LP208S", "T3D", "T3DS208P")), size = 6, aes(x = PercI, y = PercR, color = Virus)) +
  scale_color_manual(values = c(T1Lcolor, T1LP208Scolor, T3Dcolor, T3DS208Pcolor)) +
  theme(legend.position = c(0.8, 0.2)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "\nProportion infected cells",
       y = "Proportion reassortant\n")
Figure_4F_AllSerotypeReassortment

ggsave("4F_ReassortmentvsPercI_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)

#Plot Figure 5##################################################################################

#Get Data for Figure 5
AllSerotypes_PercI_Replicates_Narrow <- read_csv("data/AllSerotypes_PercI_Replicates_Narrow.csv")

AllSerotypes_Averages <- aggregate(AllSerotypes_PercI_Replicates_Narrow$PercI, 
                                   by = list(AllSerotypes_PercI_Replicates_Narrow$MOI, AllSerotypes_PercI_Replicates_Narrow$Virus), 
                                   mean, na.rm = TRUE)

T3DPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T3D")
T3DS208PPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T3DS208P")
T1LPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T1L")
T1LP208SPerci <- subset(AllSerotypes_PercI_Replicates_Narrow, Virus == "T1LP208S")

#Get predicted infection levels based on Poisson
moilist <- seq(-1, 2.5, by = 0.05)
moilist <- 10^moilist 
percipredicted <- 1 - exp(-moilist)
percipredicted_adjusted <- 1 - exp(-0.4 * moilist)
prediction.df <- data.frame(moi = moilist,
                            perci = percipredicted,
                            perciadj = percipredicted_adjusted)

prediction.df <- melt(prediction.df, id.vars = "moi")

#Plot the Data
Figure_5A_PercentInfectionData <- BasePlot_legend +
  geom_line(data = AllSerotypes_Averages, size = 1.5, aes(x = Group.1, y = x, color = Group.2)) +
  scale_color_manual(values = c(T1Lcolor, T1LP208Scolor, T3Dcolor, T3DS208Pcolor)) +
  new_scale_color() +
  geom_point(data = AllSerotypes_PercI_Replicates_Narrow, size = 3, aes(x = MOI, y = PercI, color = Virus)) +
  scale_color_manual(values = c(T1Lcolor, T1LP208Scolor, T3Dcolor, T3DS208Pcolor), guide = 'none') +
  geom_point(color = "black", pch = 21, size = 3.3, data = AllSerotypes_PercI_Replicates_Narrow, aes(x = MOI, y = PercI)) +
  new_scale_color() +
  geom_line(size = 1.5, data = prediction.df, aes(x = moi, y = value, color = variable)) +
  scale_color_manual(labels = c("Prediction\n", "Adjusted\nPrediction"), values = c("#AB2208", "#F76F56")) +
  theme(legend.text = element_text(size = 25),
        legend.position = c(0.15, 0.75),)+
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion infected cells\n",
       color = "Virus")  
Figure_5A_PercentInfectionData  

ggsave("5A_PercIvsPrediction_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)

#Get the Data for 5B
AllSerotypes_Invsout_All_Replicates_NarrowGate <- read_csv("data/AllSerotypes_Invsout_All Replicates_NarrowGate.csv")

Figure_5B_InputvsOutput <- BasePlot_legend +
  geom_point(data = AllSerotypes_Invsout_All_Replicates_NarrowGate, aes(x = MOI, y = Output, color = Virus), size = 3) +
  scale_color_manual(values = c(T1Lcolor, T1LP208Scolor, T3Dcolor, T3DS208Pcolor)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(.x))) +
  labs(title = "", 
       x = expression("MOI (log"[10]*"PFU/cell)"),
       y = expression("Virus output (log"[10]*"PFU)"),
       #y = expression("log"[10]*"(Virus output (PFU))"),
       color = "Virus")
Figure_5B_InputvsOutput

ggsave("5B_InputvsOutput_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)

#Plot Figure 6######################################################################################

Figure_6C_T3DvsT3DNocoReassortment <- BasePlot_legend +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% c("T3D", "T3D Noco")), size = 6, aes(x = MOI, y = PercR, color = Virus)) +
  scale_color_manual(values = c(T3Dcolor, Nococolor)) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "T3D", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_6C_T3DvsT3DNocoReassortment     

ggsave("6C_T3DvsNoco_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)

Figure_6D_T3DS208PvsT3DS208PNocoReassortment <- BasePlot_legend +
  geom_point(data = subset(AllSerotypes_ReassortmentData_NarrowGate, Virus %in% c("T3DS208P", "T3DS208P Noco")), size = 6, aes(x = MOI, y = PercR, color = Virus)) +
  scale_color_manual(values = c(T3DS208Pcolor, Nococolor)) +
  scale_x_log10(limits = c(NA, 320)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = c(0.3, 0.9)) +
  labs(title = "T3DS208P", 
       x = "\nMOI (PFU/cell)",
       y = "Proportion reassortant\n")
Figure_6D_T3DS208PvsT3DS208PNocoReassortment  

ggsave("6D_T3DS208PvsT3DS208PNocoReassortmen_ReviewerRevision.tiff", width = 10, height = 10, dpi = 1200)
