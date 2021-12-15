# Meg Reo Pairwise Association Analysis

require(MASS)
require(tidyverse)
require(lmerTest)
require(robustlmm)
require(latex2exp)
require(ggthemes)
require(ggridges)
require(reshape2)

theme_set(theme_grey())

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
Proj.Home = "" #Set Home Location here
Data.Dir = file.path(Proj.Home,"data")
setwd(Data.Dir)

MRH3 = read.csv(file = file.path(Proj.Home,"data", "T3D_Reassortment_PairwiseFormat.csv"), header = TRUE, fileEncoding = "UTF-8-BOM") %>%
  na.omit() %>%
  mutate(L1 = L1 %>% as.numeric,
         L2 = L2 %>% as.numeric,
         L3 = L3 %>% as.numeric,
         M1 = M1 %>% as.numeric,
         M2 = M2 %>% as.numeric,
         M3 = M3 %>% as.numeric,
         S1 = S1 %>% as.numeric,
         S2 = S2 %>% as.numeric,
         S3 = S3 %>% as.numeric,
         S4 = S4 %>% as.numeric)

L1.Index = which(colnames(MRH3) == "L1")

Reps = table(MRH3$Rep) %>% names 
MOIs = table(MRH3$MOI) %>% names

#Generate empty data frame with desired columns
Corr = MRH3 %>% filter(Rep == 1, 
                       MOI == 100) %>%
  select(L1:S4) %>% cor %>% melt %>%
  mutate()

Corr = Corr[0,] #Isolate headers

for (rep in 1:length(Reps)) {     # Iterate through Replicates
  for (moi in 1:length(MOIs)) {   # Iterate through MOIs
    Corr = rbind(Corr, #Append data from current MOI and Replicate to "Corr" data frame
                 MRH3 %>% filter(Rep == Reps[rep],
                                 MOI == MOIs[moi]) %>%
                   select(L1:S4) %>% #Select the columns containing presence/absence flags for each segment
                   cor %>% #Calculate correlations between each segment pair
                   melt %>% # Convert correlation matrix to long form (Segment 1, Segment 2, Correlation, and MOI)
                   mutate(MOI = MOIs[moi],
                          Rep = Reps[rep]))
  }
}

Reo.Corr = Corr %>% 
  na.omit() %>%
  rename(X = Var1,
         Y = Var2,
         r = value) %>%
  mutate(MOI = MOI %>% as.numeric,
         r2 = r * r)

Corr = Corr[0,] #Empty data frame, retain only headers

for (moi in 1:length(MOIs)) { # Iterate through MOIs, but data from all replicates are pooled
  Corr = rbind(Corr,
               MRH3 %>% filter(MOI == MOIs[moi]) %>%
                 select(L1:S4) %>% 
                 cor %>%
                 melt %>%
                 mutate(MOI = MOIs[moi]))
}

Pooled.df = Corr %>% rename(X = Var1,
                            Y = Var2,
                            r = value) %>%
  mutate(MOI = MOI %>% as.numeric,
         r2 = r * r)

Pooled.df2 = Pooled.df %>% filter(X != Y) %>% na.omit()

# Export data ----
write.csv(Reo.Corr, file = file.path(Proj.Home,"data", "T3D_Segment_Associations.csv"), row.names = FALSE)
write.csv(Pooled.df2, file = file.path(Proj.Home,"data","T3D_Pooled_Associations.csv"), row.names = FALSE)

# Base plots ----
Base.Plot = ggplot() +
  theme(text=element_text(size=14,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=0,vjust=0,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

Stair.Plot = ggplot() +
  theme(text=element_text(size=10,face="bold"),
        plot.title=element_text(size=rel(2),face="bold"),
        legend.title.align=3,
        legend.title=element_text(size=rel(1.5),face="bold",vjust=0),
        legend.text=element_text(size=rel(1.2),face="bold",vjust=0),
        legend.key.width=unit(2,"cm"),
        axis.text.x = element_text(size=rel(1),vjust = .5,hjust=.5,face="bold",color = "black", margin=margin(5,0,0,0), angle=90),
        axis.text.y = element_text(size=rel(1),face="bold",color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line=element_blank(),
        strip.text.x = element_text(size=rel(1),angle=0,face="bold"),
        strip.text.y = element_text(size=rel(1),angle=0,face="bold"),
        strip.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank())

# Stair Plot ----

Stair.Plot +
  geom_tile(data = Pooled.df,
            aes(x = X,
                y = Y,
                fill = r2)) +
  facet_wrap(~MOI) +
  coord_equal() +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL,
       y = NULL,
       fill = TeX("\\textbf{$\\r^2$}")) +
  theme(legend.position = "NA")

ggsave(file = "T3D_Association_Stair_Plot.tiff",
       width = 5,
       height = 4,
       dpi = 1200)

Stair.Plot +
  geom_tile(data = Reo.Corr,
            aes(x = X,
                y = Y,
                fill = r2)) +
  facet_grid(Rep~MOI) +
  coord_equal() +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL,
       y = NULL,
       fill = TeX("\\textbf{$\\r^2$}")) +
  theme(legend.position = "NA")

ggsave(file = "T3D_Association_Stair_Plot_Individual.tiff",
       width = 8,
       height = 4,
       dpi = 1200)

Legend.Plot = Stair.Plot +
  geom_tile(data = Pooled.df,
            aes(x = X,
                y = Y,
                fill = r2)) +
  facet_wrap(~MOI) +
  coord_equal() +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL,
       y = NULL,
       fill = TeX("\\textbf{$\\r^2$}")) +
  theme(legend.position = "bottom")

ggsave(g_legend(Legend.Plot),
       file = "T3D_Association_Stair_Plot_Legend.tiff",
       width = 5,
       height = 1,
       dpi = 1200)

# Ridge Plot ----
Base.Plot +
  geom_density_ridges(data = Pooled.df2,
                      aes(y = factor(MOI),
                          x = r2),
                      alpha = 0.8,
                      fill = "navy") +
  labs(x = TeX("\\textbf{Pairwise association ($\\r^2 $)}"),
       y = "Multiplicity of infection (PFU/cell)") +
  coord_cartesian(xlim = c(0,1))

ggsave(file = "T3D_Association_Ridge_Plot.tiff",
       width = 6,
       height = 4.5,
       dpi = 1200)

Base.Plot +
  geom_density_ridges(data = Reo.Corr,
                      aes(y = factor(MOI),
                          x = r2),
                      alpha = 0.8,
                      fill = "navy") +
  facet_wrap(~Rep) +
  labs(x = TeX("\\textbf{Pairwise association ($\\r^2 $)}"),
       y = "Multiplicity of infection (PFU/cell)") +
  coord_cartesian(xlim = c(0,1)) +
  theme(panel.spacing = unit(0.2, "in"))

ggsave(file = "T3D_Association_Ridge_Plot_Facet.tiff",
       width = 10,
       height = 4.5,
       dpi = 1200)

Base.Plot +
  geom_density_ridges(data = Reo.Corr,
                      aes(y = factor(MOI),
                          x = r2,
                          fill = Rep),
                      alpha = 0.6) +
  scale_fill_tableau() +
  labs(x = TeX("\\textbf{Pairwise association ($\\r^2 $)}"),
       y = "Multiplicity of infection (PFU/cell)") +
  coord_cartesian(xlim = c(0,1)) +
  theme(panel.spacing = unit(0.2, "in"))

ggsave(file = "T3D_Association_Ridge_Plot_Individual.tiff",
       width = 6,
       height = 4.5,
       dpi = 1200)

# Scatter Plot ----

Base.Plot +
  geom_jitter(data = Pooled.df %>% filter(X != Y),
              aes(x = MOI %>% log10,
                  y = r2)) +
  scale_color_tableau() +
  labs(x = TeX("\\textbf{MOI ($\\log_1_0$ PFU/cell)}"),
       y = TeX("\\textbf{Pairwise association ($\\r^2 $)}"))

ggsave(file = "T3D_Association_Scatter_Plot_Pooled.tiff",
       width = 6,
       height = 4.5,
       dpi = 1200)

Base.Plot +
  geom_jitter(data = Reo.Corr %>% filter(X != Y) %>% mutate(Pair = str_c(X,Y)),
              aes(x = MOI %>% log10,
                  y = r2)) +
  geom_smooth(data = Reo.Corr %>% 
                filter(X != Y) %>% 
                mutate(Pair = str_c(X,Y)) %>%
                filter(Pair %in% c("L1L2","L1M2")),
              aes(x = MOI %>% log10,
                  y = r2,
                  color = Pair),
              fill = NA) +
  geom_point(data = Reo.Corr %>% 
               filter(X != Y) %>% 
               mutate(Pair = str_c(X,Y)) %>%
               filter(Pair %in% c("L1L2","L1M2")),
             aes(x = MOI %>% log10,
                 y = r2,
                 color = Pair),
             size = 2) +
  scale_color_tableau() +
  labs(x = TeX("\\textbf{MOI ($\\log_1_0$ PFU/cell)}"),
       y = TeX("\\textbf{Pairwise association ($\\r^2 $)}"))

ggsave(file = "T3D_Association_Scatter_Plot_Individual.tiff",
       width = 6,
       height = 4.5,
       dpi = 1200)

Sig.Test = Pooled.df %>%
  filter(X != Y) %>%
  mutate(Pair = str_c(X,Y)) %>%
  group_by(X, Y) %>%
  summarize(Cumulative = sum(r2)) %>% #Sum of r2. Higher sum -> association decays more slowly as MOI increases. Could also use area under curve.
  ungroup %>%
  mutate(Z = (Cumulative - mean(Cumulative)) / sd(Cumulative),
         p = pmin(2*pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE),rep(1,90)), #2*p for two-tailed test
         index = 1:90) 
# Note, these p-values are uncorrected. Of 45 comparisons, 3 with p < 0.05 is not necessarily unexpected

Sig.Test %>% filter (p< 0.05)
Sig.Test

Stair.Plot +
  geom_tile(data = Sig.Test %>% mutate(X = as.character(X),
                                       Y = as.character(Y)) %>%
              filter(X > Y),
            aes(x = X,
                y = Y,
                fill = abs(Z))) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL,
       y = NULL,
       fill = TeX("Z")) +
  theme(legend.position = "NA")

ggsave(file = "T3D_Association_Stair_Plot_Cumulative.tiff",
       width = 5,
       height = 4,
       dpi = 1200)

Legend.Plot = Stair.Plot +
  geom_tile(data = Sig.Test %>% mutate(X = as.character(X),
                                       Y = as.character(Y)) %>%
              filter(X > Y),
            aes(x = X,
                y = Y,
                fill = abs(Z))) +
  coord_equal() +
  scale_y_discrete(position = "right") +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL,
       y = NULL,
       fill = TeX("Z")) +
  theme(legend.position = "bottom")

ggsave(g_legend(Legend.Plot),
       file = "T3D_Association_Stair_Plot_Z_Legend.tiff",
       width = 5,
       height = 1,
       dpi = 1200)
