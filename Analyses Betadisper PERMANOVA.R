
# Packages
library(vegan)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(xtable)
library(tidyverse)
library(broom)
library(MetBrewer)
library(patchwork) 

# Colors
as.vector(met.brewer("Hiroshige", 6))
as.vector(met.brewer("Java", 6))

# Data PCAs

carib.pca.cwm.values <- readRDS("Completeness_data_Islands/carib.pca.cwm.values.rds")
carib.pca.cwm.loadings <- readRDS("Completeness_data_Islands/carib.pca.cwm.loadings.rds")
orinma.pca.cwm.values <- readRDS("Completeness_data_Islands/orinma.pca.cwm.values.rds")
orinma.pca.cwm.loadings <- readRDS("Completeness_data_Islands/orinma.pca.cwm.loadings.rds")

summary(carib.pca.cwm.values) # 4 localities without predators, remove
carib.pca.cwm.values <- subset(carib.pca.cwm.values, !is.na(pred_richness))

summary(orinma.pca.cwm.values) # 3 localities without predators, remove
orinma.pca.cwm.values <- subset(orinma.pca.cwm.values, !is.na(pred_richness))

# Analyses

# Quantify dispersion - Expansion ####
## Indo-Pacific Betadisper ####
orinma.bd <- betadisper(dist(cbind(orinma.pca.cwm.values[,8:10])), 
                       orinma.pca.cwm.values$subregion, type = "centroid")
TukeyHSD(orinma.bd)

orinma.disp_df <- data.frame(distances = orinma.bd$distances,
                      spp_richness = orinma.pca.cwm.values$spp_richness,
                      pred_richness = orinma.pca.cwm.values$pred_richness,
                      subregion = orinma.pca.cwm.values$subregion)

# compare models
orinma.mod.Null <- glm(distances ~ 1, 
                       data = orinma.disp_df)
orinma.mod.Subr <- glm(distances ~ subregion, 
                      data = orinma.disp_df)
orinma.mod.Sm <- glm(distances ~ spp_richness, 
                    data = orinma.disp_df)
orinma.mod.Pm <- glm(distances ~ pred_richness, 
                    data = orinma.disp_df)
orinma.mod.Sm.Subr <- glm(distances ~ spp_richness*subregion, 
                         data = orinma.disp_df)
orinma.mod.Pm.Subr <- glm(distances ~ pred_richness*subregion, 
                         data = orinma.disp_df)
AICcmodavg::aictab(cand.set = list(orinma.mod.Null, 
                                   orinma.mod.Subr,
                                   orinma.mod.Sm, 
                                   orinma.mod.Pm,
                                   orinma.mod.Sm.Subr,
                                   orinma.mod.Pm.Subr),
                   modnames = c("Null", "Subregion","Sm","Pm","Sm_Subregion","Pm_Subregion"))

summary(orinma.mod.Sm.Subr)#AICcWt 0.96!

#∆AIC = 6.55:
summary(orinma.mod.Subr) 

# Extract effects of each variable
library(emmeans)
orinma.mod.Sm.Subr.eff <- emtrends(orinma.mod.Sm.Subr, ~ subregion, var = "spp_richness") |>
  as.data.frame()

orinma.mod.Sm.Subr.eff |>
  ggplot(aes(x = subregion)) +
  geom_pointrange(aes(y = spp_richness.trend, ymin = lower.CL, ymax = upper.CL)) +
  coord_flip()



## Figures Euclidean distance (Variance - dispersion) ####
ggplot(orinma.disp_df, aes(x = subregion,
                           y = spp_richness,
                           fill = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21, aes(size = spp_richness))+
  labs(x = "Subregion",
       y = "Species richness",
       title = "Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray")+
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "gray")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Sm_Region_IndoP.jpg",
       width = 8, height = 4, units = "in")

summary(orinma.mod.Sm.Subr)
Dispersion_Sm_IndoP <- ggplot(orinma.disp_df, aes(x = spp_richness,
                                                  y = distances,
                                                  fill = subregion,
                                                  color = subregion,
                                                  linetype = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#BA3841",
                                "#DF5C2D",
                                "#E6824F",
                                "#B7917F",
                                "#0C7156"))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  scale_linetype_manual(values = c("solid",
                                   "dashed",
                                   
                                   "dotted",
                                   "dotted",
                                   "solid",
                                   "dotted",
                                   "dotted"))+
  geom_point(alpha = 0.1, shape = 21)+
  labs(x = "Bird species richness",
       y = "Euclidean distance from centroid",
       title = "Expansion - Eastern Indo-Pacific",
       fill = "",
       color = "",
       linetype = "")+
  theme(legend.position = "right",
        legend.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE),
         color=guide_legend(ncol=2,byrow=TRUE),
         linetype=guide_legend(ncol=2,byrow=TRUE))
Dispersion_Sm_IndoP
ggsave(filename = "figures jpg/Dispersion_Sm_IndoP.jpg",Dispersion_Sm_IndoP,
       width = 14, height = 4, units = "in", dpi = 300)

ggplot(orinma.disp_df, aes(x = spp_richness,
                           y = distances,
                           fill = subregion,
                           color = subregion,
                           linetype = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#BA3841",
                                "#DF5C2D",
                                "#E6824F",
                                "#B7917F",
                                "#0C7156"))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  scale_linetype_manual(values = c("solid",
                                   "dashed",
                                   
                                   "dotted",
                                   "dotted",
                                   "solid",
                                   "dotted",
                                   "dotted"))+
  geom_point(alpha = 0.1, shape = 21)+
  labs(x = "Bird species richness",
       y = "Euclidean distance from centroid",
       title = "Expansion - Eastern Indo-Pacific",
       fill = "",
       color = "",
       linetype = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75,0.85),
        legend.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE),
         color=guide_legend(ncol=2,byrow=TRUE),
         linetype=guide_legend(ncol=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Sm_IndoP.jpg",
       width = 8, height = 4, units = "in")


summary(orinma.mod.Subr)
Dispersion_Region_IndoP <- ggplot(orinma.disp_df, aes(x = subregion,
                           y = distances,
                           fill = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21, aes(size = spp_richness))+
  geom_vline(xintercept = c(5.5, 6.5), linetype = "dashed", color = "gray")+
  labs(x = "Subregion",
       y = "Euclidean distance from centroid",
       title = "Expansion - Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  coord_cartesian(ylim = c(0,10))+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank())+
  guides(fill=element_blank(), color = element_blank())
Dispersion_Region_IndoP
ggsave(filename = "figures jpg/Dispersion_Region_IndoP.jpg",Dispersion_Region_IndoP,
       width = 8, height = 4, units = "in")


ggplot(orinma.disp_df, aes(x = pred_richness,
                           y = distances,
                           fill = subregion,
                           color = subregion, 
                           linetype = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#BA3841",
                                "#DF5C2D",
                                "#E6824F",
                                "#B7917F",
                                "#0C7156"))+
  geom_point(alpha = 0.15, shape = 21, position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  scale_linetype_manual(values = c("dotted",
                                   "dotted",
                                   
                                   "dotted",
                                   "dotted",
                                   "dotted",
                                   "dashed",
                                   "dotted"))+
  labs(x = "Bird predators species richness",
       y = "Euclidean distance from centroid",
       title = "Expansion - Eastern Indo-Pacific",
       fill = "",
       color = "", 
       linetype = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75,0.85),
        legend.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE),
         color=guide_legend(ncol=2,byrow=TRUE),
         linetype=guide_legend(ncol=2,byrow=TRUE))

ggsave(filename = "figures jpg/Dispersion_Pm_IndoP.jpg",
       width = 8, height = 4, units = "in")

## Caribbean Betadisper ####
carib.bd <- betadisper(dist(cbind(carib.pca.cwm.values[,8:10])), 
                       carib.pca.cwm.values$subregion, type = "centroid")
TukeyHSD(carib.bd) # no difference in distance

carib.disp_df <- data.frame(distances = carib.bd$distances,
                             spp_richness = carib.pca.cwm.values$spp_richness,
                             pred_richness = carib.pca.cwm.values$pred_richness,
                             subregion = carib.pca.cwm.values$subregion)

# compare models
carib.mod.Null <- glm(distances ~ 1, 
                      data = carib.disp_df)
carib.mod.Subr <- glm(distances ~ subregion, 
                    data = carib.disp_df)
carib.mod.Sm <- glm(distances ~ spp_richness, 
                    data = carib.disp_df)
carib.mod.Pm <- glm(distances ~ pred_richness, 
                    data = carib.disp_df)
carib.mod.Sm.Subr <- glm(distances ~ spp_richness*subregion + pred_richness, 
                    data = carib.disp_df)
carib.mod.Pm.Subr <- glm(distances ~ pred_richness*subregion + spp_richness, 
                    data = carib.disp_df)
AICcmodavg::aictab(cand.set = list(carib.mod.Null, 
                                   carib.mod.Subr,
                                   carib.mod.Sm, 
                                   carib.mod.Pm,
                                   carib.mod.Sm.Subr,
                                   carib.mod.Pm.Subr),
                   modnames = c("Null", "Subregion","Sm","Pm","Sm_Subregion","Pm_Subregion"))

summary(carib.mod.Pm.Subr)
summary(carib.mod.Sm.Subr)

carib.mod.Pm.Subr.eff <- emtrends(carib.mod.Pm.Subr, ~ subregion, var = "pred_richness") |>
  as.data.frame()

carib.mod.Pm.Subr.eff |>
  ggplot(aes(x = subregion)) +
  geom_pointrange(aes(y = pred_richness.trend, ymin = lower.CL, ymax = upper.CL)) +
  coord_flip()


carib.mod.Sm.Subr.eff <- emtrends(carib.mod.Sm.Subr, ~ subregion, var = "spp_richness") |>
  as.data.frame()

carib.mod.Sm.Subr.eff |>
  ggplot(aes(x = subregion)) +
  geom_pointrange(aes(y = spp_richness.trend, ymin = lower.CL, ymax = upper.CL)) +
  coord_flip()



## Figures Euclidean distance (Variance - Dispersion) ####
ggplot(carib.disp_df, aes(x = subregion,
                          y = spp_richness,
                          fill = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21, aes(size = spp_richness))+
  labs(x = "Subregion",
       y = "Bird species richness",
       title = "Caribbean",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Sm_Region_Carib.jpg",
       width = 8, height = 4, units = "in")

summary(carib.mod.Pm.Subr)
Dispersion_Pm_Carib <- ggplot(carib.disp_df, aes(x = pred_richness,
                          y = distances,
                          fill = subregion,
                          color = subregion, 
                          linetype = subregion))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#aadce0",
                                "#376795",
                                "#ef8a47"))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  geom_point(alpha = 0.1, shape = 21, position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  scale_linetype_manual(values = c("solid",
                                   "solid",
                                   
                                   "solid",
                                   "dotted",
                                   "dotted"))+
  labs(x = "Bird predators species richness",
       y = "Euclidean distance from centroid",
       title = "Expansion - Caribbean",
       fill = "",
       color = "",
       linetype = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75,0.9),
        legend.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(ncol =2,byrow=TRUE),
         color=guide_legend(ncol=2,byrow=TRUE),
         linetype=guide_legend(ncol=2,byrow=TRUE))
Dispersion_Pm_Carib
ggsave(filename = "figures jpg/Dispersion_Pm_Carib.jpg",Dispersion_Pm_Carib,
       width = 7, height = 4, units = "in", dpi = 300)

summary(carib.mod.Sm.Subr)
Dispersion_Sm_Carib <- ggplot(carib.disp_df, aes(x = spp_richness,
                                                 y = distances,
                                                 fill = subregion,
                                                 color = subregion,
                                                 linetype = subregion))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#aadce0",
                                "#376795",
                                "#ef8a47"))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  geom_point(alpha = 0.1, shape = 21)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  scale_linetype_manual(values = c("solid",
                                   "dashed",
                                   
                                   "solid",
                                   "solid",
                                   "solid"))+
  labs(x = "Bird species richness",
       y = "Euclidean distance from centroid",
       title = "Expansion - Caribbean",
 #      tag = expression(bold("A")),
       fill = "",
       color = "",
       linetype = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75,0.9),
        legend.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE),
         color=guide_legend(nrow=3,byrow=TRUE),
         linetype=guide_legend(nrow=3,byrow=TRUE))
Dispersion_Sm_Carib
ggsave(filename = "figures jpg/Dispersion_Sm_Carib.jpg",Dispersion_Sm_Carib,
       width = 7, height = 4, units = "in", dpi = 300)

ggplot(carib.disp_df, aes(x = subregion,
                           y = distances,
                           fill = subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.1, width = 0.15, shape = 21, aes(size = spp_richness))+
  labs(x = "Subregion",
       y = "Euclidean distance from centroid",
       title = "Expansion - Caribbean",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Region_Carib.jpg",
       width = 8, height = 4, units = "in")

# 



# PERMANOVA - shifts in mean CWM values ####


## Pairwise adonis - PERMANOVA IndoPacific ####
# Combine the two factors into a single interaction variable

orinma.PERMA <- pairwise.adonis2(orinma.pca.cwm.values[,8:10] ~ subregion + spp_richness + pred_richness, 
                                 data = orinma.pca.cwm.values, 
                                 by = "margin", # understanding the independent contribution of each variable (marginal effect)
                                 method = "euclidean", nperm = 999)
orinma.PERMA 

# using broom
orinma_pairwise_result <- map_dfr(
  orinma.PERMA,
  ~ .x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    as_tibble(),
  .id = "Comparison"
) |>
  as.data.frame() |>
  dplyr::select(c(Comparison, Term, R2, `F`, `Pr(>F)`)) |>
  filter(Comparison != "parent_call",
         Term != "Residual",
         Term != "Total") |>
  mutate(R2 = round(R2, 3),
         F = round(F, 3))

print(xtable(orinma_pairwise_result, digits = 3), include.rownames = FALSE)

saveRDS(orinma.PERMA, "Completeness_data_Islands/orinma.PERMA.rds")

orinma.PERMA <- readRDS("Completeness_data_Islands/orinma.PERMA.rds")

# compare models
orinma.PC1.mod.Subr <- glm(PC1 ~ subregion, 
                          data = orinma.pca.cwm.values)
orinma.PC1.mod.Sm <- glm(PC1 ~ spp_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC1.mod.Pm <- glm(PC1 ~ pred_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC1.mod.Sm.Subr <- glm(PC1 ~ spp_richness*subregion, 
                             data = orinma.pca.cwm.values)
orinma.PC1.mod.Pm.Subr <- glm(PC1 ~ pred_richness*subregion, 
                             data = orinma.pca.cwm.values)

orinma.PC2.mod.Subr <- glm(PC2 ~ subregion, 
                          data = orinma.pca.cwm.values)
orinma.PC2.mod.Sm <- glm(PC2 ~ spp_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC2.mod.Pm <- glm(PC2 ~ pred_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC2.mod.Sm.Subr <- glm(PC2 ~ spp_richness*subregion, 
                             data = orinma.pca.cwm.values)
orinma.PC2.mod.Pm.Subr <- glm(PC2 ~ pred_richness*subregion, 
                             data = orinma.pca.cwm.values)

orinma.PC3.mod.Subr <- glm(PC3 ~ subregion, 
                          data = orinma.pca.cwm.values)
orinma.PC3.mod.Sm <- glm(PC3 ~ spp_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC3.mod.Pm <- glm(PC3 ~ pred_richness, 
                        data = orinma.pca.cwm.values)
orinma.PC3.mod.Sm.Subr <- glm(PC3 ~ spp_richness*subregion, 
                             data = orinma.pca.cwm.values)
orinma.PC3.mod.Pm.Subr <- glm(PC3 ~ pred_richness*subregion, 
                             data = orinma.pca.cwm.values)

AICcmodavg::aictab(cand.set = list(orinma.PC1.mod.Subr,
                                   orinma.PC1.mod.Sm, 
                                   orinma.PC1.mod.Pm,
                                   orinma.PC1.mod.Sm.Subr,
                                   orinma.PC1.mod.Pm.Subr),
                   modnames = c("Subregion.PC1","Sm.PC1","Pm.PC1","Sm_Subregion.PC1","Pm_Subregion.PC1"))

AICcmodavg::aictab(cand.set = list(orinma.PC2.mod.Subr,
                                   orinma.PC2.mod.Sm, 
                                   orinma.PC2.mod.Pm,
                                   orinma.PC2.mod.Sm.Subr,
                                   orinma.PC2.mod.Pm.Subr),
                   modnames = c("Subregion.PC2","Sm.PC2","Pm.PC2","Sm_Subregion.PC2","Pm_Subregion.PC2"))

AICcmodavg::aictab(cand.set = list(orinma.PC3.mod.Subr,
                                   orinma.PC3.mod.Sm, 
                                   orinma.PC3.mod.Pm,
                                   orinma.PC3.mod.Sm.Subr,
                                   orinma.PC3.mod.Pm.Subr),
                   modnames = c("Subregion.PC3","Sm.PC3","Pm.PC3","Sm_Subregion.PC3","Pm_Subregion.PC3"))

summary(orinma.PC1.mod.Sm.Subr)
summary(orinma.PC2.mod.Sm.Subr)
summary(orinma.PC3.mod.Sm.Subr)

## Figures mean per region PC1-PC3 IndoPacific ####
names(orinma.pca.cwm.values)

# Define new labels
PC_labels_IndoP <- c("PC1" = " PC1\n(59.37%)",
                     "PC2" = "PC2\n(16.66%)",
                     "PC3" = "PC3\n(8.48%)")


Shift_Region_IndoP <- orinma.pca.cwm.values |>
  dplyr::select(cell, Meta_Archipelago, region, subregion, fig_group, 
                spp_richness, pred_richness, PC1, PC2, PC3) |>
  pivot_longer(cols = c(PC1, PC2, PC3)) |>
  mutate(linetype = case_when(
    # PC1 rules
    name == "PC1" & subregion %in% c("Mainland", "Andaman & Nicobar", "Sunda islands", "Philippines") ~ "solid",
    name == "PC1" & subregion %in% c("Continental islands", "New Guinea") ~ "dotted",
    name == "PC1" & subregion == "Wallacea" ~ "dashed",
    
    # PC2 rules
    name == "PC2" & subregion %in% c("Mainland", "Sunda islands", "New Guinea") ~ "solid",
    name == "PC2" & subregion %in% c("Continental islands", "Andaman & Nicobar", "Philippines", "Wallacea") ~ "dotted",
    
    # PC3 rules
    name == "PC3" & subregion %in% c("Mainland", "Sunda islands") ~ "solid",
    name == "PC3" & subregion %in% c("Continental islands", "Andaman & Nicobar", "Philippines", "Wallacea", "New Guinea") ~ "dotted",
    
    TRUE ~ "solid"  # fallback
  )) |> 
  ggplot(aes(x = spp_richness,
             y = value,
             fill = subregion,
             color = subregion)) +
  facet_wrap(~name, scales = "free_y", ncol = 1, 
             labeller = as_labeller(PC_labels_IndoP)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.1) +
  geom_smooth(aes(linetype = linetype), method = "lm", se = TRUE, alpha = 0.1) +
  scale_linetype_identity() +
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#BA3841",
                               "#DF5C2D",
                               "#E6824F",
                               "#B7917F",
                               "#0C7156"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#BA3841",
                                "#DF5C2D",
                                "#E6824F",
                                "#B7917F",
                                "#0C7156"))+
  labs(x = "Species richness",
       y = "Score", 
       title = "Shift of PC score with species richness",
       subtitle = "Eastern Indo-Pacific")+
  theme_classic() +
  theme(legend.position = "left",
        legend.title = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"),
        strip.background = element_blank())

Shift_Region_IndoP
ggsave(filename = "figures jpg/Shift_Region_IndoP.jpg",Shift_Region_IndoP,
       width = 6, height = 7, units = "in", dpi = 300)
## Pairwise adonis - PERMANOVA Caribbean ####
Carib.PERMA <- readRDS("Completeness_data_Islands/Carib.PERMA.rds")

Carib.PERMA <- pairwise.adonis2(carib.pca.cwm.values[,8:10] ~ subregion + spp_richness + pred_richness, 
                                data = carib.pca.cwm.values, 
                                by = "margin", # understanding the independent contribution of each variable (marginal effect)
                                method = "euclidean", nperm = 999)
Carib.PERMA # Oceanic archipelagos are way different!

Carib_pairwise_result <- map_dfr(
  Carib.PERMA,
  ~ .x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    as_tibble(),
  .id = "Comparison"
) |>
  as.data.frame() |>
  dplyr::select(c(Comparison, Term, R2, `F`, `Pr(>F)`)) |>
  filter(Comparison != "parent_call",
         Term != "Residual",
         Term != "Total") 

print(xtable(Carib_pairwise_result, digits = 3), include.rownames = FALSE)

saveRDS(Carib.PERMA, "Completeness_data_Islands/Carib.PERMA.rds")

# compare models
Carib.PC1.mod.Null <- glm(PC1 ~ 1, 
                          data = carib.pca.cwm.values)
Carib.PC1.mod.Subr <- glm(PC1 ~ subregion, 
                       data = carib.pca.cwm.values)
Carib.PC1.mod.Sm <- glm(PC1 ~ spp_richness, 
                     data = carib.pca.cwm.values)
Carib.PC1.mod.Pm <- glm(PC1 ~ pred_richness, 
                     data = carib.pca.cwm.values)
Carib.PC1.mod.Sm.Subr <- glm(PC1 ~ spp_richness*subregion, 
                          data = carib.pca.cwm.values)
Carib.PC1.mod.Pm.Subr <- glm(PC1 ~ pred_richness*subregion, 
                          data = carib.pca.cwm.values)

Carib.PC2.mod.Null <- glm(PC2 ~ 1, 
                          data = carib.pca.cwm.values)
Carib.PC2.mod.Subr <- glm(PC2 ~ subregion, 
                          data = carib.pca.cwm.values)
Carib.PC2.mod.Sm <- glm(PC2 ~ spp_richness, 
                        data = carib.pca.cwm.values)
Carib.PC2.mod.Pm <- glm(PC2 ~ pred_richness, 
                        data = carib.pca.cwm.values)
Carib.PC2.mod.Sm.Subr <- glm(PC2 ~ spp_richness*subregion, 
                             data = carib.pca.cwm.values)
Carib.PC2.mod.Pm.Subr <- glm(PC2 ~ pred_richness*subregion, 
                             data = carib.pca.cwm.values)

Carib.PC3.mod.Null <- glm(PC3 ~ 1, 
                          data = carib.pca.cwm.values)
Carib.PC3.mod.Subr <- glm(PC3 ~ subregion, 
                          data = carib.pca.cwm.values)
Carib.PC3.mod.Sm <- glm(PC3 ~ spp_richness, 
                        data = carib.pca.cwm.values)
Carib.PC3.mod.Pm <- glm(PC3 ~ pred_richness, 
                        data = carib.pca.cwm.values)
Carib.PC3.mod.Sm.Subr <- glm(PC3 ~ spp_richness*subregion, 
                             data = carib.pca.cwm.values)
Carib.PC3.mod.Pm.Subr <- glm(PC3 ~ pred_richness*subregion, 
                             data = carib.pca.cwm.values)

AICcmodavg::aictab(cand.set = list(Carib.PC1.mod.Null,
                                   Carib.PC1.mod.Subr,
                                   Carib.PC1.mod.Sm, 
                                   Carib.PC1.mod.Pm,
                                   Carib.PC1.mod.Sm.Subr,
                                   Carib.PC1.mod.Pm.Subr),
                   modnames = c("Null", "Subregion.PC1","Sm.PC1","Pm.PC1","Sm_Subregion.PC1","Pm_Subregion.PC1"))

AICcmodavg::aictab(cand.set = list(Carib.PC2.mod.Null, 
                                   Carib.PC2.mod.Subr,
                                   Carib.PC2.mod.Sm, 
                                   Carib.PC2.mod.Pm,
                                   Carib.PC2.mod.Sm.Subr,
                                   Carib.PC2.mod.Pm.Subr),
                   modnames = c("Null", "Subregion.PC2","Sm.PC2","Pm.PC2","Sm_Subregion.PC2","Pm_Subregion.PC2"))

AICcmodavg::aictab(cand.set = list(Carib.PC3.mod.Null, 
                                   Carib.PC3.mod.Subr,
                                   Carib.PC3.mod.Sm, 
                                   Carib.PC3.mod.Pm,
                                   Carib.PC3.mod.Sm.Subr,
                                   Carib.PC3.mod.Pm.Subr),
                   modnames = c("Null", "Subregion.PC3","Sm.PC3","Pm.PC3","Sm_Subregion.PC3","Pm_Subregion.PC3"))

summary(Carib.PC1.mod.Sm.Subr)
summary(Carib.PC2.mod.Sm.Subr)
summary(Carib.PC3.mod.Sm.Subr)

## Figures mean per region PC1-PC3 Caribbean ####
names(carib.pca.cwm.values)

# Define new labels
PC_labels_Carib <- c("PC1" = " PC1\n(55.63%)",
                     "PC2" = "PC2\n(14.42%)",
                     "PC3" = "PC3\n(12.92%)")


Shift_Region_Carib <- carib.pca.cwm.values |>
  dplyr::select(cell, Meta_Archipelago, region, subregion, fig_group, 
                spp_richness, pred_richness, PC1, PC2, PC3) |>
  pivot_longer(cols = c(PC1, PC2, PC3)) |>
  mutate(linetype = case_when(
    # PC1 rules
    name == "PC1" & subregion %in% c("Mainland", "Lesser Antilles (Kalinago)", "Greater Antilles", "Bahamas (Lucayan)") ~ "solid",
    name == "PC1" & subregion == "Continental islands" ~ "dotted",
    
    # PC2 rules
    name == "PC2" & subregion %in% c("Mainland", "Lesser Antilles (Kalinago)", "Greater Antilles") ~ "solid",
    name == "PC2" & subregion %in% c("Continental islands", "Bahamas (Lucayan)") ~ "dotted",
    
    # PC3 rules
    name == "PC3" & subregion == "Mainland" ~ "solid",
    name == "PC3" & subregion %in% c("Continental islands", "Bahamas (Lucayan)", "Lesser Antilles (Kalinago)", "Greater Antilles") ~ "dotted",
    
    TRUE ~ "solid"  # fallback
  )) |>
  ggplot(aes(x = spp_richness,
             y = value,
             fill = subregion,
             color = subregion)) +
  facet_wrap(~name, scales = "free_y", ncol = 1, 
             labeller = as_labeller(PC_labels_Carib)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.1) +
  geom_smooth(aes(linetype = linetype), method = "lm", se = TRUE, alpha = 0.1) +
  scale_linetype_identity() +
  scale_fill_manual(values = c("#A1A1A1", "#663171", "#aadce0", "#376795", "#ef8a47")) +
  scale_color_manual(values = c("#A1A1A1", "#663171", "#aadce0", "#376795", "#ef8a47")) +
  labs(x = "Species richness",
       y = "Score",
       title = "Shift – Caribbean") +
  theme_classic() +
  theme(legend.position = "left",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"),
        strip.background = element_blank())

Shift_Region_Carib
ggsave(filename = "figures jpg/Shift_Region_Carib.jpg",Shift_Region_Carib,
       width = 6, height = 7, units = "in", dpi = 300)

