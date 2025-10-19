
# Packages
library(vegan)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(xtable)

# Data PCAs

carib.pca.cwm.values <- readRDS("Completeness_data_Islands/carib.pca.cwm.values.rds")
carib.pca.cwm.loadings <- readRDS("Completeness_data_Islands/carib.pca.cwm.loadings.rds")
orinma.pca.cwm.values <- readRDS("Completeness_data_Islands/orinma.pca.cwm.values.rds")
orinma.pca.cwm.loadings <- readRDS("Completeness_data_Islands/orinma.pca.cwm.loadings.rds")

summary(carib.pca.cwm.values) # 4 localities without predators, remove
carib.pca.cwm.values <- subset(carib.pca.cwm.values, !is.na(pred_richness))

summary(orinma.pca.cwm.values) # 4 localities without predators, remove
orinma.pca.cwm.values <- subset(orinma.pca.cwm.values, !is.na(pred_richness))

# Analyses

# Quantify dispersion ####
orinma.bd <- betadisper(dist(cbind(orinma.pca.cwm.values$PC1, orinma.pca.cwm.values$PC2)), 
                       orinma.pca.cwm.values$subregion)
TukeyHSD(orinma.bd)

orinma.disp_df <- data.frame(distances = orinma.bd$distances,
                      spp_richness = orinma.pca.cwm.values$spp_richness,
                      pred_richness = orinma.pca.cwm.values$pred_richness,
                      subregion = orinma.pca.cwm.values$subregion)

lm_disp_orinma <- lm(distances ~ spp_richness + pred_richness + subregion, data = orinma.disp_df)
summary(lm_disp_orinma)

# Figures Euclidean distance (Variance - dispersion) ####
ggplot(orinma.disp_df, aes(x = factor(subregion, 
                                      levels = c("Mainland",
                                                 "Continental islands",
                                                 "Andaman & Nicobar",
                                                 "Sunda islands",
                                                 "Philippines",
                                                 "Wallacea",
                                                 "Papua",
                                                 "Solomons",
                                                 "Vanuatu"),
                                      labels = c("Mainland",
                                                 "Continental islands",
                                                 "Andaman & Nicobar",
                                                 "Sunda islands",
                                                 "Philippines",
                                                 "Wallacea",
                                                 "New Guinea",
                                                 "Solomon islands",
                                                 "Vanuatu")),
                           y = spp_richness,
                           fill = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21)+
  labs(x = "Subregion",
       y = "Species richness",
       title = "Betadispersion analysis - Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Sm_Region_IndoP.jpg",
       width = 8, height = 4, units = "in")


ggplot(orinma.disp_df, aes(x = factor(subregion, 
                                      levels = c("Mainland",
                                                 "Continental islands",
                                                 "Andaman & Nicobar",
                                                 "Sunda islands",
                                                 "Philippines",
                                                 "Wallacea",
                                                 "Papua",
                                                 "Solomons",
                                                 "Vanuatu"),
                                      labels = c("Mainland",
                                                 "Continental islands",
                                                 "Andaman & Nicobar",
                                                 "Sunda islands",
                                                 "Philippines",
                                                 "Wallacea",
                                                 "New Guinea",
                                                 "Solomon islands",
                                                 "Vanuatu")),
                           y = distances,
                           fill = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21)+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray")+
  labs(x = "Subregion",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Region_IndoP.jpg",
       width = 8, height = 4, units = "in")


ggplot(orinma.disp_df, aes(x = spp_richness,
                           y = distances,
                           fill = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu")),
                           color = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  scale_color_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15)+
  geom_point(alpha = 0.25, shape = 21)+
  labs(x = "Bird species richness",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Sm_IndoP.jpg",
       width = 8, height = 4, units = "in")

ggplot(orinma.disp_df, aes(x = pred_richness,
                           y = distances,
                           fill = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu")),
                           color = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "Papua",
                                                    "Solomons",
                                                    "Vanuatu"),
                                         labels = c("Mainland",
                                                    "Continental islands",
                                                    "Andaman & Nicobar",
                                                    "Sunda islands",
                                                    "Philippines",
                                                    "Wallacea",
                                                    "New Guinea",
                                                    "Solomon islands",
                                                    "Vanuatu"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  scale_color_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#43a2ca",
                               "#7bccc4",
                               "#a8ddb5"))+
  geom_point(alpha = 0.25, shape = 21, position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15)+
  labs(x = "Bird predators species richness",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Eastern Indo-Pacific",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Pm_IndoP.jpg",
       width = 8, height = 4, units = "in")

# Caribbean Betadisper ####
carib.bd <- betadisper(dist(cbind(carib.pca.cwm.values$PC1, carib.pca.cwm.values$PC2)), 
                       carib.pca.cwm.values$subregion)
TukeyHSD(carib.bd)

carib.disp_df <- data.frame(distances = carib.bd$distances,
                             spp_richness = carib.pca.cwm.values$spp_richness,
                             pred_richness = carib.pca.cwm.values$pred_richness,
                             subregion = carib.pca.cwm.values$subregion)

lm_disp_carib <- lm(distances ~ spp_richness + pred_richness + subregion, data = carib.disp_df)
summary(lm_disp_carib)

# Figures Euclidean distance (Variance - Dispersion) ####
ggplot(carib.disp_df, aes(x = factor(subregion, 
                                     levels = c("Mainland",
                                                "Continental islands",
                                                "Bahamas (Lucayan)",
                                                "Greater Antilles",
                                                "Lesser Antilles (Kalinago)")),
                          y = spp_richness,
                          fill = factor(subregion, 
                                        levels = c("Mainland",
                                                   "Continental islands",
                                                   "Bahamas (Lucayan)",
                                                   "Greater Antilles",
                                                   "Lesser Antilles (Kalinago)"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21)+
  labs(x = "Subregion",
       y = "Bird species richness",
       title = "Betadispersion analysis - Caribbean",
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

ggplot(carib.disp_df, aes(x = factor(subregion, 
                                      levels = c("Mainland",
                                                 "Continental islands",
                                                 "Bahamas (Lucayan)",
                                                 "Greater Antilles",
                                                 "Lesser Antilles (Kalinago)")),
                           y = distances,
                           fill = factor(subregion, 
                                         levels = c("Mainland",
                                                    "Continental islands",
                                                    "Bahamas (Lucayan)",
                                                    "Greater Antilles",
                                                    "Lesser Antilles (Kalinago)"))))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  geom_boxplot(alpha = 0.5, outliers = FALSE)+
  geom_jitter(alpha = 0.15, width = 0.15, shape = 21)+
  labs(x = "Subregion",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Caribbean",
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

ggplot(carib.disp_df, aes(x = spp_richness,
                          y = distances,
                          fill = factor(subregion, 
                                        levels = c("Mainland",
                                                   "Continental islands",
                                                   "Bahamas (Lucayan)",
                                                   "Greater Antilles",
                                                   "Lesser Antilles (Kalinago)")),
                          color = factor(subregion, 
                                        levels = c("Mainland",
                                                   "Continental islands",
                                                   "Bahamas (Lucayan)",
                                                   "Greater Antilles",
                                                   "Lesser Antilles (Kalinago)"))))+
  scale_color_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15)+
  geom_point(alpha = 0.25, shape = 21)+
  labs(x = "Bird species richness",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Caribbean",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Sm_Carib.jpg",
       width = 8, height = 4, units = "in")

ggplot(carib.disp_df, aes(x = pred_richness,
                          y = distances,
                          fill = factor(subregion, 
                                        levels = c("Mainland",
                                                   "Continental islands",
                                                   "Bahamas (Lucayan)",
                                                   "Greater Antilles",
                                                   "Lesser Antilles (Kalinago)")),
                          color = factor(subregion, 
                                        levels = c("Mainland",
                                                   "Continental islands",
                                                   "Bahamas (Lucayan)",
                                                   "Greater Antilles",
                                                   "Lesser Antilles (Kalinago)"))))+
  scale_color_manual(values = c("#A1A1A1",
                               "#481F70",
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               "#E3E418",
                               "#35B779",
                               "#21908C"))+
  geom_point(alpha = 0.25, shape = 21, position = position_jitter(width = 0.2))+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15)+
  labs(x = "Bird predators species richness",
       y = "Euclidean distance",
       title = "Betadispersion analysis - Caribbean",
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         color=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "figures jpg/Dispersion_Pm_Carib.jpg",
       width = 8, height = 4, units = "in")


# PERMANOVA - shifts in mean CWM values ####
names(orinma.pca.cwm.values)
orinma.PERMA <- pairwise.adonis2(orinma.pca.cwm.values[,8:9] ~ factor(subregion)+pred_richness+spp_richness, 
                                 data = orinma.pca.cwm.values, 
                                 by = "margin", # understanding the independent contribution of each variable (marginal effect)
                                 method = "euclidean", nperm = 100)
orinma.PERMA 

library(broom)
tidy_pairwise_result <- map_dfr(
  orinma.PERMA,
  ~ .x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    as_tibble(),
  .id = "Comparison"
) |>
  dplyr::select(c(Comparison, Term, R2, `F`, `Pr(>F)`)) |>
  filter(Comparison != "parent_call",
         Term != "Residual",
         Term != "Total")

print(xtable(tidy_pairwise_result), include.rownames = FALSE)

saveRDS(orinma.PERMA, "Completeness_data_Islands/orinma.PERMA.rds")

pair.orinma.PERMA <- pairwise.adonis(orinma.pca.cwm.values[,8:9], 
                                     factors = factor(orinma.pca.cwm.values$subregion),
                                     sim.method = "euclidean",
                                     perm = 100)
pair.orinma.PERMA


Carib.PERMA <- pairwise.adonis2(carib.pca.cwm.values[,8:9] ~ factor(subregion)+pred_richness+spp_richness, 
                                data = carib.pca.cwm.values, 
                                by = "margin", # understanding the independent contribution of each variable (marginal effect)
                                method = "euclidean", nperm = 100)
Carib.PERMA # Oceanic archipelagos are way different!
saveRDS(Carib.PERMA, "Completeness_data_Islands/Carib.PERMA.rds")

