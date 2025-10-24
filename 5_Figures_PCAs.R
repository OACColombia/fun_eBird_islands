# Figure of the PCAs Functional Biogeography in three meta-archipelagos

library(tidyverse)
library(ggExtra) # for ggMarginal() of the density plot
library(gridExtra) #for grid.arrange()
library(grid)
#library(ggpubr) #ggarrange() with common legend
library(ggrepel) # for geom_text_repel()
library(funspace) # functional trait space
library(ggstar) # points as hexagons
library(GGally)


# Call the data ####

## Global and Island functional traits ####
functraits_global <- readRDS("FunctionalTraits/birds_functraits_imputed.rds")
functraits_islands <- readRDS("FunctionalTraits/FuncTraits_birds_Islands.rds")

## Species Richness per cell and Predator Species Richness per cell ####
Caribbean_Cell_list <- readRDS("Completeness_data_Islands/Caribbean_Cell_list.rds")
IndoPacific_Cell_list <- readRDS("Completeness_data_Islands/IndoPacific_Cell_list.rds")

C_SR <- Caribbean_Cell_list |>
  ungroup() |>
  left_join(functraits_global, join_by(scientific_name == Species2_eBird)) |>
  group_by(cell) |>
  summarise(SppRichnessCell = n()) |>
  as.data.frame()

C_PSR <- Caribbean_Cell_list |>
  ungroup() |>
  left_join(functraits_global, join_by(scientific_name == Species2_eBird)) |>
  filter(predatory == 1) |>
  group_by(cell) |>
  summarise(PredSppRichnessCell = n()) |>
  as.data.frame()

Caribbean_SppR_Pred_SppR <- C_SR |>
  left_join(C_PSR)

IP_SR <- IndoPacific_Cell_list |>
  ungroup() |>
  left_join(functraits_global, join_by(scientific_name == Species2_eBird)) |>
  group_by(cell) |>
  summarise(SppRichnessCell = n()) |>
  as.data.frame()

IP_PSR <- IndoPacific_Cell_list |>
  ungroup() |>
  left_join(functraits_global, join_by(scientific_name == Species2_eBird)) |>
  filter(predatory == 1) |>
  group_by(cell) |>
  summarise(PredSppRichnessCell = n()) |>
  as.data.frame()

IndoPacific_SppR_Pred_SppR <- IP_SR |>
  left_join(IP_PSR) 

islands_SppR_Pred_SppR <- bind_rows(Caribbean_SppR_Pred_SppR, IndoPacific_SppR_Pred_SppR) |>
  mutate(cell = as.character(cell))

## Community Weighted Means and Community Weighted Variation of traits per cell ####
islands_cwm_mu <- readRDS("Completeness_data_Islands/Islands_CWM_mu.rds") |> # summary()
  filter(spp_richness >= 17) |> # first quartile of # of species with OUSS fitted model
  mutate(`verticality cwv` = `verticality cwv`+0.01)

islands_cwm_mu_sem <- islands_cwm_mu |>
  mutate(Meta_Archipelago = ifelse(Meta.Archipelago %in% "Indo.Pacific",
                                   "Indo-Pacific",
                                          "Caribbean")) |> 
  left_join(islands_SppR_Pred_SppR)

spp_cell_df <- readRDS("Completeness_data_Islands/Species_in_cells_MetaArchipelagos.rds")

# PCA Caribbean CWM ####

caribbean_cwm <- islands_cwm_mu_sem |>
  filter(Meta_Archipelago %in% "Caribbean") |>
  # bring back the mean value of longitude for each cell
  left_join(spp_cell_df |> 
              dplyr::select(cell,longitude, latitude) |> 
              group_by(cell) |> 
              summarize(longitude = mean(longitude),
                        latitude = mean(latitude))) |>
  filter(latitude < 24 | longitude > -79.5)

names(caribbean_cwm)

caribbean_cwm[,c(3,5:14,17:20)] |>
  pivot_longer(cols = !c(subregion), names_to = "Trait",values_to = "CWM") |>
  ggplot(aes(x = scale(log10(CWM)))) + geom_histogram() +facet_wrap(~Trait, scales = "free")+theme_classic()


funspaceDim(scale(log10(caribbean_cwm[,c(5:14,17:20)]))) # 3 dimensions

carib.pca.cwm <- prcomp(scale(log10(caribbean_cwm[,c(5:14,17:20)])))
carib.pca.cwm.values <- data.frame(cell = caribbean_cwm$cell,
                               Meta_Archipelago = caribbean_cwm$Meta_Archipelago,
                               region = caribbean_cwm$region,
                               subregion = caribbean_cwm$subregion,
                               fig_group = caribbean_cwm$fig_group,
                               spp_richness = caribbean_cwm$SppRichnessCell,
                               pred_richness = caribbean_cwm$PredSppRichnessCell,
                               carib.pca.cwm$x)

carib.pca.cwm.loadings <- data.frame(Variables = rownames(carib.pca.cwm$rotation), 
                                     carib.pca.cwm$rotation)

aload.car <- abs(carib.pca.cwm$rotation)
sweep(aload.car, 2, colSums(aload.car), "/")


table(carib.pca.cwm.values$subregion)

carib.pca.cwm.values$subregion <- factor(carib.pca.cwm.values$subregion,
                                          levels = c("Mainland",
                                                     "Continental islands",
                                                     "Lesser Antilles (Kalinago)",
                                                     "Greater Antilles",
                                                     "Bahamas (Lucayan)"))

caribe <- ggpairs(carib.pca.cwm.values,
        columns = 8:10,
        aes(color = subregion,
            fill = subregion,
            alpha = 0.001),
        upper = list(continuous = wrap("density", alpha = 0.1, bins = 8)))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#aadce0",
                                "#376795",
                                "#ef8a47"))+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        strip.background = element_blank())

saveRDS(caribe, "figures jpg/Caribe_PC1PC3.RDS")

## Figure CWM-a - PC1-PC2 Caribbean ####
summary(carib.pca.cwm.values[,8:10])
summary(carib.pca.cwm)

cwm.a <- ggplot(carib.pca.cwm.values) +
    geom_density2d(aes(x = PC1, y = PC2,
                       color = subregion), 
                          alpha = 0.5, bins = 10)+
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC1, 
                y = PC2, 
                fill = subregion,
                color = subregion, 
                size = spp_richness), 
            alpha = 0.1,
            starshape = 6,
            starstroke = 0.5) +
#  facet_wrap(vars(subregion))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#aadce0",
                                "#376795",
                                "#ef8a47"))+
  scale_x_continuous(limits = c(-9.1,9.1))+
  scale_y_continuous(limits = c(-6.1,6.1))+
  coord_fixed()+
#  geom_segment(data = carib.pca.cwm.loadings |>
#                 filter(Variables %in% c("hand wing index",
#                                         "beak width",
#                                         "beak depth",
#                                         "verticality")), 
#               linewidth = 0.25,
#               aes(x = 0, xend = PC1*10, 
#                   y = 0, yend = PC2*10),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwm.loadings, 
                  aes(x = PC1*20, y = PC2*12, 
                      label = Variables),
                  size = 3, color = "black") + 
  labs(x = "PC1 (55.63%)",
       y = "PC2 (14.42%)",
#       title = expression(bold("a")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        strip.background = element_blank())+
  guides(fill=guide_legend(ncol=2,byrow=TRUE, override.aes = list(alpha = 0.6)))

cwm.density.a <- ggMarginal(cwm.a,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

ggsave(filename = "figures jpg/CaribbeanPCA1_PCA2_CWM.jpg", plot = cwm.density.a,
       dpi = 300, 
       width = 5, height = 5, units = "in")

## Figure CWM-c - PC2-PC3 Caribbean ####
summary(carib.pca.cwm.values[,8:10])
summary(carib.pca.cwm)

table(carib.pca.cwm.values$subregion)

cwm.c <- ggplot(carib.pca.cwm.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  facet_wrap(vars(subregion))+
  geom_star(aes(x = PC2, 
                y = PC3, 
                fill = subregion, 
                size = spp_richness), 
            alpha = 0.3,
            starshape = 6) +
  geom_density2d(aes(x = PC2, y = PC3,
                     color = subregion), 
                 alpha = 0.5, bins = 5)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("#A1A1A1",
                               "#663171",
                               
                               "#aadce0",
                               "#376795",
                               "#ef8a47"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#663171",
                                
                                "#aadce0",
                                "#376795",
                                "#ef8a47"))+
  scale_x_continuous(limits = c(-6,6))+
  scale_y_continuous(limits = c(-6.1,6.1))+
  coord_fixed()+
#  geom_segment(data = carib.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC2*5, 
#                   y = 0, yend = PC3*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwm.loadings, 
                  aes(x = PC2*6, y = PC3*6, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC2 (14.42%)",
       y = "PC3 (12.92%)",
 #      title = expression(bold("c")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        strip.background = element_blank())+
  guides(fill=guide_legend(ncol = 2,byrow=TRUE, override.aes = list(alpha = 0.6)))

cwm.density.c <- ggMarginal(cwm.c,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

ggsave(filename = "figures jpg/CaribbeanPCA2_PCA3_CWM.jpg", plot = cwm.density.c,
       dpi = 300, 
       width = 5, height = 5, units = "in")

# PCA Indo-Pacific CWM ####
table(islands_cwm_mu_sem$Meta_Archipelago)

orinma_cwm <- islands_cwm_mu_sem |> 
  filter(Meta_Archipelago == "Indo-Pacific") |>
  # bring back the mean value of longitude for each cell
  left_join(spp_cell_df |> 
              dplyr::select(cell,longitude, latitude) |> 
              group_by(cell) |> 
              summarize(longitude = mean(longitude),
                        latitude = mean(latitude))) |>
  filter(subregion != "Vanuatu",
         subregion != "Solomons",
         latitude > -10)

funspaceDim(scale(log10(orinma_cwm[,c(5:14,17:20)]))) # 3 dimensions

orinma.pca.cwm <- prcomp(scale(log10(orinma_cwm[,c(5:14,17:20)])))

orinma.pca.cwm.values <- data.frame(cell = orinma_cwm$cell,
                                Meta_Archipelago = orinma_cwm$Meta_Archipelago,
                                region = orinma_cwm$region,
                                subregion = orinma_cwm$subregion,
                                fig_group = orinma_cwm$fig_group,
                                spp_richness = orinma_cwm$SppRichnessCell,
                                pred_richness = orinma_cwm$PredSppRichnessCell,
                                orinma.pca.cwm$x)

orinma.pca.cwm.loadings <- data.frame(Variables = rownames(orinma.pca.cwm$rotation), orinma.pca.cwm$rotation)

aload <- abs(orinma.pca.cwm$rotation)
sweep(aload, 2, colSums(aload), "/")

table(orinma.pca.cwm.values$subregion)

orinma.pca.cwm.values$subregion <- factor(orinma.pca.cwm.values$subregion,
                                          levels = c("Mainland",
                                                     "Continental islands",
                                                     "Andaman & Nicobar",
                                                     "Sunda islands",
                                                     "Philippines",
                                                     "Wallacea",
                                                     "Papua"),
                                          labels = c("Mainland",
                                                     "Continental islands",
                                                     "Andaman & Nicobar",
                                                     "Sunda islands",
                                                     "Philippines",
                                                     "Wallacea",
                                                     "New Guinea"))

indop <- ggpairs(orinma.pca.cwm.values,
                  columns = 8:10,
                  aes(color = subregion,
                      fill = subregion,
                      alpha = 0.001),
                  upper = list(continuous = wrap("density", alpha = 0.1, bins = 8)))+
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
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"), 
        strip.background = element_blank())

saveRDS(indop, "figures jpg/IndoP_PC1PC3.RDS")

## Figure CWM-b - PC1-PC2 OrInMa ####
summary(orinma.pca.cwm.values[,8:10])
summary(orinma.pca.cwm)



cwm.b <- ggplot(orinma.pca.cwm.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC1, 
                y = PC2, 
                fill = subregion,
                color = subregion), 
            alpha = 0.3,
            starshape = 6) +
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                                                    
                               "#0868ac",
                               "#a8ddb5"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#481F70",
                                
                                "#fc8d62",
                                "#8da0cb",
                                "#66c2a5",
                                
                                "#0868ac",
                                "#a8ddb5"))+
  scale_x_continuous(limits = c(-8.5,9))+
  scale_y_continuous(limits = c(-6,6))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*5, 
#                   y = 0, yend = PC2*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwm.loadings, 
                  aes(x = PC1*8, y = PC2*8, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC1 (59.37%)",
       y = "PC2 (16.66%)",
 #      title = expression(bold("b")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow = 3,byrow=TRUE, override.aes = list(alpha = 0.6)))

cwm.density.b <- ggMarginal(cwm.b,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)
ggsave(filename = "figures jpg/IndoPacificPCA1_PCA2_CWM.jpg", plot = cwm.density.b,
       dpi = 300, 
       width = 5, height = 5, units = "in")

## Figure CWM-d - PC2-PC3 OrInMa ####
summary(orinma.pca.cwm.values[,8:10])
summary(orinma.pca.cwm)

cwm.d <- ggplot(orinma.pca.cwm.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC2, 
                y = PC3, 
                fill = subregion,
                color = subregion), 
            alpha = 0.3,
            starshape = 6) +
  scale_fill_manual(values = c("#A1A1A1",
                               "#481F70",
                               
                               "#fc8d62",
                               "#8da0cb",
                               "#66c2a5",
                               
                               "#0868ac",
                               "#a8ddb5"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#481F70",
                                
                                "#fc8d62",
                                "#8da0cb",
                                "#66c2a5",
                                
                                "#0868ac",
                                "#a8ddb5"))+
  scale_x_continuous(limits = c(-6,6))+
  scale_y_continuous(limits = c(-6.1,6.1))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC2*5, 
#                   y = 0, yend = PC3*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwm.loadings, 
                  aes(x = PC2*8, y = PC3*8, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC2 (16.66%)",
       y = "PC3 (8.48%)",
#       title = expression(bold("d")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow = 3,byrow=TRUE, override.aes = list(alpha = 0.6)))

cwm.density.d <- ggMarginal(cwm.d,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

ggsave(filename = "figures jpg/IndoPacificPCA2_PCA3_CWM.jpg", plot = cwm.density.d,
       dpi = 300, 
       width = 5, height = 5, units = "in")

# Combine figure CWM a-b ####

Figure.CWM <- grid.arrange(cwm.density.a, cwm.density.b,
                          # cwm.density.c, cwm.density.d,
                           ncol = 2)

ggsave(filename = "figures jpg/PCA1_PCA2_CWM.jpg", plot = Figure.CWM,
       dpi = 300, 
       width = 10, height = 5, units = "in")

# Export pcas values and loadings

saveRDS(carib.pca.cwm.values, "Completeness_data_Islands/carib.pca.cwm.values.rds")
saveRDS(carib.pca.cwm.loadings, "Completeness_data_Islands/carib.pca.cwm.loadings.rds")
saveRDS(orinma.pca.cwm.values, "Completeness_data_Islands/orinma.pca.cwm.values.rds")
saveRDS(orinma.pca.cwm.loadings, "Completeness_data_Islands/orinma.pca.cwm.loadings.rds")

# PCA Caribbean CWV ####
names(caribbean_cwm)

caribbean_cwm[,c(3,21:30,33:36)] |>
  pivot_longer(cols = !c(subregion), names_to = "Trait",values_to = "CWV") |>
  ggplot(aes(x = scale(log10(CWV)))) + geom_histogram() +facet_wrap(~Trait, scales = "free")+theme_classic()

funspaceDim(scale(log10(caribbean_cwm[,c(21:30,33:36)]))) # 4 dimensions!

carib.pca.cwv <- prcomp(scale(log10(caribbean_cwm[,c(21:30,33:36)])))
carib.pca.cwv.values <- data.frame(cell = caribbean_cwm$cell,
                                   Meta_Archipelago = caribbean_cwm$Meta_Archipelago,
                                   region = caribbean_cwm$region,
                                   subregion = caribbean_cwm$subregion,
                                   fig_group = caribbean_cwm$fig_group,
                                   spp_richness = caribbean_cwm$SppRichnessCell,
                                   pred_richness = caribbean_cwm$PredSppRichnessCell,
                                   carib.pca.cwv$x)

carib.pca.cwv.loadings <- data.frame(Variables = rownames(carib.pca.cwv$rotation), carib.pca.cwv$rotation)

## Figure CWV-a - PC1-PC2 Caribbean ####
summary(carib.pca.cwv.values[,8:11])
summary(carib.pca.cwv)

cwv.a <- ggplot(carib.pca.cwv.values) +
  geom_point(aes(x = PC1, 
                y = PC2, 
                fill = subregion,
                color = subregion), 
            alpha = 0.0001,
            shape = 21) +
  geom_star(aes(x = PC1, 
                y = PC2, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#E3E418",
                               "#481F70",
                               "#35B779",
                               "#21908C",
                               "#A1A1A1"))+
  scale_color_manual(values = c("#E3E418",
                                "#481F70",
                                "#35B779",
                                "#21908C",
                                "#A1A1A1"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
#  geom_segment(data = carib.pca.cwv.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*10, 
#                   y = 0, yend = PC2*8),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC1*12, y = PC2*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC1 (44.93%)",
       y = "PC2 (13.61%)",
 #             title = expression(bold("a")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.a <- ggMarginal(cwv.a,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure CWV-c - PC2-PC3 Caribbean ####
summary(carib.pca.cwv)

cwv.c <- ggplot(carib.pca.cwv.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC2, 
                y = PC3, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#E3E418",
                               "#481F70",
                               "#35B779",
                               "#21908C",
                               "#A1A1A1"))+
  scale_color_manual(values = c("#E3E418",
                                "#481F70",
                                "#35B779",
                                "#21908C",
                                "#A1A1A1"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
  #  geom_segment(data = carib.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC2*12, y = PC3*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC2 (13.61%)",
       y = "PC3 (11.92%)",
#              title = expression(bold("c")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.c <- ggMarginal(cwv.c,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure CWV-e - PC3-PC4 Caribbean ####
summary(carib.pca.cwv)

cwv.e <- ggplot(carib.pca.cwv.values) +
  geom_point(aes(x = PC3, 
                 y = PC4, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC3, 
                y = PC4, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#E3E418",
                               "#481F70",
                               "#35B779",
                               "#21908C",
                               "#A1A1A1"))+
  scale_color_manual(values = c("#E3E418",
                                "#481F70",
                                "#35B779",
                                "#21908C",
                                "#A1A1A1"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+ #removing cell 750988 Cont Island -12 on PC4
  coord_fixed()+
  #  geom_segment(data = carib.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC3*12, y = PC4*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC3 (11.92%)",
       y = "PC4 (7.94%)",
 #             title = expression(bold("e")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.e <- ggMarginal(cwv.e,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)


# PCA Oriental-Indo-Malayan CWV ####

funspaceDim(scale(log10(orinma_cwm[,c(21:30,33:36)]))) # 4 dimensions

orinma.pca.cwv <- prcomp(scale(log10(orinma_cwm[,c(21:30,33:36)])))

orinma.pca.cwv.values <- data.frame(cell = orinma_cwm$cell,
                                    Meta_Archipelago = orinma_cwm$Meta_Archipelago,
                                    region = orinma_cwm$region,
                                    subregion = orinma_cwm$subregion,
                                    fig_group = orinma_cwm$fig_group,
                                    spp_richness = orinma_cwm$SppRichnessCell,
                                    pred_richness = orinma_cwm$PredSppRichnessCell,
                                    orinma.pca.cwv$x)

orinma.pca.cwv.loadings <- data.frame(Variables = rownames(orinma.pca.cwv$rotation), orinma.pca.cwv$rotation)

## Figure CWV-b - PC1-PC2 OrInMa ####
summary(orinma.pca.cwv.values[,8:11])
summary(orinma.pca.cwv)

cwv.b <- ggplot(orinma.pca.cwv.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC1, 
                y = PC2, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#D44292",
                               "#481F70",
                               "#A1A1A1",
                               "#F6A97A",
                               "#F66D7A"))+
  scale_color_manual(values = c("#D44292",
                                "#481F70",
                                "#A1A1A1",
                                "#F6A97A",
                                "#F66D7A"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwv.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*10, 
#                   y = 0, yend = PC2*8),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC1*12, y = PC2*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC1 (46.06%)",
       y = "PC2 (14.23%)",
#              title = expression(bold("b")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.b <- ggMarginal(cwv.b,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure CWV-d -  PC2-PC3 OrInMa ####
summary(orinma.pca.cwv)

cwv.d <- ggplot(orinma.pca.cwv.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC2, 
                y = PC3, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#D44292",
                               "#481F70",
                               "#A1A1A1",
                               "#F6A97A",
                               "#F66D7A"))+
  scale_color_manual(values = c("#D44292",
                                "#481F70",
                                "#A1A1A1",
                                "#F6A97A",
                                "#F66D7A"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
  #  geom_segment(data = orinma.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC2*12, y = PC3*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC2 (14.23%)",
       y = "PC3 (11.03%)",
  #            title = expression(bold("d")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.d <- ggMarginal(cwv.d,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure CWV-f - PC3-PC4 OrInMa ####
summary(orinma.pca.cwv)

cwv.f <- ggplot(orinma.pca.cwv.values) +
  geom_point(aes(x = PC3, 
                 y = PC4, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.0001,
             shape = 21) +
  geom_star(aes(x = PC3, 
                y = PC4, 
                fill = subregion,
                color = subregion), 
            alpha = 0.4,
            starshape = 6) +
  scale_fill_manual(values = c("#D44292",
                               "#481F70",
                               "#A1A1A1",
                               "#F6A97A",
                               "#F66D7A"))+
  scale_color_manual(values = c("#D44292",
                                "#481F70",
                                "#A1A1A1",
                                "#F6A97A",
                                "#F66D7A"))+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
  #  geom_segment(data = orinma.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC3*12, y = PC4*12, 
                      label = Variables),
                  size = 3) + 
  labs(x = "PC3 (11.03%)",
       y = "PC4 (8.48%)",
 #             title = expression(bold("f")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwv.density.f <- ggMarginal(cwv.f,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Combine figure ####
Figure.CWV <- grid.arrange(cwv.density.a, cwv.density.b, 
                            cwv.density.c, cwv.density.d, 
                            cwv.density.e, cwv.density.f,
                           ncol = 2)

ggsave(filename = "PCAs_CWV.pdf", plot = Figure.CWV,
       dpi = 600, 
       width = 180, height = 300, units = "mm")

# Richness and Predator richness
islands_cwm_mu_sem |> 
  ggplot(aes(x = SppRichnessCell, 
             y = PredSppRichnessCell,
             color = factor(subregion,
                            levels = c("Mainland",
                                       "Continental islands",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Wallacea",
                                       "Papua",
                                       "Solomons",
                                       "Vanuatu", 
                                       "Bahamas (Lucayan)",
                                       "Greater Antilles",
                                       "Lesser Antilles (Kalinago)")),
             fill = factor(subregion,
                            levels = c("Mainland",
                                       "Continental islands",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Wallacea",
                                       "Papua",
                                       "Solomons",
                                       "Vanuatu", 
                                       "Bahamas (Lucayan)",
                                       "Greater Antilles",
                                       "Lesser Antilles (Kalinago)"))))+
  facet_grid(Meta_Archipelago~fig_group) +
  geom_point(aes(y = 0), shape = "|")+
  geom_point(aes(x = 0), shape = "_")+
  geom_point(alpha = 0.25, shap = 21) +
  labs(x = "Bird species richness",
       y = "Bird predator species richness",
       color = "Subregion \n (Archipelago)",
       fill = "Subregion \n (Archipelago)") +
  scale_x_continuous(limits = c(0,550))+
  scale_y_continuous(limits = c(0,35))+
  scale_color_manual(values = c("#A1A1A1",
                                "#481F70",
                                
                                "#fc8d62",
                                "#8da0cb",
                                "#66c2a5",
                                
                                "#0868ac",
                                "#43a2ca",
                                "#7bccc4",
                                "#a8ddb5", 
                                
                                "#E3E418",
                                "#35B779",
                                "#21908C"))+
  scale_fill_manual(values = c("#A1A1A1",
                                "#481F70",
                                
                                "#fc8d62",
                                "#8da0cb",
                                "#66c2a5",
                                
                                "#0868ac",
                                "#43a2ca",
                                "#7bccc4",
                                "#a8ddb5", 
                                
                                "#E3E418",
                                "#35B779",
                                "#21908C"))+
  geom_smooth(method = "lm", se = TRUE) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))
ggsave(filename = "figures jpg/Pm_Sm.jpg",
       width = 8, height = 6, units = "in")


# Extract latitude and longitude for communities used -wrapped files####

wrapped_gridIndoPacific <- readRDS("Completeness_data_Islands/wrapped_gridIndoPacific.rds")
wrapped_gridCaribbean <- readRDS("Completeness_data_Islands/wrapped_gridCaribbean.rds")

wrapped_gridIndoPacific <- wrapped_gridIndoPacific[wrapped_gridIndoPacific$seqnum %in% orinma_cwm$cell, ]
wrapped_gridCaribbean <- wrapped_gridCaribbean[wrapped_gridCaribbean$seqnum %in% caribbean_cwm$cell, ]

wrapped_gridIndoPacific <- wrapped_gridIndoPacific |>
  mutate(cell = as.character(seqnum)) |>
  left_join(orinma_cwm) 

wrapped_gridCaribbean <- wrapped_gridCaribbean |>
  mutate(cell = as.character(seqnum)) |>
  left_join(caribbean_cwm) 

world1 <- sf::st_as_sf(maps::map(database = 'world', plot = FALSE, fill = TRUE))
world1

# OrInMa
mapOrInMa <- ggplot()+
  geom_sf(data=world1)+
  geom_sf(data=wrapped_gridIndoPacific, 
          aes(fill = factor(subregion,
                            levels = c("Mainland",
                                       "Continental islands",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Wallacea",
                                       "Papua")), 
              color = factor(subregion,
                             levels = c("Mainland",
                                        "Continental islands",
                                        "Andaman & Nicobar",
                                        "Sunda islands",
                                        "Philippines",
                                        "Wallacea",
                                        "Papua"))))+
  coord_sf(ylim=c(-11,19),xlim=c(94,150.5))+
  #coord_sf(ylim=c(-20,20),xlim=c(92,170))+
  scale_fill_manual(values = c("#A1A1A190",
                               "#481F7070",
                               
                               "#fc8d6270",
                               "#8da0cb70",
                               "#66c2a570",
                               
                               "#0868ac70",
                               
                               "#a8ddb570"),
                    labels = c("Mainland",
                               "Continental islands",
                               "Andaman & Nicobar",
                               "Sunda islands",
                               "Philippines",
                               "Wallacea",
                               "New Guinea"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#481F70",
                                
                                "#fc8d62",
                                "#8da0cb",
                                "#66c2a5",
                                
                                "#0868ac",
                                
                                "#a8ddb5"),
                     labels = c("Mainland",
                                "Continental islands",
                                "Andaman & Nicobar",
                                "Sunda islands",
                                "Philippines",
                                "Wallacea",
                                "New Guinea"))+
  labs(x = "Longitude",
       y = "Latitude",
       title = "East Indo-Pacific",
       subtitle = "n = 447",
       fill = "",
       color = "")+
  theme_classic()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.9))+
  guides(color = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 2))
# Save maps
ggsave(filename = "figures jpg/OrInMa_Map_Hexagons.jpg", plot = mapOrInMa,
       dpi = 300, 
       width = 7, height = 5, units = "in")

# Caribbean
mapCarib <- ggplot()+
  geom_sf(data=world1)+
  geom_sf(data=wrapped_gridCaribbean, 
          aes(fill = factor(subregion,
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
  coord_sf(ylim=c(8,26.5),
           xlim=c(-92,-59))+
#  geom_hline(yintercept = 24)+geom_vline(xintercept = -79.5)+
  scale_fill_manual(values = c("#A1A1A190",
                               "#481F7070",
                               "#E3E41870",
                               "#35B77970",
                               "#21908C70"),
                    labels = c("Mainland",
                               "Continental islands",
                               "Bahamas \n(Lucayan)",
                               "Greater Antilles",
                               "Lesser Antilles \n(Kalinago)"))+
  scale_color_manual(values = c("#A1A1A1",
                                "#481F70",
                                "#E3E418",
                                "#35B779",
                                "#21908C"),
                     labels = c("Mainland",
                                "Continental islands",
                                "Bahamas \n(Lucayan)",
                                "Greater Antilles",
                                "Lesser Antilles \n(Kalinago)"))+
  labs(x = "Longitude",
       y = "Latitude",
       title = "Caribbean",
       subtitle = "n = 1682",
       fill = "",
       color = "")+
  theme_classic()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.9))+
  guides(color = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 2))


ggsave(filename = "Carib_Map_Hexagons.jpg", plot = mapCarib,
       dpi = 300, 
       width = 7, height = 5, units = "in")

