# Figure of the PCAs Functional Biogeography in three meta-archipelagos

library(tidyverse)
library(ggExtra) # for ggMarginal() of the density plot
library(gridExtra) #for grid.arrange()
library(ggrepel) # for geom_text_repel()

# Call the data 

islands_cwm_mu <- readRDS("Completeness_data_Islands/Islands_CWM_mu.rds")

spp_cell_df <- readRDS("Completeness_data_Islands/Species_in_cells_MetaArchipelagos.rds")

## Caribbean ####

caribbean_cwm <- islands_cwm_mu |>
  filter(Meta.Archipelago %in% "Neotropical")


funspaceDim(scale(log10(caribbean_cwm[,c(5:20)]))) # 4 dimensions

carib.pca.cwm <- prcomp(scale(log10(caribbean_cwm[,c(5:20)])))

summary(carib.pca.cwm)

carib.pca.values <- data.frame(cell = caribbean_cwm$cell,
                               region = caribbean_cwm$region,
                               subregion = caribbean_cwm$subregion,
                               fig_group = caribbean_cwm$fig_group,
                               sp_richness = caribbean_cwm$spp_richness,
                               carib.pca.cwm$x)

carib.pca.loadings <- data.frame(Variables = rownames(carib.pca.cwm$rotation), carib.pca.cwm$rotation)

#figure(s) ####
# Figure e - PC1-PC2 Caribbean ####
fig.e <- ggplot(carib.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 21) +
  scale_shape_manual(values = c(21,22))+
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
  scale_x_continuous(limits = c(-11,12.5))+
  scale_y_continuous(limits = c(-11,7))+
  coord_fixed()+
  geom_segment(data = carib.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC1*10, 
                   y = 0, yend = PC2*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = carib.pca.loadings, 
                  aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC1 (44.91%)",
       y = "PC2 (12.45%)",
#       title = expression(bold("e")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.e <- ggMarginal(fig.e,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Figure h - PC1-PC2 Caribbean ####
fig.h <- ggplot(carib.pca.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 21) +
  scale_shape_manual(values = c(21,22))+
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
  scale_x_continuous(limits = c(-11,7))+
  scale_y_continuous(limits = c(-6.5,8.5))+
  coord_fixed()+
  geom_segment(data = carib.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC2*10, 
                   y = 0, yend = PC3*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = carib.pca.loadings, 
                  aes(x = PC2*12, y = PC3*10, label = Variables)) + 
  labs(x = "PC2 (12.45%)",
       y = "PC3 (9.72%)",
#       title = expression(bold("h")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.h <- ggMarginal(fig.h,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Oriental-Indo-Malayan ####

orinma_cwm <- islands_cwm_mu |>
  filter(Meta.Archipelago %in% "Indo.Pacific",
         subregion %in% c("Andaman & Nicobar",
                          "Continental islands",
                          "Mainland",
                          "Philippines",
                          "Sunda islands")) |>
  # bring back the mean value of longitude for each cell
  left_join(spp_cell_df |> 
              dplyr::select(cell,longitude, latitude) |> 
              group_by(cell) |> 
              summarize(longitude = mean(longitude),
                        latitude = mean(latitude))) |>
  # extract the cells to the left of longitude 140 and north of -11 latitude
  filter(longitude < 140,
         latitude > -11)

funspaceDim(scale(log10(orinma_cwm[,c(5:20)]))) # 3 dimensions

orinma.pca.cwm <- prcomp(scale(log10(orinma_cwm[,c(5:20)])))

summary(orinma.pca.cwm)

orinma.pca.values <- data.frame(cell = orinma_cwm$cell,
                                region = orinma_cwm$region,
                                subregion = orinma_cwm$subregion,
                                fig_group = orinma_cwm$fig_group,
                                sp_richness = orinma_cwm$spp_richness,
                                orinma.pca.cwm$x)

orinma.pca.loadings <- data.frame(Variables = rownames(orinma.pca.cwm$rotation), orinma.pca.cwm$rotation)

# Figure f - PC1-PC2 OrInMa ####

fig.f <- ggplot(orinma.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 22) +
  scale_fill_manual(values = c("#D44292",
                               "#952EA0",
                               "#A1A1A1",
                               "#F6A97A",
                               "#F66D7A"))+
  scale_color_manual(values = c("#D44292",
                                "#952EA0",
                                "#A1A1A1",
                                "#F6A97A",
                                "#F66D7A"))+
  scale_x_continuous(limits = c(-11,12.5))+
  scale_y_continuous(limits = c(-11,7))+
  coord_fixed()+
  geom_segment(data = orinma.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC1*10, 
                   y = 0, yend = PC2*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = orinma.pca.loadings, 
                  aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC1 (47.53%)",
       y = "PC2 (14.62%)",
#       title = expression(bold("f")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.f <- ggMarginal(fig.f,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Figure i - PC2-PC3 OrInMa ####

fig.i <- ggplot(orinma.pca.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 22) +
  scale_fill_manual(values = c("#D44292",
                               "#952EA0",
                               "#A1A1A1",
                               "#F6A97A",
                               "#F66D7A"))+
  scale_color_manual(values = c("#D44292",
                                "#952EA0",
                                "#A1A1A1",
                                "#F6A97A",
                                "#F66D7A"))+
  scale_x_continuous(limits = c(-11,7))+
  scale_y_continuous(limits = c(-6.5,8.5))+
  coord_fixed()+
  geom_segment(data = orinma.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC2*10, 
                   y = 0, yend = PC3*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = orinma.pca.loadings, 
                  aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC2 (14.62%)",
       y = "PC3 (8.49%)",
#       title = expression(bold("i")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.i <- ggMarginal(fig.i,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Papuan-Melanesian ####

papuan_cwm <- islands_cwm_mu |>
  # extract subregions of Papuan-Melanesian region - including mainland for Australia
  filter(subregion %in% c("Bismarcks",
                          "Papua",
                          "Solomons",
                          "Vanuatu",
                          "Mainland",
                          "Wallacea")) |> 
  # bring back the mean value of longitude for each cell
  left_join(spp_cell_df |> 
              dplyr::select(cell,longitude, latitude) |> 
              group_by(cell) |> 
              summarize(longitude = mean(longitude),
                        latitude = mean(latitude))) |>
  # extract the cells to the right of longitude 120 (~ Sulawesi)
  filter(longitude >=110)

funspaceDim(scale(log10(papuan_cwm[,c(5:20)]))) # 3 dimensions

papuan.pca.cwm <- prcomp(scale(log10(papuan_cwm[,c(5:20)])))

summary(papuan.pca.cwm)

papuan.pca.values <- data.frame(cell = papuan_cwm$cell,
                                region = papuan_cwm$region,
                                subregion = papuan_cwm$subregion,
                                fig_group = papuan_cwm$fig_group,
                                sp_richness = papuan_cwm$spp_richness,
                                papuan.pca.cwm$x)

papuan.pca.loadings <- data.frame(Variables = rownames(papuan.pca.cwm$rotation), papuan.pca.cwm$rotation)

# Figure g - PC1-PC2 Papuan ####

fig.g <- ggplot(papuan.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 22) +
  scale_fill_manual(values = c("#0072B2",
                               "#A1A1A1",
                               "#CC79A7",
                               "#F0E442",
                               "#009E73",
                               "#D55E00"))+
  scale_color_manual(values = c("#0072B2",
                                "#A1A1A1",
                                "#CC79A7",
                                "#F0E442",
                                "#009E73",
                                "#D55E00"))+
  scale_x_continuous(limits = c(-11,12.5))+
  scale_y_continuous(limits = c(-11,7))+
  coord_fixed()+
  geom_segment(data = papuan.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC1*10, 
                   y = 0, yend = PC2*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = papuan.pca.loadings, 
                  aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC1 (48.90%)",
       y = "PC2 (15.12%)",
#       title = expression(bold("g")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.g <- ggMarginal(fig.g,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Figure j - PC2-PC3 Papuan ####

fig.j <- ggplot(papuan.pca.values) +
  geom_point(aes(x = PC2, 
                 y = PC3, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.3,
             shape = 22) +
  scale_fill_manual(values = c("#0072B2",
                               "#A1A1A1",
                               "#CC79A7",
                               "#F0E442",
                               "#009E73",
                               "#D55E00"))+
  scale_color_manual(values = c("#0072B2",
                                "#A1A1A1",
                                "#CC79A7",
                                "#F0E442",
                                "#009E73",
                                "#D55E00"))+
  scale_x_continuous(limits = c(-11,7))+
  scale_y_continuous(limits = c(-6.5,8.5))+
  coord_fixed()+
  geom_segment(data = papuan.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC1*10, 
                   y = 0, yend = PC2*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = papuan.pca.loadings, 
                  aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC2 (15.12%)",
       y = "PC3 (9.06%)",
#       title = expression(bold("j")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.j <- ggMarginal(fig.j,
                            type = "density",
                            size = 15,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Combine figure ####
Figure.PCAs <- grid.arrange(fig.density.e, fig.density.f, fig.density.g,
             fig.density.h, fig.density.i, fig.density.j, 
             ncol = 3)

ggsave(filename = "PCAs_FunctionalBiogeography.pdf", plot = Figure.PCAs,
       dpi = 600, 
       width = 330, height = 200, units = "mm")
