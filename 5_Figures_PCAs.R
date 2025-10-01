# Figure of the PCAs Functional Biogeography in three meta-archipelagos

library(tidyverse)
library(ggExtra) # for ggMarginal() of the density plot
library(gridExtra) #for grid.arrange()
library(ggrepel) # for geom_text_repel()
library(funspace) # functional trait space
library(ggstar) # points as hexagons

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
  mutate(Meta_Archipelago = ifelse(Meta.Archipelago %in% "Indo.Pacific" &
                                     subregion %in% c("Andaman & Nicobar",
                                                      "Continental islands",
                                                      "Mainland",
                                                      "Philippines",
                                                      "Sunda islands"),
                                   "Or-In-Ma",
                                   ifelse(Meta.Archipelago %in% "Indo.Pacific" &
                                            subregion %in% c("Bismarcks",
                                                             "Papua",
                                                             "Solomons",
                                                             "Vanuatu",
                                                             "Mainland",
                                                             "Wallacea"),
                                          "Papuan",
                                          "Caribbean"))) |> 
  left_join(islands_SppR_Pred_SppR)

spp_cell_df <- readRDS("Completeness_data_Islands/Species_in_cells_MetaArchipelagos.rds")

# PCA Caribbean CWM ####

caribbean_cwm <- islands_cwm_mu_sem |>
  filter(Meta_Archipelago %in% "Caribbean")

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

carib.pca.cwm.loadings <- data.frame(Variables = rownames(carib.pca.cwm$rotation), carib.pca.cwm$rotation)

## Figure CWM-a - PC1-PC2 Caribbean ####
summary(carib.pca.cwm.values[,8:10])
summary(carib.pca.cwm)

cwm.a <- ggplot(carib.pca.cwm.values) +
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
#  geom_segment(data = carib.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*5, 
#                   y = 0, yend = PC2*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwm.loadings, 
                  aes(x = PC1*12, y = PC2*12, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC1 (57.83%)",
       y = "PC2 (13.47%)",
       title = expression(bold("a")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwm.density.a <- ggMarginal(cwm.a,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure CWM-c - PC2-PC3 Caribbean ####
summary(carib.pca.cwm.values[,8:10])
summary(carib.pca.cwm)


cwm.c <- ggplot(carib.pca.cwm.values) +
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
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
#  geom_segment(data = carib.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC2*5, 
#                   y = 0, yend = PC3*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwm.loadings, 
                  aes(x = PC2*12, y = PC3*12, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC2 (13.47%)",
       y = "PC3 (12.04%)",
       title = expression(bold("c")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwm.density.c <- ggMarginal(cwm.c,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# PCA Oriental-Indo-Malayan CWM ####

orinma_cwm <- islands_cwm_mu_sem |>
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
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*5, 
#                   y = 0, yend = PC2*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwm.loadings, 
                  aes(x = PC1*12, y = PC2*12, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC1 (60.46%)",
       y = "PC2 (15.24%)",
       title = expression(bold("b")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwm.density.b <- ggMarginal(cwm.b,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

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
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwm.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC2*5, 
#                   y = 0, yend = PC3*4),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwm.loadings, 
                  aes(x = PC2*12, y = PC3*12, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC2 (15.24%)",
       y = "PC3 (8.77%)",
       title = expression(bold("d")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

cwm.density.d <- ggMarginal(cwm.d,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Combine figure CWM a-d ####

Figure.CWM <- grid.arrange(cwm.density.a, cwm.density.b, 
                            cwm.density.c, cwm.density.d, 
                            ncol = 2)

ggsave(filename = "PCAs_CWM.pdf", plot = Figure.CWM,
       dpi = 600, 
       width = 170, height = 170, units = "mm")





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

## Figure e - CWV-PC1-PC2 Caribbean ####
summary(carib.pca.cwv.values[,8:11])
summary(carib.pca.cwv)

fig.e <- ggplot(carib.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-8.6,8.9))+
  scale_y_continuous(limits = c(-6.3,9.4))+
  coord_fixed()+
#  geom_segment(data = carib.pca.cwv.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*10, 
#                   y = 0, yend = PC2*8),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC1*10, y = PC2*10, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC1 (44.93%)",
       y = "PC2 (13.61%)",
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
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure g - CWV-PC2-PC3 Caribbean ####
summary(carib.pca.cwv)

fig.g <- ggplot(carib.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-6.3,9.4))+
  scale_y_continuous(limits = c(-5.8,5.6))+
  coord_fixed()+
  #  geom_segment(data = carib.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC2*10, y = PC3*10, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC2 (13.61%)",
       y = "PC3 (11.92%)",
       #       title = expression(bold("e")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.g <- ggMarginal(fig.g,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure i - CWV-PC3-PC4 Caribbean ####
summary(carib.pca.cwv)

fig.i <- ggplot(carib.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-5.8,5.6))+
  scale_y_continuous(limits = c(-5.1,4.2))+
  coord_fixed()+
  #  geom_segment(data = carib.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = carib.pca.cwv.loadings, 
                  aes(x = PC3*6, y = PC4*6, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC3 (11.92%)",
       y = "PC4 (7.94%)",
       #       title = expression(bold("e")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.i <- ggMarginal(fig.i,
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

## Figure f - CWV-PC1-PC2 OrInMa ####
summary(orinma.pca.cwv.values[,8:11])
summary(orinma.pca.cwv)

fig.f <- ggplot(orinma.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-8.6,8.9))+
  scale_y_continuous(limits = c(-6.3,9.4))+
  coord_fixed()+
#  geom_segment(data = orinma.pca.cwv.loadings, 
#               size = 0.25,
#               aes(x = 0, xend = PC1*10, 
#                   y = 0, yend = PC2*8),
#               arrow = arrow(length = unit(0.1, "cm")),
#               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC1*10, y = PC2*10, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC1 (46.06%)",
       y = "PC2 (14.23%)",
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
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure h - CWV-PC2-PC3 OrInMa ####
summary(orinma.pca.cwv)

fig.h <- ggplot(orinma.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-6.3,9.4))+
  scale_y_continuous(limits = c(-5.8,5.6))+
  coord_fixed()+
  #  geom_segment(data = orinma.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC1*10, y = PC2*10, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC2 (14.23%)",
       y = "PC3 (11.03%)",
       #       title = expression(bold("f")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.h <- ggMarginal(fig.h,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

## Figure j - CWV-PC3-PC4 OrInMa ####
summary(orinma.pca.cwv)

fig.j <- ggplot(orinma.pca.cwv.values) +
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
  scale_x_continuous(limits = c(-5.8,5.6))+
  scale_y_continuous(limits = c(-5.1,4.2))+
  coord_fixed()+
  #  geom_segment(data = orinma.pca.cwv.loadings, 
  #               size = 0.25,
  #               aes(x = 0, xend = PC1*10, 
  #                   y = 0, yend = PC2*8),
  #               arrow = arrow(length = unit(0.1, "cm")),
  #               colour = "black") +
  geom_text_repel(data = orinma.pca.cwv.loadings, 
                  aes(x = PC1*6, y = PC2*6, 
                      label = Variables),
                  size = 3.5) + 
  labs(x = "PC3 (11.03%)",
       y = "PC4 (8.48%)",
       #       title = expression(bold("f")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

fig.density.j <- ggMarginal(fig.j,
                            type = "density",
                            size = 10,
                            groupColour = TRUE, 
                            groupFill = TRUE)

# Combine figure ####
Figure.PCAs <- grid.arrange(fig.density.a, fig.density.b, 
                            fig.density.c, fig.density.d, 
                            fig.density.e, fig.density.f, 
                            fig.density.g, fig.density.h,
                            fig.density.i, fig.density.j,
             ncol = 2)

ggsave(filename = "PCAs_FunctionalBiogeography.pdf", plot = Figure.PCAs,
       dpi = 600, 
       width = 178, height = 340, units = "mm")

# SEM for mechanisms ####

islands_cwm_mu_sem |> 
  ggplot(aes(x = SppRichnessCell, 
             y = PredSppRichnessCell,
             color = factor(subregion,
                            levels = c("Bahamas (Lucayan)",
                                       "Lesser Antilles (Kalinago)",
                                       "Greater Antilles",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Papua",
                                       "Wallacea",
                                       #"Bismarcks",
                                       "Solomons",
                                       "Vanuatu",
                                       "Continental islands",
                                       "Mainland"))))+
  facet_wrap(fig_group~Meta_Archipelago, scale = "free") +
  geom_point(aes(y = -0.5), shape = ".")+
  geom_point(alpha = 0.4) +
  labs(x = "Species richness",
       y = "(bird) Predator species richness",
       color = "Subregion \n (Archipelago)") +
  scale_color_manual(values = c("#E3E418",
                                "#35B779",
                                "#21908C",
                                "#D44292",
                                "#F66D7A",
                                "#F6A97A",
                                "#CC79A7",
                                "#D55E00",
                                #"#0072B2",
                                "#F0E442",
                                "#009E73",
                                "#952EA0",
                                "#A1A1A1"))+
  geom_smooth(method = "lm", se = TRUE) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))

# bring back the PCA values
pca_values <- bind_rows(carib.pca.values,
                        orinma.pca.values,
                        papuan.pca.values)

head(pca_values[,1:10])

pca_values <- pca_values |>
  filter(!is.na(pred_richness))

pca_values |> 
  dplyr::select(cell, Meta_Archipelago, subregion, sp_richness, pred_richness, PC1, PC2, PC3) |>
  filter(Meta_Archipelago != "Papuan") |>
  pivot_longer(cols = c(PC1, PC2, PC3), names_to = "PC dimension", values_to = "PC value") |>
  ggplot(aes(x = sp_richness, 
             y = `PC value`,
             color = factor(subregion,
                            levels = c("Bahamas (Lucayan)",
                                       "Lesser Antilles (Kalinago)",
                                       "Greater Antilles",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Papua",
                                       "Wallacea",
                                       #"Bismarcks",
                                       #"Solomons",
                                       "Vanuatu",
                                       "Continental islands",
                                       "Mainland"))))+
  facet_grid(`PC dimension`~Meta_Archipelago, 
             scales = "free") +
  geom_point(alpha = 0.2)+
  scale_color_manual(values = c("#E3E418",
                                "#35B779",
                                "#21908C",
                                "#D44292",
                                "#F66D7A",
                                "#F6A97A",
                                #"#CC79A7",
                                #"#D55E00",
                                #"#0072B2",
                                #"#F0E442",
                                #"#009E73",
                                "#952EA0",
                                "#A1A1A1"))+
  labs(x = "Species richness",
       y = "PC values",
       color = "Subregion \n (Archipelago)")+
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))

pca_values |> 
  dplyr::select(cell, Meta_Archipelago, subregion, sp_richness, pred_richness, PC1, PC2, PC3) |>
  filter(Meta_Archipelago != "Papuan") |>
  pivot_longer(cols = c(PC1, PC2, PC3), names_to = "PC dimension", values_to = "PC value") |>
  ggplot(aes(x = pred_richness, 
             y = `PC value`,
             color = factor(subregion,
                            levels = c("Bahamas (Lucayan)",
                                       "Lesser Antilles (Kalinago)",
                                       "Greater Antilles",
                                       "Andaman & Nicobar",
                                       "Sunda islands",
                                       "Philippines",
                                       "Papua",
                                       "Wallacea",
                                       #"Bismarcks",
                                       #"Solomons",
                                       "Vanuatu",
                                       "Continental islands",
                                       "Mainland"))))+
  facet_grid(`PC dimension`~Meta_Archipelago, 
             scales = "free") +
  geom_point(alpha = 0.2)+
  scale_color_manual(values = c("#E3E418",
                                "#35B779",
                                "#21908C",
                                "#D44292",
                                "#F66D7A",
                                "#F6A97A",
                                #"#CC79A7",
                                #"#D55E00",
                                #"#0072B2",
                                #"#F0E442",
                                #"#009E73",
                                "#952EA0",
                                "#A1A1A1"))+
  labs(x = "Predator species richness (birds)",
       y = "PC values",
       color = "Subregion \n (Archipelago)")+
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))





# this should be for each meta archipelago independently

library(piecewiseSEM)
library(lme4)

pca_carib_sem <- pca_values |> 
  filter(Meta_Archipelago == "Caribbean")

pca_orinma_sem <- pca_values |> 
  filter(Meta_Archipelago == "Or-In-Ma")

mod_PC1_car <- lmer(PC1 ~ sp_richness + pred_richness + (1 | subregion), data = pca_carib_sem)
mod_PC2_car <- lmer(PC2 ~ sp_richness + pred_richness + (1 | subregion), data = pca_carib_sem)
mod_PC3_car <- lmer(PC3 ~ sp_richness + pred_richness + (1 | subregion), data = pca_carib_sem)

sem_pca_car <- psem(mod_PC1_car, mod_PC2_car, mod_PC3_car)
summary_sem_pca_car <- summary(sem_pca_car)

saveRDS(summary_sem_pca_car, "Summary_SEM_PCA_Caribbean.rds")

mod_PC1_oim <- lmer(PC1 ~ sp_richness + pred_richness + (1 | subregion), data = pca_orinma_sem)
mod_PC2_oim <- lmer(PC2 ~ sp_richness + pred_richness + (1 | subregion), data = pca_orinma_sem)
mod_PC3_oim <- lmer(PC3 ~ sp_richness + pred_richness + (1 | subregion), data = pca_orinma_sem)

sem_pca_oim <- psem(mod_PC1_oim, mod_PC2_oim, mod_PC3_oim)
summary_sem_pca_oim <- summary(sem_pca_oim)

saveRDS(summary_sem_pca_oim, "Summary_SEM_PCA_OrInMa.rds")
