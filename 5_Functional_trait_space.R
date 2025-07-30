# Functional trait space (PCAs)

# Packages ####

# Data manipulation and figures
library(tidyverse)
library(ggExtra) # for ggMarginal() of the density plot
library(gridExtra) #for grid.arrange()
library(ggrepel) # for geom_text_repel()
# gis
library(sf)
# phylogeny
library(phytools); library(ape)
# functional trait space
library(funspace)
# correlation plot
library(corrplot)

### bring back the data
spp_cell_df_func <- readRDS("Completeness_data_Islands/Records_eBird_func_traits.rds")

names(spp_cell_df_func)

# select only columns of traits to do PCA of the Community Weighted Means & Variance

# function to compute the community weighted variance
cwv <- function(x, wt, CWM){
  sum(wt*(x-CWM)^2)
}

# community weighted mean can be computed with `weighted.means()`

# Calculating the CWM and CWV

islands_cwm_mu_summary <- spp_cell_df_func |>
  group_by(cell, Meta.Archipelago, subregion) |>
  summarise(spp_richness = n(),
            # CWM
            # dispersal morphology
            `body mass` = weighted.mean(bod.mas, mu_hat, na.rm = T),
            `hand wing index` = weighted.mean(hwi, mu_hat, na.rm = T),
            # foraging morphology
            `tarsus length` = weighted.mean(tar.len, mu_hat, na.rm = T),
            `wing length` = weighted.mean(win.len, mu_hat, na.rm = T),
            # dietary morphology
            `beak length` = weighted.mean(bea.len, mu_hat, na.rm = T), 
            `beak width` = weighted.mean(bea.wid, mu_hat, na.rm = T),
            `beak depth` = weighted.mean(bea.dep, mu_hat, na.rm = T), 
            # ecological niche
            `body length` = weighted.mean(bod.len, mu_hat, na.rm = T), #highly correlated
            `habitat breadth` = weighted.mean(hab.bre, mu_hat, na.rm = T),
            `range size` = weighted.mean(ran.siz, mu_hat, na.rm = T),
            # converting 0-1 index to 1-2 (log10 transformation)
            `endemic insularity` = (weighted.mean(end.ins, mu_hat, na.rm = T)+1),
            # foraging niche 
            # converting 0-1 index to 1-2 (log10 transformation)
            pred = (weighted.mean(predatory, mu_hat, na.rm = T)+1),
            # converting 0-1 index to 1-2 (log10 transformation)
            `diet generalist` = (weighted.mean(generalist, mu_hat, na.rm = T)+1),
            # behavioral / reproductive niche
            verticality = weighted.mean(vertical, mu_hat, na.rm = T),
            `nest structure breadth` = weighted.mean(nes.str.bre, mu_hat, na.rm = T),
            `nest location breadth` = weighted.mean(nes.loc.bre, mu_hat, na.rm = T),
            # CWV  
            # dispersal morphology 
            `bmass cwv` = cwv(bod.mas, mu_hat, `body mass`), 
            `hwi cwv` = cwv(hwi, mu_hat, `hand wing index`),
            # foraging morphology
            `tarsus cwv` = cwv(tar.len, mu_hat, `tarsus length`),
            `wing cwv` = cwv(win.len, mu_hat, `wing length`),
            # dietary morphology
            `beak l cwv` = cwv(bea.len, mu_hat, `beak length`), 
            `beak w cwv` = cwv(bea.wid, mu_hat, `beak width`),
            `beak d cwv` = cwv(bea.dep, mu_hat, `beak depth`), 
            # ecological niche
            `body l cwv` = cwv(bod.len, mu_hat, `body length`),
            `habitat b cwv` = cwv(hab.bre, mu_hat, `habitat breadth`),
            `range cwv` = cwv(ran.siz, mu_hat, `range size`),
            # converting 0-1 index to 1-2 (log10 transformation)
            `insularity cwv` = (cwv(end.ins, mu_hat, `endemic insularity`)),
            # foraging niche 
            # converting 0-1 index to 1-2 (log10 transformation)
            `predatory cwv` = (cwv(predatory, mu_hat, pred)),
            # converting 0-1 index to 1-2 (log10 transformation)
            `generalist cwv` = (cwv(generalist, mu_hat, `diet generalist`)),
            # behavioral / reproductive niche
            `verticality cwv` = cwv(vertical, mu_hat, verticality),
            `nest str cwv` = cwv(nes.str.bre, mu_hat, `nest structure breadth`),
            `nest loc cwv` = cwv(nes.loc.bre, mu_hat, `nest location breadth`),
  ) # 17 spp is the 1st quartile, but this reduces a lot 

islands_cwm_omega_summary <- spp_cell_df_func |>
  group_by(cell, Meta.Archipelago, subregion) |>
  summarise(spp_richness = n(),
            # CWM
            # dispersal morphology
            `body mass` = weighted.mean(bod.mas, omega_hat, na.rm = T),
            `hand wing index` = weighted.mean(hwi, omega_hat, na.rm = T),
            # foraging morphology
            `tarsus length` = weighted.mean(tar.len, omega_hat, na.rm = T),
            `wing length` = weighted.mean(win.len, omega_hat, na.rm = T),
            # dietary morphology
            `beak length` = weighted.mean(bea.len, omega_hat, na.rm = T), 
            `beak width` = weighted.mean(bea.wid, omega_hat, na.rm = T),
            `beak depth` = weighted.mean(bea.dep, omega_hat, na.rm = T), 
            # ecological niche
            `body length` = weighted.mean(bod.len, omega_hat, na.rm = T), #highly correlated
            `habitat breadth` = weighted.mean(hab.bre, omega_hat, na.rm = T),
            `range size` = weighted.mean(ran.siz, omega_hat, na.rm = T),
            # converting 0-1 index to 1-2 (log10 transformation)
            `endemic insularity` = (weighted.mean(end.ins, omega_hat, na.rm = T)+1),
            # foraging niche 
            # converting 0-1 index to 1-2 (log10 transformation)
            pred = (weighted.mean(predatory, omega_hat, na.rm = T)+1),
            # converting 0-1 index to 1-2 (log10 transformation)
            `diet generalist` = (weighted.mean(generalist, omega_hat, na.rm = T)+1),
            # behavioral / reproductive niche
            verticality = weighted.mean(vertical, omega_hat, na.rm = T),
            `nest structure breadth` = weighted.mean(nes.str.bre, omega_hat, na.rm = T),
            `nest location breadth` = weighted.mean(nes.loc.bre, omega_hat, na.rm = T),
            # CWV  
            # dispersal morphology 
            `bmass cwv` = cwv(bod.mas, omega_hat, `body mass`), 
            `hwi cwv` = cwv(hwi, omega_hat, `hand wing index`),
            # foraging morphology
            `tarsus cwv` = cwv(tar.len, omega_hat, `tarsus length`),
            `wing cwv` = cwv(win.len, omega_hat, `wing length`),
            # dietary morphology
            `beak l cwv` = cwv(bea.len, omega_hat, `beak length`), 
            `beak w cwv` = cwv(bea.wid, omega_hat, `beak width`),
            `beak d cwv` = cwv(bea.dep, omega_hat, `beak depth`), 
            # ecological niche
            `body l cwv` = cwv(bod.len, omega_hat, `body length`),
            `habitat b cwv` = cwv(hab.bre, omega_hat, `habitat breadth`),
            `range cwv` = cwv(ran.siz, omega_hat, `range size`),
            # converting 0-1 index to 1-2 (log10 transformation)
            `insularity cwv` = (cwv(end.ins, omega_hat, `endemic insularity`)),
            # foraging niche 
            # converting 0-1 index to 1-2 (log10 transformation)
            `predatory cwv` = (cwv(predatory, omega_hat, pred)),
            # converting 0-1 index to 1-2 (log10 transformation)
            `generalist cwv` = (cwv(generalist, omega_hat, `diet generalist`)),
            # behavioral / reproductive niche
            `verticality cwv` = cwv(vertical, omega_hat, verticality),
            `nest str cwv` = cwv(nes.str.bre, omega_hat, `nest structure breadth`),
            `nest loc cwv` = cwv(nes.loc.bre, omega_hat, `nest location breadth`),
  ) # 17 spp is the 1st quartile

# group for figure
islands_cwm_mu <- islands_cwm_mu_summary |>
  mutate(subregion = ifelse(subregion == "Mainland", "Mainland",
                            ifelse(subregion == "Caribbean - Continental islands", "Continental islands",
                                   ifelse(subregion == "Caribbean - Greater Antilles", "Greater Antilles",
                                          ifelse(subregion == "Caribbean - Kalinago", "Lesser Antilles (Kalinago)",
                                                 ifelse(subregion == "Caribbean - Lucayan", "Bahamas (Lucayan)",
                                                        ifelse(subregion == "Oriental-Indo-Malayan - Andaman & Nicobar", "Andaman & Nicobar",
                                                               ifelse(subregion == "Oriental-Indo-Malayan - Continental islands", "Continental islands",
                                                                      ifelse(subregion == "Oriental-Indo-Malayan - Philippines", "Philippines",
                                                                             ifelse(subregion == "Oriental-Indo-Malayan - Sunda islands", "Sunda islands",
                                                                                    ifelse(subregion == "Oriental-Indo-Malayan - Wallacea", "Wallacea",
                                                                                           ifelse(subregion == "Papuan-Melanesian - Bismarcks", "Bismarcks",
                                                                                                  ifelse(subregion == "Papuan-Melanesian - Papua", "Papua",
                                                                                                         ifelse(subregion == "Papuan-Melanesian - Solomons", "Solomons",
                                                                                                                ifelse(subregion == "Papuan-Melanesian - Vanuatu", "Vanuatu",
                                                                                                                       subregion)))))))))))))),
         region = ifelse(Meta.Archipelago == "Indo.Pacific", "Indo Pacific",
                         "Neotropical"),
         fig_group = ifelse(subregion == "Mainland", "Mainland",
                            ifelse(subregion == "Continental islands", "Islands",
                                   ifelse(subregion == "Greater Antilles", "Islands",
                                          ifelse(subregion == "Lesser Antilles (Kalinago)", "Islands",
                                                 ifelse(subregion == "Bahamas (Lucayan)", "Islands",
                                                        ifelse(subregion == "Andaman & Nicobar", "Islands",
                                                               ifelse(subregion == "Philippines", "Islands",
                                                                      ifelse(subregion == "Sunda islands", "Islands",
                                                                             ifelse(subregion == "Wallacea", "Islands",
                                                                                    ifelse(subregion == "Bismarcks", "Islands",
                                                                                           ifelse(subregion == "Papua", "Islands",
                                                                                                  ifelse(subregion == "Solomons", "Islands",
                                                                                                         ifelse(subregion == "Vanuatu", "Islands",
                                                                                                                subregion))))))))))))))

islands_cwm_omega <- islands_cwm_omega_summary |>
  mutate(subregion = ifelse(subregion == "Mainland", "Mainland",
                            ifelse(subregion == "Caribbean - Continental islands", "Continental islands",
                                   ifelse(subregion == "Caribbean - Greater Antilles", "Greater Antilles",
                                          ifelse(subregion == "Caribbean - Kalinago", "Lesser Antilles (Kalinago)",
                                                 ifelse(subregion == "Caribbean - Lucayan", "Bahamas (Lucayan)",
                                                        ifelse(subregion == "Oriental-Indo-Malayan - Andaman & Nicobar", "Andaman & Nicobar",
                                                               ifelse(subregion == "Oriental-Indo-Malayan - Continental islands", "Continental islands",
                                                                      ifelse(subregion == "Oriental-Indo-Malayan - Philippines", "Philippines",
                                                                             ifelse(subregion == "Oriental-Indo-Malayan - Sunda islands", "Sunda islands",
                                                                                    ifelse(subregion == "Oriental-Indo-Malayan - Wallacea", "Wallacea",
                                                                                           ifelse(subregion == "Papuan-Melanesian - Bismarcks", "Bismarcks",
                                                                                                  ifelse(subregion == "Papuan-Melanesian - Papua", "Papua",
                                                                                                         ifelse(subregion == "Papuan-Melanesian - Solomons", "Solomons",
                                                                                                                ifelse(subregion == "Papuan-Melanesian - Vanuatu", "Vanuatu",
                                                                                                                       subregion)))))))))))))),
         region = ifelse(Meta.Archipelago == "Indo.Pacific", "Indo Pacific",
                         "Neotropical"),
         fig_group = ifelse(subregion == "Mainland", "Mainland",
                            ifelse(subregion == "Continental islands", "Islands",
                                   ifelse(subregion == "Greater Antilles", "Islands",
                                          ifelse(subregion == "Lesser Antilles (Kalinago)", "Islands",
                                                 ifelse(subregion == "Bahamas (Lucayan)", "Islands",
                                                        ifelse(subregion == "Andaman & Nicobar", "Islands",
                                                               ifelse(subregion == "Philippines", "Islands",
                                                                      ifelse(subregion == "Sunda islands", "Islands",
                                                                             ifelse(subregion == "Wallacea", "Islands",
                                                                                    ifelse(subregion == "Bismarcks", "Islands",
                                                                                           ifelse(subregion == "Papua", "Islands",
                                                                                                  ifelse(subregion == "Solomons", "Islands",
                                                                                                         ifelse(subregion == "Vanuatu", "Islands",
                                                                                                                subregion))))))))))))))

names(islands_cwm_mu)
names(islands_cwm_omega)

corrplot(round(cor(scale(log10(islands_cwm_mu[,c(5:20)]))),2), 
         type="upper", order="hclust", 
         tl.col="black", tl.srt=45)

corrplot(round(cor(scale(log10(islands_cwm_omega[,c(5:20)]))),2), 
         type="upper", order="hclust", 
         tl.col="black", tl.srt=45)

# Global CWM pca ####

funspaceDim(scale(log10(islands_cwm_mu[,c(5:20)]))) # 4 dimensions
funspaceDim(scale(log10(islands_cwm_omega[,c(5:20)]))) # 4 dimensions

pca.cwm.mu <- prcomp(scale(log10(islands_cwm_mu[,c(5:20)])))
pca.cwm.omega <- prcomp(scale(log10(islands_cwm_omega[,c(5:20)])))

summary(pca.cwm.mu)
summary(pca.cwm.omega)

pca.values.mu <- data.frame(cell = islands_cwm_mu$cell,
                            region = islands_cwm_mu$region,
                            subregion = islands_cwm_mu$subregion,
                            fig_group = islands_cwm_mu$fig_group,
                            sp_richness = islands_cwm_mu$spp_richness,
                            pca.cwm.mu$x)

pca.values.omega <- data.frame(cell = islands_cwm_omega$cell,
                               region = islands_cwm_omega$region,
                               subregion = islands_cwm_omega$subregion,
                               fig_group = islands_cwm_omega$fig_group,
                               sp_richness = islands_cwm_omega$spp_richness,
                               pca.cwm.omega$x)

pca.loadings.mu <- data.frame(Variables = rownames(pca.cwm.mu$rotation), 
                              pca.cwm.mu$rotation)

pca.loadings.omega <- data.frame(Variables = rownames(pca.cwm.omega$rotation), 
                                 pca.cwm.omega$rotation)

# to make sense in the ordination
pca.values$PC2 <- (pca.values$PC2*-1)
pca.loadings$PC2 <- (pca.loadings$PC2*-1)

#figure
pca.fig.mu.a <- ggplot(pca.values.mu) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 shape = region,
                 fill = fig_group,
                 color = fig_group), 
             alpha = 0.25) +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("#2980B9","gray"))+
  scale_color_manual(values = c("#2980B9","gray"))+
  scale_x_continuous(limits = c(-9.5,8.5))+
  scale_y_continuous(limits = c(-5.5,5.5))+
  coord_fixed()+
  geom_segment(data = pca.loadings.mu, 
               linewidth = 0.25,
               aes(x = 0, xend = PC1*6, 
                   y = 0, yend = PC2*4),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = pca.loadings.mu, 
                  aes(x = PC1*7, y = PC2*5, 
                      label = Variables)) + 
  labs(x = "PC1 (48.12%)",
       y = "PC2 (16.02%)",
       title = expression(bold("a ")~hat(mu)),
       fill = "Geographic setting",
       color = "Geographic setting",
       shape = "Region")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape=guide_legend(nrow=2,byrow=TRUE))

pca.density.fig.mu.a <- ggMarginal(pca.fig.mu.a,
                                   type = "density",
                                   groupColour = TRUE, 
                                   groupFill = TRUE)

pca.fig.mu.b <- ggplot(pca.values.mu) +
  geom_point(aes(x = PC3, 
                 y = PC4, 
                 shape = region,
                 fill = fig_group,
                 color = fig_group), 
             alpha = 0.25) +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("#2980B9","gray"))+
  scale_color_manual(values = c("#2980B9","gray"))+
  scale_x_continuous(limits = c(-9.5,8.5))+
  scale_y_continuous(limits = c(-5.5,5.5))+
  coord_fixed()+
  geom_segment(data = pca.loadings.mu, 
               linewidth = 0.25,
               aes(x = 0, xend = PC3*6, 
                   y = 0, yend = PC4*4),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = pca.loadings.mu, 
                  aes(x = PC3*7, y = PC4*5, 
                      label = Variables)) + 
  labs(x = "PC3 (8.34%)",
       y = "PC4 (7.31%)",
       title = expression(bold("b ")~hat(mu)),
       fill = "Geographic setting",
       color = "Geographic setting",
       shape = "Region")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape=guide_legend(nrow=2,byrow=TRUE))

pca.density.fig.mu.b <- ggMarginal(pca.fig.mu.b,
                                   type = "density",
                                   groupColour = TRUE, 
                                   groupFill = TRUE)


pca.fig.omega.a <- ggplot(pca.values.omega) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 shape = region,
                 fill = fig_group,
                 color = fig_group), 
             alpha = 0.25) +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("#2980B9","gray"))+
  scale_color_manual(values = c("#2980B9","gray"))+
  scale_x_continuous(limits = c(-9.5,8.5))+
  scale_y_continuous(limits = c(-5.5,5.5))+
  coord_fixed()+
  geom_segment(data = pca.loadings.omega, 
               linewidth = 0.25,
               aes(x = 0, xend = PC1*6, 
                   y = 0, yend = PC2*4),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = pca.loadings.omega, 
                  aes(x = PC1*7, y = PC2*5, 
                      label = Variables)) + 
  labs(x = "PC1 (46.6%)",
       y = "PC2 (15.82%)",
       title = expression(bold("a ")~hat(omega)),
       fill = "Geographic setting",
       color = "Geographic setting",
       shape = "Region")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape=guide_legend(nrow=2,byrow=TRUE))

pca.density.fig.omega.a <- ggMarginal(pca.fig.omega.a,
                                      type = "density",
                                      groupColour = TRUE, 
                                      groupFill = TRUE)

pca.fig.omega.b <- ggplot(pca.values.omega) +
  geom_point(aes(x = PC3, 
                 y = PC4, 
                 shape = region,
                 fill = fig_group,
                 color = fig_group), 
             alpha = 0.25) +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("#2980B9","gray"))+
  scale_color_manual(values = c("#2980B9","gray"))+
  scale_x_continuous(limits = c(-9.5,8.5))+
  scale_y_continuous(limits = c(-5.5,5.5))+
  coord_fixed()+
  geom_segment(data = pca.loadings.omega, 
               linewidth = 0.25,
               aes(x = 0, xend = PC3*6, 
                   y = 0, yend = PC4*4),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = pca.loadings.omega, 
                  aes(x = PC3*7, y = PC4*5, 
                      label = Variables)) + 
  labs(x = "PC3 (8.54%)",
       y = "PC4 (7.19%)",
       title = expression(bold("b ")~hat(omega)),
       fill = "Geographic setting",
       color = "Geographic setting",
       shape = "Region")+
  theme(legend.position = "bottom",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape=guide_legend(nrow=2,byrow=TRUE))

pca.density.fig.omega.b <- ggMarginal(pca.fig.omega.b,
                                      type = "density",
                                      groupColour = TRUE, 
                                      groupFill = TRUE)

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

#figure

global.carib.pca.fig <- ggplot(carib.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.25,
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
  scale_x_continuous(limits = c(-10.8,12.3))+
  scale_y_continuous(limits = c(-10.9,7))+
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
       title = expression(bold("(d)")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

global.carib.pca.density.fig <- ggMarginal(global.carib.pca.fig,
                                           type = "density",
                                           size = 15,
                                           groupColour = TRUE, 
                                           groupFill = TRUE)

### Caribbean per group ####

bahamas.pca.fig <- ggplot(carib.pca.values |> 
                            filter(subregion %in% "Bahamas (Lucayan)")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#FDE725",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#FDE725",
                               "#73D055",
                               "#33638D",
                               "#440154",
                               "gray65"))+
  scale_color_manual(values = c("#FDE725",
                                "#73D055",
                                "#33638D",
                                "#440154",
                                "gray65"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (28.93%)",
       y = "PC2 (23.59%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

bahamas.pca.density.fig <- ggMarginal(bahamas.pca.fig,
                                      type = "density",
                                      groupColour = TRUE, 
                                      groupFill = TRUE)

coniscarib.pca.fig <- ggplot(carib.pca.values |> 
                               filter(subregion %in% "Continental islands")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#73D055",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#73D055"))+
  scale_color_manual(values = c("#73D055",
                                "#33638D",
                                "#440154",
                                "gray65"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (28.93%)",
       y = "PC2 (23.59%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

coniscarib.pca.density.fig <- ggMarginal(coniscarib.pca.fig,
                                         type = "density",
                                         groupColour = TRUE, 
                                         groupFill = TRUE)

#greater Antilles

greant.pca.fig <- ggplot(carib.pca.values |> 
                           filter(subregion %in% "Greater Antilles")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#33638D",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#33638D"))+
  scale_color_manual(values = c("#33638D",
                                "#440154",
                                "gray65"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (28.93%)",
       y = "PC2 (23.59%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

greant.pca.density.fig <- ggMarginal(greant.pca.fig,
                                     type = "density",
                                     groupColour = TRUE, 
                                     groupFill = TRUE)

# lesser
lesant.pca.fig <- ggplot(carib.pca.values |> 
                           filter(subregion %in% "Lesser Antilles (Kalinago)")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#440154",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#440154"))+
  scale_color_manual(values = c("#440154",
                                "gray65"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (28.93%)",
       y = "PC2 (23.59%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

lesant.pca.density.fig <- ggMarginal(lesant.pca.fig,
                                     type = "density",
                                     groupColour = TRUE, 
                                     groupFill = TRUE)

mainlcarib.pca.fig <- ggplot(carib.pca.values |> 
                               filter(subregion %in% "Mainland")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "gray95",
                 alpha = 0.4)+
  scale_fill_manual(values = c("gray65"))+
  scale_color_manual(values = c("gray65"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (28.93%)",
       y = "PC2 (23.59%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

mainlcarib.pca.density.fig <- ggMarginal(mainlcarib.pca.fig,
                                         type = "density",
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
  left_join(spp_cell_df_func |> 
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

#figure

global.orinma.pca.fig <- ggplot(orinma.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.25,
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
  scale_x_continuous(limits = c(-10.8,12.3))+
  scale_y_continuous(limits = c(-10.9,7))+
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
       title = expression(bold("(e)")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

global.orinma.pca.density.fig <- ggMarginal(global.orinma.pca.fig,
                                            type = "density",
                                            size = 15,
                                            groupColour = TRUE, 
                                            groupFill = TRUE)

### Oriental-Indo-Malayan per group ####

an.pca.fig <- ggplot(orinma.pca.values |> 
                       filter(subregion %in% "Andaman & Nicobar")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#f1c18e",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#f1c18e",
                               "#f98477",
                               "gray65",
                               "#ea4f88",
                               "#b1339e",
                               "#692a99"))+
  scale_color_manual(values = c("#f1c18e",
                                "#f98477",
                                "gray65",
                                "#ea4f88",
                                "#b1339e",
                                "#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

an.pca.density.fig <- ggMarginal(an.pca.fig,
                                 type = "density",
                                 groupColour = TRUE, 
                                 groupFill = TRUE)

# continental islands in orinma
conisorinma.pca.fig <- ggplot(orinma.pca.values |> 
                                filter(subregion %in% "Continental islands")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.2,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#f98477",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#f98477",
                               "gray65",
                               "#ea4f88",
                               "#b1339e",
                               "#692a99"))+
  scale_color_manual(values = c("#f98477",
                                "gray65",
                                "#ea4f88",
                                "#b1339e",
                                "#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

conisorinma.pca.density.fig <- ggMarginal(conisorinma.pca.fig,
                                          type = "density",
                                          groupColour = TRUE, 
                                          groupFill = TRUE)

#mainland orinma
mainlorinma.pca.fig <- ggplot(orinma.pca.values |> 
                                filter(subregion %in% "Mainland")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "gray95",
                 alpha = 0.4)+
  scale_fill_manual(values = c("gray65",
                               "#ea4f88",
                               "#b1339e",
                               "#692a99"))+
  scale_color_manual(values = c("gray65",
                                "#ea4f88",
                                "#b1339e",
                                "#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

mainlorinma.pca.density.fig <- ggMarginal(mainlorinma.pca.fig,
                                          type = "density",
                                          groupColour = TRUE, 
                                          groupFill = TRUE)

# Philippines orinma
phili.pca.fig <- ggplot(orinma.pca.values |> 
                          filter(subregion %in% "Philippines")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#ea4f88",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#ea4f88",
                               "#b1339e",
                               "#692a99"))+
  scale_color_manual(values = c("#ea4f88",
                                "#b1339e",
                                "#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

phili.pca.density.fig <- ggMarginal(phili.pca.fig,
                                    type = "density",
                                    groupColour = TRUE, 
                                    groupFill = TRUE)


# Sunda orinma
sunda.pca.fig <- ggplot(orinma.pca.values |> 
                          filter(subregion %in% "Sunda islands")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#b1339e",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#b1339e",
                               "#692a99"))+
  scale_color_manual(values = c("#b1339e",
                                "#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

sunda.pca.density.fig <- ggMarginal(sunda.pca.fig,
                                    type = "density",
                                    groupColour = TRUE, 
                                    groupFill = TRUE)

# Wallacea orinma
wallace.pca.fig <- ggplot(orinma.pca.values |> 
                            filter(subregion %in% "Wallacea")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#692a99",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#692a99"))+
  scale_color_manual(values = c("#692a99"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (29.47%)",
       y = "PC2 (24.58%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

wallace.pca.density.fig <- ggMarginal(wallace.pca.fig,
                                      type = "density",
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
  left_join(spp_cell_df_func |> 
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

#figure

global.papuan.pca.fig <- ggplot(papuan.pca.values) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 fill = subregion,
                 color = subregion), 
             alpha = 0.25,
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
  scale_x_continuous(limits = c(-10.8,12.3))+
  scale_y_continuous(limits = c(-10.9,7))+
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
       title = expression(bold("(f)")),
       fill = "Archipelago \n(subregion)",
       color = "Archipelago \n(subregion)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

global.papuan.pca.density.fig <- ggMarginal(global.papuan.pca.fig,
                                            type = "density",
                                            size = 15,
                                            groupColour = TRUE, 
                                            groupFill = TRUE)

### Papuan - per group ####

# Bismarcks
bismarck.pca.fig <- ggplot(papuan.pca.values |> 
                             filter(subregion %in% "Bismarcks")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#d20231",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#d20231",
                               "gray65",
                               "#ffd033",
                               "#62f476",
                               "#218be7"))+
  scale_color_manual(values = c("#d20231",
                                "gray65",
                                "#ffd033",
                                "#62f476",
                                "#218be7"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (37.26%)",
       y = "PC2 (19.78%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

bismarck.pca.density.fig <- ggMarginal(bismarck.pca.fig,
                                       type = "density",
                                       groupColour = TRUE, 
                                       groupFill = TRUE)

# mainland from Australia
# Bismarcks
mainlpapuan.pca.fig <- ggplot(papuan.pca.values |> 
                                filter(subregion %in% "Mainland")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "gray95",
                 alpha = 0.4)+
  scale_fill_manual(values = c("gray65",
                               "#ffd033",
                               "#62f476",
                               "#218be7"))+
  scale_color_manual(values = c("gray65",
                                "#ffd033",
                                "#62f476",
                                "#218be7"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (37.26%)",
       y = "PC2 (19.78%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

mainlpapuan.pca.density.fig <- ggMarginal(mainlpapuan.pca.fig,
                                          type = "density",
                                          groupColour = TRUE, 
                                          groupFill = TRUE)

# Papua
papua.pca.fig <- ggplot(papuan.pca.values |> 
                          filter(subregion %in% "Papua")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#ffd033",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#ffd033",
                               "#62f476",
                               "#218be7"))+
  scale_color_manual(values = c("#ffd033",
                                "#62f476",
                                "#218be7"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (37.26%)",
       y = "PC2 (19.78%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

papua.pca.density.fig <- ggMarginal(papua.pca.fig,
                                    type = "density",
                                    groupColour = TRUE, 
                                    groupFill = TRUE)

# Solomon Archipelago
solomon.pca.fig <- ggplot(papuan.pca.values |> 
                            filter(subregion %in% "Solomons")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#62f476",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#62f476",
                               "#218be7"))+
  scale_color_manual(values = c("#62f476",
                                "#218be7"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (37.26%)",
       y = "PC2 (19.78%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

solomon.pca.density.fig <- ggMarginal(solomon.pca.fig,
                                      type = "density",
                                      groupColour = TRUE, 
                                      groupFill = TRUE)

# Vanuatu
vanuatu.pca.fig <- ggplot(papuan.pca.values |> 
                            filter(subregion %in% "Vanuatu")) +
  geom_point(aes(x = PC1, 
                 y = PC2, 
                 color = subregion,
                 fill = subregion), 
             alpha = 0.4,
             shape = 21) +
  geom_density2d(aes(x = PC1, 
                     y = PC2),
                 color = "#218be7",
                 alpha = 0.4)+
  scale_fill_manual(values = c("#218be7"))+
  scale_color_manual(values = c("#218be7"))+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(-8,8))+
  coord_fixed()+
  labs(x = "PC1 (37.26%)",
       y = "PC2 (19.78%)")+
  theme(legend.position = "none",
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

vanuatu.pca.density.fig <- ggMarginal(vanuatu.pca.fig,
                                      type = "density",
                                      groupColour = TRUE, 
                                      groupFill = TRUE)

# Figure of global PCAs with density ####

#Combine

grid.arrange(global.carib.pca.density.fig,
             global.orinma.pca.density.fig, 
             global.papuan.pca.density.fig, 
             ncol = 3)


# woohoo


# Species approach ####

Caribbean_funbiog <- spp_cell_df_func |>
  filter(Meta.Archipelago %in% "Neotropical") |>
  dplyr::select(phylo_name, scientific_name, cell, mu_count, 
                subregion, latitude, longitude, 
                area_km2, 
                # dispersal morphology
                bod.mas, hwi,
                # foraging morphology
                tar.len, win.len,
                # dietary morphology
                bea.len, bea.wid, bea.dep,
                # ecological niche
                bod.len, hab.bre, ran.siz, end.ins,
                # behavioral/reproductive niche
                vertical,
                nes.str.bre, nes.loc.bre, 
                P.nes.str.cat, P.nes.loc.cat)

# We have to impute values

# traits per species (1319 spp)
Caribbean_spp_traits <- unique(Caribbean_funbiog[,c(1,9:24)])
rownames(Caribbean_spp_traits) <- Caribbean_spp_traits$phylo_name

str(Caribbean_spp_traits)

# confirm categorical variables as factor
Caribbean_spp_traits <- Caribbean_spp_traits |>
  mutate(P.nes.loc.cat = factor(P.nes.loc.cat),
         P.nes.str.cat = factor(P.nes.str.cat))

table(Caribbean_spp_traits$P.nes.loc.cat)

# scale (z-score) the numeric values
Caribbean_spp_traits_scl <- as.data.frame(scale(Caribbean_spp_traits[,c(2:15)]))

#add categorical variables
Caribbean_spp_traits_scl <- cbind(Caribbean_spp_traits[,16:17],Caribbean_spp_traits_scl)

# and we can prune the phylogenetic tree for the two regions
phylo.Caribb <- ape::drop.tip(birds_phylo, 
                              birds_phylo$tip.label[-match(Caribbean_spp_traits$phylo_name,
                                                           birds_phylo$tip.label)])
plotTree(phylo.Caribb, ftype = "i",
         fsize = 0.2, lwd = 1, type = "fan")

length(phylo.Caribb$tip.label)

#impute
Caribbean_spp_traits_imputed <- impute(traits = Caribbean_spp_traits_scl,
                                       phylo = phylo.Caribb)

#compare the effect of imputation
summary(Caribbean_spp_traits_scl)
summary(Caribbean_spp_traits_imputed$imputed)

imputed.traits <- Caribbean_spp_traits_imputed$imputed

funspaceDim(imputed.traits[,3:16]) # 4 dimensions

# PCA
pca.traits <- princomp(imputed.traits[,3:16], cor = T)

# building the functional trait space
trait_space_Caribbean_global_PC1_PC2 <- funspace(x = pca.traits, PCs = c(1,2))
trait_space_Caribbean_global_PC2_PC3 <- funspace(x = pca.traits, PCs = c(2,3))

par(mfrow = c(1,2))
plot(x = trait_space_Caribbean_global_PC1_PC2, type = "global",
     quant.plot = T, arrows = T, arrows.length = 0.9)
plot(x = trait_space_Caribbean_global_PC2_PC3, type = "global",
     quant.plot = T, arrows = T, arrows.length = 0.9)

# Grouping
Caribbean_subregion <- Caribbean_funbiog |>
  group_by(phylo_name, subregion) |>
  count() |> 
  group_by(phylo_name) |>
  summarise(subregion = paste(unique(subregion), collapse = ", "),
            n.subregions = str_count(paste(unique(subregion),
                                           collapse = ", "), ", ") +1) |>
  mutate(subregion = gsub("Caribbean - ", "", subregion),
         subregion = ifelse(n.subregions >= 3, "widespread (≥3)",
                            subregion))

trait_space_Caribbean_subregion_PC1_PC2 <- funspace(pca.traits,
                                                    PCs = c(1,2), 
                                                    group.vec = Caribbean_subregion$subregion)


plot(x = trait_space_Caribbean_subregion_PC1_PC2, type = "groups",
     quant.plot = T, globalContour = T, pnt = T)

summary(trait_space_Caribbean_subregion_PC1_PC2)