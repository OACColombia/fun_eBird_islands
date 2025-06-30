# Organizing eBird data for stationary distribution population models

#### Packages ####
library(tidyverse)

#### Data from eBird ####

Caribbean_SS <- readRDS("Completeness_data_Islands/Caribbean_SS.rds") |>
  as.data.frame()

Caribbean_OUSS_Cells <- Caribbean_SS |>
  group_by(scientific_name, cell, Time.t) |> 
  count() |>
  group_by(cell, scientific_name) |>
  count() |>
  filter(n > 5) |>
  as.data.frame() |>
  rename(n.time = n)

Caribbean_OUSS <- Caribbean_SS |>
  left_join(Caribbean_OUSS_Cells) |>
  filter(!is.na(n.time)) 

Caribbean_OUSS$cell <- as.character(Caribbean_OUSS$cell)

saveRDS(Caribbean_OUSS, "Completeness_data_Islands/Caribbean_preOUSS.rds")

spp <- unique(Caribbean_OUSS$scientific_name) # 1752
sites <- unique(Caribbean_OUSS$cell) # 2991

# and now with the Indo Pacific region

IndoPacific_SS <- readRDS("Completeness_data_Islands/Caribbean_SS.rds") |>
  as.data.frame()

IndoPacific_OUSS_Cells <- IndoPacific_SS |>
  group_by(scientific_name, cell, Time.t) |> 
  count() |>
  group_by(cell, scientific_name) |>
  count() |>
  filter(n > 5) |>
  as.data.frame() |>
  rename(n.time = n)

IndoPacific_OUSS <- IndoPacific_SS |>
  left_join(IndoPacific_OUSS_Cells) |>
  filter(!is.na(n.time)) 

IndoPacific_OUSS$cell <- as.character(IndoPacific_OUSS$cell)

saveRDS(IndoPacific_OUSS, "Completeness_data_Islands/IndoPacific_preOUSS.rds")

sppIP <- unique(IndoPacific_OUSS$scientific_name) # 1752
sitesIP <- unique(IndoPacific_OUSS$cell) # 2991

# End of this code