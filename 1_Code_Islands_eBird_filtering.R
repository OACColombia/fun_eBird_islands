########################################
###~ eBird data - Island filtering ~####
########################################

#Packages
library(auk); #eBird data filters
library(tidyverse); #manipulate data
library(dggridR) #Spatiotemporal Subsampling
library(sf) # Simple Features, to encode spatial vector data

#~ eBird data ~####
#Data downloaded from eBird, up to May 2024, saved in the folder of the HPG

#Define the range
CaribbeanRange <- c(-91,8,-59,27)
IndoMalayRange <- c(90,-20,180,20)

#select the columns to extract (based on x$col_idx$name)
colsE <- c("observer_id", "sampling_event_identifier",
           "group identifier",
           "common_name", "scientific_name",
           "observation_count",
           "country", "state_code", "locality_id", "latitude", "longitude",
           "protocol_type", "all_species_reported",
           "observation_date",
           "time_observations_started",
           "duration_minutes", "effort_distance_km",
           "number_observers")

#For Caribbean
f_ebdCa <- "ebd_IslandsCa.txt" #Temporal file to save the filtering eBird data (records)

ebd_filtCarib<-auk_ebd("ebd_relMay-2024.txt") %>%
  auk_bbox(CaribbeanRange) %>% #W, S, E, N
  auk_year(c(2015:2024)) %>%
  auk_protocol(c("Traveling", "Stationary")) %>%
  auk_distance(distance = c(0,5)) %>%
  auk_duration(duration = c(0,300))%>%
  auk_complete() %>%
  auk_filter(f_ebdCa, overwrite=T, keep = colsE)

f_ebdIM <- "ebd_IslandsIM.txt" #Temporal file to save the filtering eBird data (records)

ebd_filtIndoP<-auk_ebd("ebd_relMay-2024.txt") %>%
  auk_bbox(IndoMalayRange) %>% #W, S, E, N
  auk_year(c(2015:2024)) %>%
  auk_protocol(c("Traveling", "Stationary")) %>%
  auk_distance(distance = c(0,5)) %>%
  auk_duration(duration = c(0,300))%>%
  auk_complete() %>%
  auk_filter(f_ebdIM, overwrite=T, keep = colsE)

# and with read_ebd I apply another filter to do not repeat records from groups
ebd_Carib <- read_ebd(f_ebdCa)
ebd_IndoM <- read_ebd(f_ebdIM)

# Function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

#Caribbean data
Carib_eff <- ebd_Carib |>
  mutate(
    # I don't want count in X, to convert to NA
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type == "Stationary",
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    hour_sampling = round(time_observations_started, 0),
    # split date into year, month, week, and day of year
    year = year(observation_date),
    month = month(observation_date),
    week = week(observation_date),
    day_of_year = yday(observation_date),
    # temporal aggregation - month
    Time.t = case_when(year == 2015 ~ month,
                       year > 2015 ~ month+(12*(year-2015)))) |>
  filter(number_observers <= 10,         #Only list with less than 10 observers
         effort_distance_km <= 5,        #be sure of distance effort
         duration_minutes %in% (0:300), #be sure of duration effort
         !is.na(observation_count))      #only records with abundance estimation

Checklists_with_more_SppC <- Carib_eff %>%
  dplyr::select(checklist_id, scientific_name) %>%
  group_by(checklist_id) %>%
  summarise(N_species = n()) %>%
  filter(N_species >= 10) # minimum number of species per sampling

Caribbean_cooc <- Carib_eff[Carib_eff$checklist_id %in% Checklists_with_more_SppC$checklist_id, ]

#And now with the Oriental-Indo_Malayan-Papuan-Melanesian
IndoM_eff <- ebd_IndoM |>
  mutate(
    # I don't want here count in X, to convert to NA
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type == "Stationary",
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    hour_sampling = round(time_observations_started, 0),
    # split date into year, month, week, and day of year
    year = year(observation_date),
    month = month(observation_date),
    week = week(observation_date),
    day_of_year = yday(observation_date),
    # temporal aggregation - month
    Time.t = case_when(year == 2015 ~ month,
                       year > 2015 ~ month+(12*(year-2015)))) |>
  filter(number_observers <= 10,         #Only list with less than 10 observers
         effort_distance_km <= 5,        #be sure of distance effort
         duration_minutes %in% (0:300), #be sure of duration effort
         !is.na(observation_count))      #remove records without abundance estimation

Checklists_with_more_SppIP <- IndoM_eff %>%
  dplyr::select(checklist_id, scientific_name) %>%
  group_by(checklist_id) %>%
  summarise(N_species = n()) %>%
  filter(N_species >= 10)

IndoM_eff_cooc <- IndoM_eff[IndoM_eff$checklist_id %in% Checklists_with_more_SppIP$checklist_id, ]

##Spatial Grids ####

# Spatiotemporal Subsampling in diameters of ~11 km (area of ~100 km^2)
  # spacing 1 correspond to Characteristic Length Scale, or diameter of spherical cell
dggs_island <- dgconstruct(spacing = 11) 

#add a new variable that identify cell
Caribbean_cooc_cells <- Caribbean_cooc |>
  mutate(cell = dgGEO_to_SEQNUM(dggs_island, #id for cells
                                longitude, latitude)$seqnum) |>
  group_by(cell, scientific_name, Time.t) |>
  # extract the minimum number of individual detected per month per cell of each species
  mutate(Observed.y = round(max(observation_count))) |> 
  ungroup()

# filter cells with more than 10 localities within
cells_10_Ca <- Caribbean_cooc_cells |>
  group_by(cell, locality_id) |>
  count() |>
  group_by(cell) |>
  count() |>
  ungroup() |>
  filter(n >= 10)

Caribbean_SS = Caribbean_cooc_cells[Caribbean_cooc_cells$cell %in% cells_10_Ca$cell, ]

saveRDS(Caribbean_SS, "Completeness_data_Islands/Caribbean_SS.rds")

# now with the IndoPacific

IndoM_eff_cooc_cells <- IndoM_eff_cooc %>%
  filter(observation_count < 7000) %>% #there is an outlying count for _Eurylaimus ochromalus_ https://ebird.org/checklist/S112199250
  mutate(cell = dgGEO_to_SEQNUM(dggs_island, #id for cells
                                longitude, latitude)$seqnum) %>%
  group_by(cell, scientific_name, Time.t) |>
  # extract the minimum number of individual detected per month per cell of each species
  mutate(Observed.y = round(max(observation_count))) |> 
  ungroup()

# filter cells with more than 10 localities within
cells_10_IP <- IndoM_eff_cooc_cells |>
  group_by(cell, locality_id) |>
  count() |>
  group_by(cell) |>
  count() |>
  ungroup() |>
  filter(n >= 10)

IndoPacific_SS = IndoM_eff_cooc_cells[IndoM_eff_cooc_cells$cell %in% cells_10_IP$cell, ]

saveRDS(IndoPacific_SS, "Completeness_data_Islands/IndoPacific_SS.rds")

# Get the number of records in each cell
RichnessCellCaribbean   <- Caribbean_SS %>%
  group_by(cell, scientific_name) %>%
  summarise(count=n()) %>%
  group_by(cell) %>%
  summarise(SpRichness=n())

# Get the grid cell boundaries for cells which had quakes
gridCaribbean <- dgcellstogrid(dggs_island,RichnessCellCaribbean$cell)

# Update the grid cells' properties to include the number of lists in each cell
gridCaribbean <- merge(gridCaribbean,RichnessCellCaribbean,by.x="seqnum",by.y="cell")

# Handle cells that cross 180 degrees
wrapped_gridCaribbean = st_wrap_dateline(gridCaribbean,
                                         options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

saveRDS(wrapped_gridCaribbean, "Completeness_data_Islands/wrapped_gridCaribbean.rds")

# Get the number of records in each cell
RichnessCellIndoPacific   <- IndoPacific_SS %>%
  group_by(cell, scientific_name) %>%
  summarise(count=n()) %>%
  group_by(cell) %>%
  summarise(SpRichness=n())

#Get the grid cell boundaries for cells which had quakes
gridIndoPacific <- dgcellstogrid(dggs_island,RichnessCellIndoPacific$cell)

#Update the grid cells' properties to include the number of lists in each cell
gridIndoPacific <- merge(gridIndoPacific,RichnessCellIndoPacific,by.x="seqnum",by.y="cell")

# Handle cells that cross 180 degrees (and some weird in the map)
wrapped_gridIndoPacific = st_wrap_dateline(gridIndoPacific,
                                           options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

saveRDS(wrapped_gridIndoPacific, "Completeness_data_Islands/wrapped_gridIndoPacific.rds")

#End of this code
