# Functional traits and trait space

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

# Data from eBird records  ####
  # summarise sampling unit and convert to sf spatial object

Caribbean_Spp_sf <- readRDS("Completeness_data_Islands/Caribbean_OUSS_with_parms.rds") |>
  # C_ouss |>
  group_by(scientific_name, cell) |>
  summarise(latitude = mean(latitude, na.rm = TRUE),
            longitude = mean(longitude, na.rm = TRUE),
            mu_hat = mean(mu_hat, na.rm = TRUE),
            omega_hat = mean(omega_hat, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(Meta.Archipelago = "Neotropical") |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326)

IndoPacific_Spp_sf <- readRDS("Completeness_data_Islands/IndoPacific_OUSS_with_parms.rds") |>
  # IP_ouss |>
  group_by(scientific_name, cell) |>
  summarise(latitude = mean(latitude, na.rm = TRUE),
            longitude = mean(longitude, na.rm = TRUE),
            mu_hat = mean(mu_hat, na.rm = TRUE),
            omega_hat = mean(omega_hat, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(Meta.Archipelago = "Indo.Pacific") |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326)

islands_spp_cells <- rbind(Caribbean_Spp_sf, IndoPacific_Spp_sf)

# recover latitude and longitude
islands_spp_cells$longitude <- st_coordinates(islands_spp_cells)[, 1]
islands_spp_cells$latitude <- st_coordinates(islands_spp_cells)[, 2]

# Polygons of the meta-archipelagos ####
# Mainland in america polygon
mainland_sf <- st_read("Coastline_Islands/Mainland_coastline.shp") |>
  st_transform(crs = 4326) |>
  mutate(subregion = "Mainland")

# Caribbean
# Continental Islands polygon
cont_isln_Car <- st_read("Coastline_Islands/ContinentalIslands.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Caribbean - Continental islands")
# Greater Antilles polygon
gr_antilles <- st_read("Coastline_Islands/GreaterAntilles.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Caribbean - Greater Antilles")
# Kalinago or Lesser Antilles polygon
kalinago <- st_read("Coastline_Islands/Kalinago.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Caribbean - Kalinago")
# Lucayan or Bahamas polygon
lucayan <- st_read("Coastline_Islands/Lucayan.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Caribbean - Lucayan")

# Oriental-Indo-Malayan
# Andaman & Nicobar polygon
andaman <- st_read("Coastline_Islands/AndamanNicobar.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Oriental-Indo-Malayan - Andaman & Nicobar")
# Continental Islands - Oriental polygon
cont_isln_OrInd <- st_read("Coastline_Islands/ContinentalIslandsOriental.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Oriental-Indo-Malayan - Continental islands")
# Philippines polygon
philippines <- st_read("Coastline_Islands/Philippines.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Oriental-Indo-Malayan - Philippines")
# Sunda Islands polygon
sunda <- st_read("Coastline_Islands/SundaIslands.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Oriental-Indo-Malayan - Sunda islands")
# Wallacea polygon
wallacea <- st_read("Coastline_Islands/Wallacea.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Oriental-Indo-Malayan - Wallacea")

# Papuan-Melanesian
# Bismarck polygon
bismarcks <- st_read("Coastline_Islands/Biskmark.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Papuan-Melanesian - Bismarcks")
# New Guinea and Australian continental islands polygon
papua <- st_read("Coastline_Islands/NewGuinea_AU_ContIslands.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Papuan-Melanesian - Papua")
# Solomons polygon
solomons <- st_read("Coastline_Islands/Solomon.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Papuan-Melanesian - Solomons")
# Vanuatu polygon
vanuatu <- st_read("Coastline_Islands/Vanuatu.shp") |>
  st_transform(crs = 4326) |>
  st_make_valid() |>
  mutate(subregion = "Papuan-Melanesian - Vanuatu")

# Combine the sf objects for the islands
islands_sf <- rbind(cont_isln_Car, gr_antilles, kalinago, lucayan, 
                     cont_isln_OrInd, andaman, philippines, sunda, wallacea,  
                     papua, bismarcks, solomons, vanuatu)

# Spatial join to filter points withIN the mainland ####
mainland_cells <- st_join(islands_spp_cells, 
                           mainland_sf, 
                           join = st_within)
# 404047 to 202321 cell-species rows 

ggplot()+
  geom_sf(data = mainland_cells)+
  geom_point(data = mainland_cells |> filter(is.na(subregion)), 
             aes(x = longitude, y = latitude),
             color = "red", 
             size = 0.05)+
  geom_point(data = mainland_cells |> filter(!is.na(subregion)), 
             aes(x = longitude, y = latitude),
             color = "blue", 
             size = 0.1)+
  coord_sf(ylim=c(7,28),
           xlim=c(-58,-92))+
  theme_classic()

#Filter only points within the mainland
mainland_cells <- mainland_cells[!is.na(mainland_cells$subregion), ]
# 160839 observations

#recover latitude and longitude
mainland_cells$longitude <- st_coordinates(mainland_cells)[, 1]
mainland_cells$latitude <- st_coordinates(mainland_cells)[, 2]

#bring back the df structure
spp_cell_mainland <- mainland_cells |>
  st_drop_geometry() |>
  as.data.frame()

# Spatial join to filter points withOUT the mainland ####
island_cells <- st_join(islands_spp_cells, 
                          islands_sf, 
                          join = st_within)
# 202321 cell-species rows 

ggplot()+
  geom_sf(data = island_cells)+
  geom_point(data = island_cells |> filter(is.na(subregion)), 
             aes(x = longitude, y = latitude),
             color = "red", 
             size = 0.05)+
  geom_point(data = island_cells |> filter(!is.na(subregion)), 
             aes(x = longitude, y = latitude),
             color = "blue", 
             size = 0.1)+
  coord_sf(ylim=c(7,28),
           xlim=c(-58,-92))+
  theme_classic()

#Filter only points withOUT the mainland
island_cells <- island_cells[!is.na(island_cells$subregion), ]
# 32776 (different than expected from Mainland! where are the 8706 missed?)

#recover latitude and longitude
island_cells$longitude <- st_coordinates(island_cells)[, 1]
island_cells$latitude <- st_coordinates(island_cells)[, 2]

#bring back the df structure
spp_cell_islands <- island_cells |>
  st_drop_geometry() |>
  as.data.frame()

# Combine df with subregion #### 
spp_cell_df <- spp_cell_mainland |> 
  full_join(spp_cell_islands)

# some cells at edge or out of the subregions were removed 
  # (e.g. Fiji and Micronesia)

# Save data for PCAs figures ####

saveRDS(spp_cell_df, "Completeness_data_Islands/Species_in_cells_MetaArchipelagos.rds")

# Functional traits databases ####

## AVONET raw data ####
avonet <- read_csv("FunctionalTraits/AVONET_Raw_Data.csv") |>
  rename(Hand.Wing.Index = `Hand-wing.Index`) |>
  as.data.frame()

# the raw data of AVONET does not have body mass nor trophic niche
avonet_eBird <- read_csv("FunctionalTraits/AVONET2_eBird.csv") |> 
  dplyr::select(Species2, Mass, Trophic.Niche) |>
  rename(Species2_eBird = Species2) |>
  # seven species had NA in Trophic.Niche. We adjusted based on BoW
  mutate(Trophic.Niche = ifelse(Species2_eBird == "Vanellus macropterus","Omnivore",
                         ifelse(Species2_eBird == "Ophrysia superciliosa","Granivore",
                         ifelse(Species2_eBird == "Otus alius","Invertivore",
                         ifelse(Species2_eBird == "Otus beccarii","Invertivore",
                         ifelse(Species2_eBird == "Otus hartlaubi","Invertivore",
                         ifelse(Species2_eBird == "Otus mirus","Invertivore",
                         ifelse(Species2_eBird == "Otus siaoensis","Invertivore",
                         Trophic.Niche)))))))) |>
  as.data.frame()

length(table(avonet$Species2_eBird))
nrow(avonet_eBird)

# different number of species! - lets unify

# unify AVONET means
avonet_means <- avonet |> 
  group_by(Species2_eBird) |>
  summarise(nInd = n(),
            hwi = mean(Hand.Wing.Index, na.rm = T),
            tarsus.l = mean(Tarsus.Length, na.rm = T),
            wing.l = mean(Wing.Length, na.rm = T),
            beak.l = mean(Beak.Length_Culmen, na.rm = T),
            beak.w = mean(Beak.Width, na.rm = T),
            beak.d = mean(Beak.Depth, na.rm = T)) |>
  drop_na(Species2_eBird) |>
  as.data.frame()

# joint to the summary of AVONET2_eBird - Mass and Trophic.Niche
avonet_traits <- avonet_means |>
  left_join(avonet_eBird)  |>
  as.data.frame()

avonet_traits <- avonet_traits |>
  mutate(predatory = ifelse(Trophic.Niche == "Vertivore",1,0),
         `diet generalist` = ifelse(Trophic.Niche == "Omnivore", 1, 0)) |>
  # to be combine with tetrapod and nest traits datasets
  left_join(avonet |> 
              dplyr::select(Species1_BirdLife, 
                            Species2_eBird,
                            Species3_BirdTree)) |>
  unique()

# Tetrapod Traits - Aves (join with AVONET by Scientific.Name == Species3_BirdTree)
tetrapod <- read_csv("FunctionalTraits/TetrapodTraits_1.0.0.csv") |>
  filter(Class %in% "Aves")

# Nest traits for the world's birds - I have replaced "u" (of uncertain) for 0.5
nest.traits <- read_csv("FunctionalTraits/NestTraits_Dataset-S2.csv") |> 
  # convert to a longer format, extracting the "nest structure" and removing no
  pivot_longer(cols = starts_with("Str-"),
               names_to = "nest.str",
               values_to = "value.str") |>
  filter(value.str > 0) |> 
  # convert to a longer format, extracting the "nest structure"
  pivot_longer(cols = starts_with("Loc-"),
               names_to = "nest.loc",
               values_to = "value.loc") |> 
  filter(value.loc > 0) |>
  # remove the initial categories 
  mutate(nest.str = gsub("Str-", "", nest.str),
         nest.loc = gsub("Loc-", "", nest.loc)) |>
  # combine the categories
    # and count the number of major structure and location per species
  group_by(`Species scientific name`, HWI, `Body mass (log)`) |>
  summarise(nest.str.cat = paste(unique(nest.str[value.str > 0]), collapse = ", "),
            nest.loc.cat = paste(unique(nest.loc[value.loc > 0]), collapse = ", "),
            nest.str.bre = str_count(paste(unique(nest.str.cat),
                                         collapse = ", "), ", ") +1,
            nest.loc.bre = str_count(paste(unique(nest.loc.cat),
                                         collapse = ", "), ", ") +1) |>
  ungroup() |>
  mutate(Body.mass = 10^(`Body mass (log)`))

# Combine functional traits - Global functional traits ####
functraits <- avonet_traits |> 
  left_join(tetrapod, 
             join_by(Species3_BirdTree == Scientific.Name)) |> 
  left_join(nest.traits, 
             join_by(Species1_BirdLife == `Species scientific name`)) |>
  # use only the eBird (Clements) taxonomy
  group_by(Species2_eBird) |>
  summarise(
    n_row = n(),
    # Dispersal morphology
    bod.mas = mean(c(Mass, BodyMass_g, Body.mass), na.rm = T),   #avonet-tetrapod-nest
    sd.bod.mas = sd(c(Mass, BodyMass_g, Body.mass), na.rm = T),   #avonet-tetrapod-nest
    hwi = mean(c(hwi, HWI), na.rm = T),
    sd.hwi = sd(c(hwi, HWI), na.rm = T),
    # Foraging morphology
    tar.len = mean(tarsus.l, na.rm = T),
    win.len = mean(wing.l, na.rm = T),
    # Dietary morphology
    bea.len = mean(beak.l, na.rm = T),
    bea.wid = mean(beak.w, na.rm = T),
    bea.dep = mean(beak.d, na.rm = T),
    # Ecological niche
    bod.len = mean(BodyLength_mm, na.rm = T),
    hab.bre = mean(MajorHabitatSum, na.rm = T),
    ran.siz = mean(RangeSize, na.rm = T),
    end.ins = mean(Insularity, na.rm = T),
    # Foraging niche
    predatory = mean(predatory, na.rm = T),
    generalist = mean(`diet generalist`, na.rm = T),
    # Behavioral niche
    vertical = mean(Verticality, na.rm = T),
    nocturnal = mean(Nocturnality, na.rm = T),
    # Reproductive niche
    nes.str.bre = mean(nest.str.bre, na.rm = T),
    nes.loc.bre = mean(nest.loc.bre, na.rm = T),
    nes.str.cat = paste(unique(nest.str.cat), collapse = ", "),
    nes.loc.cat = paste(unique(nest.loc.cat), collapse = ", ")) |> 
  unique() |> 
  as.data.frame() |>
  filter(!is.na(Species2_eBird))

length(unique(functraits$Species2_eBird))
# 10521 species

summary(functraits) #traits have NAs - Impute from phylogeny - 

# Phylogenetic tree
birds_phylo <- ape::read.nexus("summary_dated_clements.nex")
length(birds_phylo$tip.label)

functraits$phylo_name <- gsub(" ", "_", functraits$Species2_eBird)

rownames(functraits) <- functraits$phylo_name

imputed_traits <- impute(functraits[,c(3,5,7:21)], phylo = birds_phylo)

# Save the imputed traits
functraits_global <- imputed_traits$imputed

functraits_global$Species2_eBird <- functraits$Species2_eBird

functraits_global <- functraits_global |> 
  # recover the SD of different mass measures and hwi estimation
  left_join(functraits) |>
  # correct name used in recent phylogeny following Clements
  mutate(phylo_name = sub(" ", "_", Species2_eBird))

saveRDS(functraits_global, "FunctionalTraits/birds_functraits_imputed.rds")

# Functional traits for the species in the islands cells - CMW ####

## bring back the data - and leverage on sister species for recent splits ####

birds_phylo <- ape::read.nexus("summary_dated_clements.nex")
functraits_global <- readRDS("FunctionalTraits/birds_functraits_imputed.rds")
spp_cell_df <- readRDS("Completeness_data_Islands/Species_in_cells_MetaArchipelagos.rds") 

list_species <- spp_cell_df |>
  group_by(scientific_name) |>
  count()

functraits_islands <- functraits_global |>
  filter(Species2_eBird %in% list_species$scientific_name)

# some (173) species in the dataset of islands do not have functraits...
no_traits_spp <- list_species |> 
  filter(!scientific_name %in% functraits_global$Species2_eBird)

# we can bring the information from available sister species or correct taxonomic
  # generate a "Sister" column for those species with data to those without data
no_traits_spp$Sister <- c("Acridotheres burmannicus", # "Acridotheres leucocephalus"
                          "Psaltria exilis", # Aegithalos exilis 
                          "Aerodramus infuscatus", # Aerodramus sororum
                          "Sericornis arfakianus", # Aethomyias arfakianus
                          "Sericornis papuensis", # Aethomyias papuensis
                          "Sericornis perspicillatus", # Aethomyias perspicillatus
                          "Sericornis rufescens", # Aethomyias rufescens
                          "Sericornis spilodera", # Aethomyias spilodera
                          "Ailuroedus crassirostris", # Ailuroedus arfakianus
                          "Ailuroedus crassirostris", # Ailuroedus maculosus 
                          "Alcedo euryzona", # Alcedo peninsulae
                          "Alcippe morrisonia", # Alcippe fratercula
                          "Alcippe morrisonia", # Alcippe hueti
                          "Alcippe morrisonia", # Alcippe striatus
                          "Charadrius alexandrinus", # Anarhynchus alexandrinus
                          "Charadrius mongolus", # Anarhynchus atrifrons
                          "Charadrius collaris", # Anarhynchus collaris
                          "Charadrius dealbatus", # Anarhynchus dealbatus
                          "Charadrius javanicus", # Anarhynchus javanicus
                          "Charadrius leschenaultii", # Anarhynchus leschenaultii
                          "Charadrius mongolus", # Anarhynchus mongolus
                          "Charadrius nivosus", # Anarhynchus nivosus
                          "Charadrius peronii", # Anarhynchus peronii
                          "Charadrius ruficapillus", # Anarhynchus ruficapillus
                          "Charadrius veredus", # Anarhynchus veredus
                          "Charadrius wilsonia", # Anarhynchus wilsonia
                          "Anthracothorax dominicus", # Anthracothorax aurulentus
                          "Anthus novaeseelandiae", # Anthus australis
                          "Anthus lutescens", # Anthus chii
                          "Antrostomus cubanensis", # Antrostomus ekmani
                          "Harpactes reinwardtii", # Apalharpactes reinwardtii
                          "Apus pacificus", # Apus cooki
                          "Ardea intermedia", # Ardea plumifera
                          "Automolus ochrolaemus", # Automolus cervinigularis
                          "Automolus subulatus", # Automolus virgatus
                          "Brachypteryx montana", # Brachypteryx erythrogyna
                          "Brachypteryx montana", # Brachypteryx poliogyna
                          "Bubulcus ibis", # Bubulcus coromandus
                          "Cecropis striolata", # Cecropis badia
                          "Elseyornis melanops", # Charadrius melanops
                          "Charmosyna papou", # Charmosyna stellae
                          "Chlorophonia musica", # Chlorophonia flavifrons
                          "Chlorophonia musica", # Chlorophonia sclateri
                          "Chloropsis cochinchinensis", # Chloropsis moluccensis
                          "Reinwardtipicus validus", # Chrysocolaptes validus
                          "Cinnyris jugularis", # Cinnyris aurora
                          "Cinnyris jugularis", # Cinnyris frenatus
                          "Cinnyris jugularis", # Cinnyris ornatus
                          "Coeligena torquata", # Coeligena conradii
                          "Coeligena bonapartei", # Coeligena conradii
                          "Collocalia esculenta", # Collocalia affinis
                          "Collocalia esculenta", # Collocalia isonata
                          "Collocalia esculenta", # Collocalia marginata
                          "Collocalia esculenta", # Collocalia natalis
                          "Collocalia esculenta", # Collocalia neglecta
                          "Collocalia esculenta", # Collocalia sumbawae
                          "Collocalia esculenta", # Collocalia uropygialis
                          "Colluricincla megarhyncha", # Colluricincla affinis
                          "Colluricincla megarhyncha", # Colluricincla fortis
                          "Colluricincla megarhyncha", # Colluricincla obscura
                          "Contopus cinereus", # Contopus bogotensis
                          "Corvus palmarum", # Corvus minutus
                          "Corvus enca", # Corvus pusillus
                          "Cuculus canorus", # Cuculus optatus
                          "Cyanoderma erythropterum", # Cyanoderma bicolor
                          "Cyornis whitei", # Cyornis montanus
                          "Milvago chimachima", # Daptrius chimachima
                          "Microeca papuana", # Devioeca papuana
                          "Dicaeum ignipectus", # Dicaeum cambodianum
                          "Dicaeum hirundinaceum", # Dicaeum keiense
                          "Dicaeum ignipectus", # Dicaeum luzoniense
                          "Dicrurus hottentottus", # Dicrurus palawanensis
                          "Dicrurus hottentottus", # Dicrurus striatus
                          "Cicinnurus magnificus", # Diphyllodes magnificus
                          "Cicinnurus respublica", # Diphyllodes respublica
                          "Eclectus roratus", # Eclectus polychloros
                          "Enicurus leschenaulti", # Enicurus borneensis
                          "Tregellasia capito", # Eopsaltria capito
                          "Tregellasia leucops", # Eopsaltria leucops
                          "Erythropitta macklotii", # Erythropitta habenichti
                          "Cyornis hoevelli", # Eumyias hoevelli
                          "Cyornis hyacinthinus", # Eumyias hyacinthinus
                          "Cyornis oscillans", # Eumyias oscillans
                          "Cyornis oscillans", # Eumyias stresemanni
                          "Falcunculus frontatus", # Falcunculus whitei
                          "Formicivora grisea", # Formicivora intermedia
                          "Furnarius leucopus", # Furnarius longirostris
                          "Pomatostomus isidorei", # Garritornis isidorei
                          "Gelochelidon nilotica", # Gelochelidon macrotarsa
                          "Gracupica contra", # Gracupica floweri
                          "Spodiornis rusticus", # Haplospiza rustica
                          "Heliangelus amethysticollis", # Heliangelus clarisse
                          "Heliangelus amethysticollis", # Heliangelus spencei
                          "Burhinus bistriatus", # Hesperoburhinus bistriatus
                          "Heteromyias albispecularis", # Heteromyias armiti
                          "Hypnelus ruficollis", # Hypnelus bicinctus
                          "Alophoixus longirostris", # Hypsipetes chloris
                          "Hypsipetes philippinus", # Hypsipetes guimarasensis
                          "Haliaeetus humilis", # Icthyophaga humilis
                          "Haliaeetus ichthyaetus", # Icthyophaga ichthyaetus
                          "Haliaeetus leucogaster", # Icthyophaga leucogaster
                          "Alophoixus finschii", # Iole finschii
                          "Irena puella", # Irena tweeddalii
                          "Microeca griseoceps", # Kempiella griseoceps
                          "Bubo sumatranus", # Ketupa sumatrana
                          "Lepidothrix coronata", # Lepidothrix velutina
                          "Locustella mandelli", # Locustella idonea
                          "Lophura ignita", # Lophura rufa
                          "Meiglyptes tristis", # Meiglyptes grammithorax
                          "Peneothello cryptoleuca", # Melanodryas cryptoleuca
                          "Peneothello cyanus", # Melanodryas cyanus
                          "Eopsaltria pulverulenta", # Melanodryas pulverulenta
                          "Peneothello sigillata", # Melanodryas sigillata
                          "Cracticus quoyi", # Melloria quoyi
                          "Melopyrrha nigra", # Melopyrrha taylori
                          "Accipiter superciliosus", # Microspizias superciliosus
                          "Brachypodius eutilotus", # Microtarsus eutilotus
                          "Brachypodius fuscoflavescens", # Microtarsus fuscoflavescens
                          "Brachypodius melanocephalos", # Microtarsus melanocephalos
                          "Brachypodius urostictus", # Microtarsus urostictus
                          "Mionectes olivaceus", # Mionectes galbinus
                          "Myiopagis caniceps", # Myiopagis parambae
                          "Polihierax insignis", # Neohierax insignis
                          "Sericornis citreogularis", # Neosericornis citreogularis
                          "Phaeomyias murina", # Nesotriccus incomtus
                          "Niltava vivida", # Niltava oatesi
                          "Ochthoeca fumicolor", # Ochthoeca superciliosa
                          "Oriolus cruentus", # Oriolus consanguineus
                          "Oriolus xanthonotus", # Oriolus consobrinus
                          "Paradisaea rudolphi", # Paradisornis rudolphi
                          "Pellorneum capistratum", # Pellorneum capistratoides
                          "Pellorneum capistratum", # Pellorneum nigrocapitatum
                          "Phacellodomus rufifrons", # Phacellodomus inornatus
                          "Pitangus lictor", # Philohydor lictor
                          "Phyllomyias burmeisteri", # Phyllomyias zeledoni
                          "Phylloscopus sarasinorum", # Phylloscopus nesophilus
                          "Phylloscopus trivirgatus", # Phylloscopus nigrorum
                          "Pitta sordida", # Pitta abbotti
                          "Pitta concinna", # Pitta elegans
                          "Pitta sordida", # Pitta novaeguineae
                          "Pitta sordida", # Pitta rosenbergii
                          "Poecilodryas albonotata", # Plesiodryas albonotata
                          "Phylloscartes ophthalmicus", # Pogonotriccus ophthalmicus
                          "Polioptila plumbea", # Polioptila albiventris
                          "Pomatorhinus montanus", # Pomatorhinus bornensis
                          "Porphyrio indicus", # Porphyrio melanotus
                          "Porphyrio indicus", # Porphyrio poliocephalus
                          "Porphyrio indicus", # Porphyrio pulverulentus
                          "Amblyornis newtoniana", # Prionodura newtoniana
                          "Psephotus chrysopterygius", # Psephotellus chrysopterygius
                          "Psephotus dissimilis", # Psephotellus dissimilis
                          "Psilopogon duvaucelii", # Psilopogon cyanotis
                          "Psittacara holochlorus", # Psittacara strenuus
                          "Pteruthius flaviscapis", # Pteruthius aeralatus
                          "Ptilinopus solomonensis", # Ptilinopus speciosus
                          "Pycnonotus blanfordi", # Pycnonotus conradi
                          "Pycnonotus finlaysoni", # Pycnonotus davisoni
                          "Pycnonotus flavescens", # Pycnonotus leucops
                          "Rhipidura dryas", # Rhipidura semicollaris
                          "Rhynchocyclus olivaceus", # Rhynchocyclus aequinoctialis
                          "Micropygia schomburgkii", # Rufirallus schomburgkii
                          "Saudareos ornatus", # Saudareos ornata
                          "Ochthoeca diadema", # Silvicultrix diadema
                          "Streptopelia chinensis", # Spilopelia chinensis
                          "Ciccaba nigrolineata", # Strix nigrolineata
                          "Ciccaba virgata", # Strix virgata
                          "Taenioptynx brodiei", # Taenioptynx sylvaticus
                          "Todiramphus chloris", # Todiramphus sacer
                          "Todiramphus chloris", # Todiramphus sordidus
                          "Tolmomyias assimilis", # Tolmomyias flavotectus
                          "Trochilus polytmus", # Trochilus scitulus
                          "Trogon rufus", # Trogon tenellus
                          "Tropicoperdix charltonii" # Tropicoperdix graydoni
                          )

# duplicate the rows
duplicated_spp <- functraits_global |>
  inner_join(no_traits_spp, by = c("Species2_eBird" = "Sister")) |> 
  mutate(Species2_eBird = scientific_name) |>
  select(-c(scientific_name, n))
  
functraits_global_update <- functraits_global |>
  bind_rows(duplicated_spp)

# check it correct this
list_species |> 
  filter(!scientific_name %in% functraits_global$Species2_eBird) |>
  nrow()
# vs
list_species |> 
  filter(!scientific_name %in% functraits_global_update$Species2_eBird) |>
  nrow()

# it worked!!!

functraits_islands <- functraits_global_update |>
  filter(Species2_eBird %in% list_species$scientific_name)

saveRDS(functraits_islands, file = "FunctionalTraits/FuncTraits_birds_Islands.rds")

## link records and functional traits ####

spp_cell_df_func <- spp_cell_df |>
  left_join(functraits_islands,
            join_by(scientific_name == Species2_eBird)) |>
  group_by(cell, phylo_name, scientific_name) |>
  # mean abundance per species per cell
  mutate(mu_hat = mean(mu_hat, na.rm = TRUE),
         omega_hat = mean(omega_hat, na.rm = TRUE)) |>
  as.data.frame() |>
  mutate(P.nes.str.cat = str_extract(nes.str.cat, "[^,]+"),
         P.nes.loc.cat = str_extract(nes.loc.cat, "[^,]+")) 

# However, there are some splits not reflected in this, like 'Troglodytes aedon'
spp_cell_df_func |> 
  filter(scientific_name == "Troglodytes aedon") |> 
  ggplot(aes(x = longitude,
             y = latitude, 
             fill = Name_USGSO)) + 
  facet_wrap(~subregion, ncol = 1) +
  geom_point(shape = 21) +
  geom_hline(yintercept = 23.5, color = "red")+
  theme_classic() +
  coord_equal()

# so, lets name these wren
# musculus at south of lat 23.5 and in mainland
spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & subregion == "Mainland",
                                  "Troglodytes musculus",
                                  scientific_name)) 

# Troglodytes beani in Cozumel island
spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO == "Cozumel",
                                  "Troglodytes beani",
                                  scientific_name)) 

# Troglodytes martinicensis in Dominica island
spp_cell_df_func <- spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO == "Dominica",
                                  "Troglodytes martinicensis",
                                  scientific_name))

# Troglodytes mesoleucus in Saint Lucia island
spp_cell_df_func <- spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO == "Saint Lucia",
                                  "Troglodytes mesoleucus",
                                  scientific_name))

# Troglodytes musicus in St Vincent island
spp_cell_df_func <- spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO == "St Vincent",
                                  "Troglodytes musicus",
                                  scientific_name))

# Troglodytes grenadensis in St Vincent island
spp_cell_df_func <- spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO == "Grenada",
                                  "Troglodytes grenadensis",
                                  scientific_name))

# musculus at south of lat 23.5 and in Continental islands, except Cozumel
spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Troglodytes aedon" & latitude < 23.5 & Name_USGSO %in% c("Barro Colorado",
                                                                                                               "Bastimentos",
                                                                                                               "Colon",
                                                                                                               "Tobago",
                                                                                                               "Trinidad"),
                                  "Troglodytes musculus",
                                  scientific_name)) 


spp_cell_df_func |> 
  filter(scientific_name %in% c("Troglodytes aedon",
                                "Troglodytes musculus",
                                "Troglodytes beani", 
                                "Troglodytes martinicensis",
                                "Troglodytes mesoleucus",
                                "Troglodytes musicus",
                                "Troglodytes grenadensis")) |> 
  ggplot(aes(x = longitude,
             y = latitude, 
             fill = scientific_name)) + 
  facet_wrap(~subregion, ncol = 1) +
  geom_point(shape = 21) +
  theme_classic() +
  coord_equal()

# And finally, the "White-breasted Thrasher (Ramphocinclus brachyurus)" split
spp_cell_df_func <- spp_cell_df_func |> 
  mutate(scientific_name = ifelse(scientific_name == "Ramphocinclus brachyurus",
                                  "Ramphocinclus sanctaeluciae",
                                  scientific_name))

length(unique(spp_cell_df_func$scientific_name))
# It increases the number of species.

# For some of these species, we should also correct the insular endemic trait (end.ins)

spp_cell_df_func <- spp_cell_df_func |>
  mutate(end.ins = ifelse(scientific_name %in% c("Troglodytes beani",
                                                 "Troglodytes grenadensis",
                                                 "Troglodytes martinicensis",
                                                 "Troglodytes mesoleucus"),
                          1, end.ins))

saveRDS(spp_cell_df_func, "Completeness_data_Islands/Records_eBird_func_traits.rds")

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

# Save data for PCAs figure!! ####

saveRDS(islands_cwm_mu, "Completeness_data_Islands/Islands_CWM_mu.rds")


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
#  scale_x_continuous(limits = c(-9,8))+
#  scale_y_continuous(limits = c(-4,6))+
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
#  scale_x_continuous(limits = c(-9,7))+
#  scale_y_continuous(limits = c(-5,4))+
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
  scale_x_continuous(limits = c(-9,8))+
  scale_y_continuous(limits = c(-5,5))+
  coord_fixed()+
  geom_segment(data = papuan.pca.loadings, 
               size = 0.25,
               aes(x = 0, xend = PC1*10, 
                   y = 0, yend = PC2*8),
               arrow = arrow(length = unit(0.1, "cm")),
               colour = "black") +
  geom_text_repel(data = papuan.pca.loadings, 
            aes(x = PC1*12, y = PC2*10, label = Variables)) + 
  labs(x = "PC1 (52.55%)",
       y = "PC2 (15.82%)",
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
grid.arrange(global.pca.density.fig, global.carib.pca.density.fig,
             global.orinma.pca.density.fig, global.papuan.pca.density.fig, 
             ncol = 2)

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
