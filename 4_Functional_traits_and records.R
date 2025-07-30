# Functional traits and the records

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



