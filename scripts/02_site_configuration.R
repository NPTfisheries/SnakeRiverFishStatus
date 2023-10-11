# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Download and configure Snake River IPTDS infrastructure for tag observation
# processing and the DABOM model
# 
# Created Date: October 10, 2023
#   Last Modified: 
#
# Notes: The output and saved file from this script is used for processing tag
# observations and for visualizing infrastructure.

#----------------------
# clear environment
rm(list = ls())

# install PITcleanr, if needed
# remotes::install_github("KevinSee/PITcleanr", ref = "develop")

# load needed libraries
library(here)
library(PITcleanr)
library(tidyverse)
library(sf)
library(magrittr)

#----------------------
# Snake River Sites of Interest

# get Snake River steelhead DPS polygon to filter area of interest
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop, spsm_pop)
sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_DPS = ESU_DPS, 
         sthd_MPG = MPG, 
         sthd_POP_NAME = POP_NAME, 
         sthd_TRT_POPID = TRT_POPID, 
         sthd_GSI_Group = GSI_Group); rm(sth_pop)

# plot sthd populations
sr_sthd_pops %>%
  ggplot() +
  geom_sf(aes(fill = sthd_MPG)) +
  theme_bw() +
  labs(title = "Snake River Steelhead DPS",
       fill = "MPG") +
  theme(legend.position = "bottom")

# query metadata for all INT and MRR sites in PTAGIS
org_config = buildConfig(node_assign = "array",
                         array_suffix = "UD")

# trim down to sites within the Snake River steelhead DPS
ptagis_sf = org_config %>%
  # clean things up a bit
  select(site_code,
         site_name,
         site_type,
         site_type_name,
         start_date,
         end_date,
         rkm,
         rkm_total,
         latitude,
         longitude,
         site_description) %>%
  # and make spatial
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", 
                      "latitude"), 
           crs = 4326) %>%
  # unique sites
  distinct(site_code,
           site_name,
           site_type,
           rkm,
           geometry,
           .keep_all = T) 

#----------------------
# Snake River INT sites
sr_int_sites_sf = ptagis_sf %>%
  # trim down to sites within the Snake River steelhead DPS
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) %>%
  # grab only INT sites for now 
  filter(site_type == "INT") %>%
  # remove some dam sites; we'll deal with those later
  filter(!str_detect(site_name, "Lower Granite")) %>%
  filter(!str_detect(site_name, "LOWER GRANITE")) %>%
  filter(!str_detect(site_name, "Little Goose")) %>%
  filter(!str_detect(site_name, "Lower Monumental")) %>%
  # remove unnecessary INT sites we don't want
  filter(!site_code %in% c("CCP", # Catherine Creek Acclimation Pond
                           "GRP", # Grande Ronde Acclimation Pond
                           "LOP", # Lostine River Acclimation Pond
                           "RPJ", # Rapid River Hatchery Pond
                           "S2I", # Lemhi Sub-reach 2 SC Inlet
                           "S2O", # Lemhi Sub-reach 2 SC Outlet
                           "S3A", # Eagle Valley Ranch S3A
                           "S3B", # Eagle Valley Ranch S3B
                           "CHN", # Challis Diversion North
                           "CHS", # Challis Diversion South
                           "CLJ", # Clearwater River Juvenile Fish Trap
                           "IMJ", # Imnaha River Juvenile Fish Trap
                           "SAJ", # Salmon River Trap
                           "SNJ", # Snake River Trap
                           "LGW")) # Lookingglass Creek weir (not sure why this wasn't in previous models)
# "PWA" (Penawawa Creek) wasn't included in previous DABOM model runs; adding in for now

#----------------------
# Snake River MRR Sites
# read in complete tag histories since SY2010
cth_df = list.files(path = here("data/complete_tag_histories/"),
                    pattern = "\\.csv$",
                    full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(.x), .id = "file_name") %>%
  # add species and spawn year
  mutate(file_name = str_replace(file_name, ".*/", ""), 
         species = str_extract(file_name, "(?<=_)[^_]+"),       
         spawn_year = str_extract(file_name, "SY[0-9]{4}")) %>%
  select(-file_name)

# summarize number of tags observed by site
tags_by_site = cth_df %>%
  select(species,
         spawn_year,
         site_code = `Event Site Code Value`,
         tag_code = `Tag Code`) %>%
  distinct() %>%
  group_by(species,
           spawn_year,
           site_code) %>%
  summarise(n_tags = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = c(species, spawn_year),
    values_from = n_tags) %>%
  mutate(n_chnk_tags = rowSums(select(., starts_with("Chinook")), na.rm = TRUE),
         n_sthd_tags = rowSums(select(., starts_with("Steelhead")), na.rm = TRUE))

# join detections to MRR sites
sr_mrr_sites_sf = ptagis_sf %>%
  # trim down to sites within the Snake River steelhead DPS
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) %>%
  # grab only MRR sites
  filter(site_type == "MRR") %>%
  left_join(tags_by_site %>%
              select(site_code,
                     n_chnk_tags,
                     n_sthd_tags)) %>%
  # filter for only sites where there have been 100 or greater detections for Chinook salmon or steelhead
  filter(!n_chnk_tags <= 100 | !n_sthd_tags <= 100) %>%
  select(-n_chnk_tags,
         -n_sthd_tags) %>%
  # remove Lower Granite Dam MRR sites
  filter(!str_detect(site_code, "LGR")) %>%
  # remove acclimation ponds 
  filter(!str_detect(site_description, "AcclimationPond")) %>%
  # Do I want to leave out all hatcheries? These detections seem unreliable i.e., what's their disposition?
  filter(!str_detect(site_description, "Hatchery")) %>%
  # finally, add back in some MRR sites that have been included in previous DABOM model runs. Do we really need all of these?
  bind_rows(ptagis_sf %>%
              filter(site_code %in% c("ALMOTC",   # Almota Creek - tributary to Snake River
                                      "BCANF",    # Big Canyon Facility
                                      "BEARVC",   # Bear Valley Creek
                                      "BIG2C",    # Big Creek, Middle Fork Salmon River
                                      "BIGBEC",   # Big Bear Creek, Potlatch River
                                      "CAMP4C",   # Camp Creek, tributary to Big Sheep Creek, Imnaha drainage
                                      "CARMEC",   # Carmen Creek - tributary to Salmon River
                                      "CHARLC",   # Charley Creek, Asotin Creek watershed
                                      "CLWH",     # Clearwater Hatchery
                                      "CROTRP",   # Crooked River Trap
                                      "DRY2C",    # Dry Creek - tributary to Imnaha River
                                      "DWOR",     # Dworshak National Fish Hatchery
                                      "EFPW",     # East Fork Potlatch River weir
                                      "FREEZC",   # Freezeout Creek - tributary to Imnaha River
                                      "GUMBTC",   # Gumboot Creek, Imnaha River Basin
                                      "HORS3C",   # Horse Creek, Imnaha River Basin
                                      "KOOS",     # Kooskia National Fish Hatchery
                                      "LBEARC",   # Little Bear Creek, Potlatch River watershed
                                      "LOOH",     # Lookingglass Hatchery
                                      "LSHEEF",   # Little Sheep Facility
                                      "LYFE",     # Lyons Ferry Hatchery
                                      "MAHOGC",   # Mahogany Creek, Imnaha River Basin
                                      "MCCA",     # McCall Hatchery
                                      "NPTH",     # Nez Perce Tribal Hatchery
                                      "OXBO",     # Oxbow Hatchery (IDFG)
                                      "PAHH",     # Pahsimeroi Hatchery
                                      "PENAWC",   # Penawawa Creek - tributary to Snake River
                                      "POTREF",   # East Fork Potlatch River
                                      "POTRWF",   # West Fork Potlatch River
                                      "RAPH",     # Rapid River Hatchery
                                      "REDR",     # Red River
                                      "REDTRP",   # Red River Trap
                                      "RPDTRP",   # Rapid River Smolt Trap
                                      "SAWT",     # Sawtooth Hatchery
                                      "TENMC2",   # Tenmile Creek, tributary to Snake River
                                      "TUCH",     # Tucannon River Hatchery
                                      "WALH",     # Wallowa Hatchery
                                      "WIMPYC"))) # Wimpey Creek (Salmon River)

#----------------------
# build configuration file
sr_sites_list = bind_rows(sr_int_sites_sf,
                          sr_mrr_sites_sf) %>%
  st_drop_geometry() %>%
  select(site_code) %>%
  distinct() %>%
  pull()

sr_config = org_config %>%
  filter(site_code %in% sr_sites_list) %>%
  # THIS IS WHERE TO FIX SNAKE RIVER SITES
  
  
dam_config = org_config %>%
  # filter our just relevant dam sites
  filter(site_code %in% c("GRA", "LGRLDR", "LGR",            # LGR Adults
                          "GRJ", "GRS", "GRX",               # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
                          "LGS", "GOJ", "GOA",               # Little Goose Dam
                          "LMA", "LMN", "LMJ",               # Lower Monumental Dam
                          "IHA", "IHR", "ICH",               # Ice Harbor Dam
                          "MCN", "MCJ", "MCX", "MC1", "MC2", # McNary Dam
                          "JDA", "JDJ", "JO1", "JO2",        # John Day Dam
                          "TDA", "TD1", "TD2",               # The Dalles Dam
                          "B1J", "B2J", "B2A", "BVJ", "BVX", # Bonneville Dam
                          "BO1", "BO2", "BO3", "BO4",
                          "BON", "BCC", "BWL", "BONAFF")) %>%
  # combine each Lower Snake and Columbia dam into a single node; upstream to downstream
  mutate(
    node = case_when(
      site_code %in% c("GRA", "LGRLDR", "LGR")            ~ "LGR", # LGR Adults
      site_code %in% c("GRJ", "GRS", "GRX")               ~ "GRS", # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
      site_code %in% c("LGS", "GOJ", "GOA")               ~ "GOA", # Little Goose Dam
      site_code %in% c("LMA", "LMN", "LMJ")               ~ "LMA", # Lower Monumental Dam
      site_code %in% c("IHA", "IHR", "ICH")               ~ "IHR", # Ice Harbor Dam
      site_code %in% c("MCN", "MCJ", "MCX", "MC1", "MC2") ~ "MCN", # McNary Dam
      site_code %in% c("JDA", "JDJ", "JO1", "JO2")        ~ "JDA", # John Day Dam
      site_code %in% c("TDA", "TD1", "TD2")               ~ "TDA", # The Dalles Dam
      site_code %in% c("B1J", "B2J", "B2A", "BVJ", "BVX",          # Bonneville Dam
                       "BO1", "BO2", "BO3", "BO4",
                       "BON", "BCC", "BWL", "BONAFF")     ~ "BON",
      TRUE ~ node))

downriver_config = org_config %>%
  mutate(
    node = case_when(
      # recode primary subbasins and tributaries outside the Snake River basin; using the "\\d+" ensures that we search for
      # all consecutive digits prior to a "." i.e., some sites in the Upper Columbia have a 4-digit start to their rkm.
      # Then "\\.\\d+" ensures we search for consecutive digits after a ".".
      as.numeric(str_extract(rkm, "\\d+")) > 539                                                             ~ "PRA", # Upper Columbia
      as.numeric(str_extract(rkm, "\\d+")) == 539 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "PRO", # Yakima
      as.numeric(str_extract(rkm, "\\d+")) == 509 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WWB", # Walla Walla
      as.numeric(str_extract(rkm, "\\d+")) == 465 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "UMW", # Umatilla
      as.numeric(str_extract(rkm, "\\d+")) == 351 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "JD1", # John Day
      as.numeric(str_extract(rkm, "\\d+")) == 328 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "DRM", # Deschutes
      as.numeric(str_extract(rkm, "\\d+")) == 290 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "KLR", # Klickitat
      as.numeric(str_extract(rkm, "\\d+")) == 273 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "FID", # Hood River
      as.numeric(str_extract(rkm, "\\d+")) == 271 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "RCX", # White Salmon
      as.numeric(str_extract(rkm, "\\d+")) == 261 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "LWL", # Little White Salmon
      as.numeric(str_extract(rkm, "\\d+")) == 251 & as.numeric(sub(".","", str_extract(rkm, "\\.\\d+"))) > 0 ~ "WRA", # Wind River
      TRUE ~ node)) %>%
  filter(node %in% c(
    "PRA", "PRO", "WWB", "UMW", "JD1", "DRM", "KLR", "FID", "RCX", "LWL", "WRA"
  ))

# sites list used to create parent-child table, make spatial
sites_sf = bind_rows(sr_int_sites_sf,
                     sr_mrr_sites_sf) %>%
  select(-contains("sthd")) %>%
  # add back in dams
  bind_rows(ptagis_sf %>%
              filter(site_code %in% c("LGR",       # LGR Adults
                                      "GRS",       # LGR Juveniles; GRJ, GRX = LGR Juvenile Bypass; GRS = LGR Spillway
                                      "GOA",       # Little Goose Dam
                                      "LMA",       # Lower Monumental Dam
                                      "IHR",       # Ice Harbor Dam
                                      "MCN",       # McNary Dam
                                      "JDA",       # John Day Dam
                                      "TDA",       # The Dalles Dam
                                      "BON"))) %>% # Bonneville Dam
  # add back in lowest ptagis sites in downriver subbasins
  bind_rows(ptagis_sf %>%
              filter(site_code %in% c("PRA",       # Upper Columbia
                                      "PRO",       # Yakima
                                      "WWB",       # Walla Walla
                                      "UMW",       # Umatilla
                                      "JD1",       # John Day
                                      "DRM",       # Deschutes
                                      "KLR",       # Klickitat
                                      "FID",       # Hood River
                                      "RCX",       # White Salmon
                                      "LWL",       # Little White Salmon
                                      "WRA"))) %>% # Wind River 
  # trim to unique sites
  distinct(site_code,
           site_name,
           site_type,
           site_type_name,
           rkm,
           site_description,
           geometry,
           .keep_all = T) %>%
  # remove an extra Prosser Dam record
  filter(!(site_code == "PRO" & site_type == "MRR")) %>%
  st_transform(crs = 5070) # NAD83

#----------------------
# download the NHDPlus v2 Flowlines
dwn_flw = T
nhd_list = queryFlowlines(sites_sf = sites_sf,
                          root_site_code = "LGR",
                          min_strm_order = 2,
                          dwnstrm_sites = dwn_flw, # do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
                          dwn_min_stream_order_diff = 4)

# compile the upstream and downstream flowlines
flowlines = nhd_list$flowlines
if(dwn_flw) {
  flowlines %<>%
    rbind(nhd_list$dwn_flowlines)
}

# now cut off areas too far upstream
upstrm_loc = "Hells Canyon Dam" # upstream extent of study area 
library(ggmap)
upstrm_comid = ggmap::geocode(upstrm_loc, output = "latlon") %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) %>%
  nhdplusTools::discover_nhdplus_id()
 
nhd_upstrm_lst = nhdplusTools::plot_nhdplus(outlets = list(upstrm_comid),
                                            streamorder = min(nhd_list$flowlines$StreamOrde),
                                            actually_plot = F) 

flowlines %<>%
  anti_join(nhd_upstrm_lst$flowline %>%
              st_drop_geometry() %>%
              select(Hydroseq))

#----------------------
# plot the flowlines and sites
library(ggrepel)
site_map = ggplot() +
  geom_sf(data = sr_sthd_pops,
          aes(fill = sthd_MPG)) +
  geom_sf(data = flowlines,
          aes(color = as.factor(StreamOrde),
              size = StreamOrde)) +
  scale_color_viridis_d(direction = -1,
                        option = "D",
                        name = "Stream\nOrder",
                        end = 0.8) +
  scale_size_continuous(range = c(0.2, 1.2),
                        guide = 'none') +
  geom_sf(data = sites_sf,
          size = 3, 
          color = "black") +
  ggrepel::geom_label_repel(
    data = sites_sf %>%
      filter(site_code != "LGR"),
    aes(label = site_code,
        geometry = geometry),
    size = 2,
    stat = "sf_coordinates",
    min.segment.length = 0,
    max.overlaps = 100) +
  geom_sf_label(data = sites_sf %>%
                  filter(site_code == "LGR"),
                aes(label = site_code),
                color = "red") +
  theme_void() +
  labs(fill = "Steelhead\nMPG") +
  theme(axis.title = element_blank(),
        legend.position = "bottom")
site_map

# save site map
ggsave(here("output/figures/full_site_map.png"),
       site_map,
       width = 14,
       height = 8.5)


