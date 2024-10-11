# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Download and configure Snake River IPTDS infrastructure for tag observation
#   processing and the DABOM model. The outputs from this script are used for processing
#   tag observations and visualizing infrastructure.
# 
# Created Date: October 10, 2023
#   Last Modified: October 2, 2024
#
# Notes: 

#----------------------
# clear environment
rm(list = ls())

# install PITcleanr, if needed
# remotes::install_github("KevinSee/PITcleanr", ref = "develop")

# load needed libraries
library(PITcleanr)
library(tidyverse)
library(here)
library(sf)
library(magrittr)
library(janitor)
library(readxl)

# set default coordinate reference system
default_crs = st_crs(32611) # WGS84 ; UTM zone 11N

#----------------------
# prep to summarize Snake River sites of interest

# get Snake River steelhead DPS polygon to filter area of interest
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop)
sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_dps = ESU_DPS,
         sthd_mpg = MPG,
         sthd_popid = TRT_POPID,
         sthd_popname = POP_NAME) %>%
  st_transform(default_crs); rm(sth_pop)

# plot sthd populations
sthd_pops %>%
  ggplot() +
  geom_sf(aes(fill = sthd_mpg)) +
  theme_bw() +
  labs(title = "Snake River Steelhead DPS",
       fill = "MPG") +
  theme(legend.position = "bottom")

# query metadata for all INT and MRR sites in PTAGIS
org_config = buildConfig(node_assign = "array",
                         array_suffix = "UD")

# trim configuration file down to unique sites and make spatial
ptagis_sf = org_config %>%
  # select only columns of interest
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
  # make spatial
  filter(!is.na(latitude) | !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", 
                      "latitude"), 
           crs = 4326) %>%
  st_transform(default_crs) %>%
  # unique sites
  distinct(site_code,
           site_name,
           site_type,
           rkm,
           geometry,
           .keep_all = T) 

#----------------------
# create list of Snake River INT sites
sr_int_sites_sf = ptagis_sf %>%
  # trim down to sites within the Snake River steelhead DPS
  st_join(sthd_pops) %>%
  filter(!is.na(sthd_dps)) %>%
  # grab only INT sites for now 
  filter(site_type == "INT") %>%
  # remove some dam sites; we'll deal with those later
  filter(!str_detect(site_name, "Lower Granite|LOWER GRANITE|Little Goose|Lower Monumental")) %>%
  # remove unnecessary INT sites we don't want
  filter(!site_code %in% c("0HR", # Henry's Inreach Array, Lemhi
                           "CCP", # Catherine Creek Acclimation Pond
                           "CHN", # Challis Diversion North
                           "CHS", # Challis Diversion South
                           "CLJ", # Clearwater River Juvenile Fish Trap
                           "COU", # Couse Creek Near Mouth
                           "CRT", # Crooked River Satellite Facility
                           # Dworshak natural-origin fish are returned back to river and so we don't 
                           # really know their destination unless they are detected elsewhere which is 
                           # not necessarily the case for other "Hatchery Return" facilities (e.g., CRT, 
                           # RRT, STL, STR) where we can be confident that they stay in that tributary
                           "DWL", # Dworshak NFH Adult Trap
                           "GRP", # Grande Ronde Acclimation Pond
                           "IMJ", # Imnaha River Juvenile Fish Trap
                           "LOP", # Lostine River Acclimation Pond
                           "RPJ", # Rapid River Hatchery Pond
                           "RRT", # Red River Satellite Facility
                           "S2I", # Lemhi Sub-reach 2 SC Inlet
                           "S2O", # Lemhi Sub-reach 2 SC Outlet
                           "S3A", # Eagle Valley Ranch S3A
                           "S3B", # Eagle Valley Ranch S3B
                           "SAJ", # Salmon River Trap
                           "SNJ"))# Snake River Trap

#----------------------
# create list of Snake River MRR Sites
# read in complete tag histories since SY2010
tags_by_site = list.files(path = here("data/complete_tag_histories/"),
                          pattern = "\\.csv$",
                          full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(.x, show_col_types = F), .id = "file_name") %>%
  clean_names() %>%
  # add species and spawn year
  mutate(file_name = str_replace(file_name, ".*/", ""), 
         species = str_extract(file_name, "(?<=_)[^_]+"),       
         spawn_year = str_extract(file_name, "SY[0-9]{4}")) %>%
  select(-file_name) %>%
  # summarize number of unique tags observed by site
  select(species,
         spawn_year,
         site_code = event_site_code_value,
         tag_code) %>%
  distinct() %>%
  group_by(species,
           spawn_year,
           site_code) %>%
  summarise(n_tags = n(), .groups = "drop") %>%
  group_by(species,
           site_code) %>%
  summarise(n_yrs = n(),
            n = sum(n_tags),
            min = round(min(n_tags, na.rm = T), 1),
            mean = round(mean(n_tags, na.rm = T), 1),
            md = round(median(n_tags, na.rm = T), 1),
            max = round(max(n_tags, na.rm = T), 1),
            yrs = toString(spawn_year),
            .groups = "drop")

# mrr sites of interest
sr_mrr_sites_sf = ptagis_sf %>%
  # trim down to sites within the Snake River steelhead DPS
  st_join(sthd_pops) %>%
  filter(!is.na(sthd_dps)) %>%
  # grab only MRR sites
  filter(site_type == "MRR") %>%
  # remove any sites below the Tucannon River mouth
  filter(str_extract(rkm, "^[0-9]{3}") == "522" & as.numeric(str_extract(rkm, "(?<=\\.)[0-9]{3}")) >= 100) %>%
  # remove any Lower Granite or Little Goose Dam MRR sites
  filter(!str_detect(site_code, "LGR|LGS")) %>%
  # remove acclimation ponds or river segments
  filter(!str_detect(site_description, "AcclimationPond|RiverSegment")) %>%
  # join unique tags observed by site
  left_join(tags_by_site %>%
              select(species, site_code, mean) %>%
              pivot_wider(names_from = species,
                          values_from = mean)) %>%
  # first, let's trim down our list of MRR sites significantly by filtering for only sites where, on average, there 
  # are 10 or greater detections for Chinook salmon, steelhead, or coho (among only years where fish were detected
  # i.e., na.rm = T)
  filter(Chinook >= 10 | Steelhead >= 10 | Coho >= 10) %>%
  select(-Chinook, -Steelhead, -Coho) %>%
  # now, remove some additional sites that pass the above filter, but are not necessarily useful for escapement analyses.
  # this can be because we can't safely assign detections to other INT or MRR locations, or that carcass recoveries maybe
  # occurred prior to an adjacent array and detections muddle interpretation of escapement estimates there
  filter(!site_code %in% c("ASOTIC",  # Asotin Creek
                           "CATHEC",  # Catherine Creek (not certain whether these could be assigned to CCW_U or not)
                           "FISHC",   # Fish Creek carcass recoveries (occurred prior to LRU being installed)
                           "FISTRP",  # Fish Creek trap (occured prior to LRU being installed)
                           "GEORGC",  # George Creek, Asotin Creek watershed
                           "HORS3C",  # Horse Creek, Imnaha River Basin
                           "PANTHC",  # Panther Creek carcass recoveries (some occurred prior to PCA being installed)
                           "SALR4",   # Salmon River - Pahsimeroi River to headwaters (km 489-650)
                           "SALRSF"   # South Fork Salmon River
                           )) %>%
  # finally, add back in some MRR sites that are useful i.e., they help alleviate a black box, can easily be grouped 
  # to a weir to help a conservative estimate, or have sufficient, regular detections that they can be used to estimate 
  # a detection probability at a downstream site, etc.
  bind_rows(ptagis_sf %>%
              filter(site_code %in% c(# UPPER SALMON
                                      "ALTULC",       # Alturas Lake Creek, Upper Salmon
                                      "BEAVEC",       # Beaver Creek, Upper Salmon
                                      "PAHSIR",       # Pahsimeroi River
                                      "SAWTRP",       # Sawtooth Trap
                                      # MIDDLE FORK SALMON
                                      "BEARVC",       # Bear Valley Creek
                                      # SOUTH FORK SALMON
                                      "BURNLC",       # Burntlog Creek, Johnson Creek
                                      "GROUSC",       # Grouse Creek
                                      "LAKEC",        # Lake Creek
                                      "MCCA",         # McCall Hatchery
                                      "SUMITC",       # Summit Creek
                                      # LITTLE SALMON
                                      "RAPH",         # Rapid River Hatchery
                                      # SOUTH FORK CLEARWATER
                                      "KOOS",         # Kooskia National Fish Hatchery
                                      # IMNAHA RIVER
                                      "LSHEEF",       # Little Sheep Facility
                                      # GRANDE RONDE
                                      "BCANF",        # Big Canyon Facility
                                      "WALH",         # Wallowa Hatchery
                                      # LOWER SNAKE
                                      "PENAWC",       # Penawawa Creek - tributary to Snake River
                                      "TENMC2"))) %>% # Tenmile Creek, tributary to Snake River
  arrange(site_code)

#----------------------
# build configuration file
# list of Snake River INT and MRR sites of interest to trim org_config
sr_sites_list = bind_rows(sr_int_sites_sf,
                          sr_mrr_sites_sf) %>%
  st_drop_geometry() %>%
  select(site_code) %>%
  distinct() %>%
  pull()

# Snake River INT and MRR Sites
sr_config = org_config %>%
  filter(site_code %in% sr_sites_list) %>%
  # Fix some INT sites
  mutate(
    node = case_when(
      site_code == "BTC"                 ~ str_replace(node, "BTC", "BTL"), # Big Timber Creek; BTC was replaced by BTL
      site_code == "18M"                 ~ str_replace(node, "18M", "HEC"), # Group 18-mile & Hawley creeks
      site_code == "IML" & antenna_id == "09"                    ~ "IML_U", # Imnaha River Work Room Antenna
      site_code == "AFC" &  grepl("MAINSTEM", antenna_group)     ~ "AFC_D", # mainstem Asotin becomes _D
      site_code == "AFC" & !grepl("MAINSTEM", antenna_group)     ~ "AFC_U", # south and north forks become _U
      TRUE ~ node)) %>%  
  # Recode and/or merge some MRR sites; often merging MRR observations into INT sites
  mutate(
    node = case_when(
      # UPPER SALMON 
      site_code %in% c("ALTULC", "BEAVEC",  
                       "SAWTRP", "STL")   ~ "SAWT",   # Group Sawtooth Hatchery/Ladder Array & Upper Salmon carcass recoveries all to SAWT
      site_code %in% c("CEY", "YANKFK")   ~ "YFK_U",  # Group Yankee Fork and Cearley Creek obs to YFK_U
      site_code == "SALREF"               ~ "SALEFT", # Group EF Salmon River obs (e.g., carcass recoveries) w/ trap
      site_code == "PAHSIR"               ~ "PAHH"  , # Group Pahsimeroi River carcass recoveries to PAHH
      # MIDDLE FORK SALMON
      site_code == "BEARVC"               ~ "BRC",    # Group Bear Valley Creek carcass recoveries to BRC
      # SOUTH FORK SALMON
      site_code %in% c("KNOXB", "MCCA", 
                       "STR")             ~ "SALSFW",   # South Fork Salmon River weir
      site_code %in% c("SECESR", "GROUSC",
                       "SUMITC", "LAKEC") ~ "ZEN_U",  # Group Secesh River obs (e.g., carcass recoveries) w/ ZEN_U
      site_code == "BURNLC"               ~ "JOHNSC", # Group Burntlog Creek to Johnson Creek
      # SOUTH FORK CLEARWATER
      site_code == "KOOS"                 ~ "CLC_U",    # Group Kooskia Hatchery to CLC_U
      # IMNAHA RIVER
      site_code == "IMNAHW"               ~ "IML_U",   # Group Imnaha Weir w/ trap array
      # GRANDE RONDE
      site_code == "GRANDW"               ~ "UGS_U",   # Group GRANDW to UGS_U (based on site configuration)
      site_code == "CATHEW"               ~ "CCW_U",   # Group CATHEW to CCW_U (based on site configuration)
      site_code %in% c("LOOKGC", "LOOH")  ~ "LGW_U",   # Group LGW, LOOKGC, and LOOH into a single node
      site_code == "JOSEPC"               ~ "JOC_U",
      # LOWER SNAKE
      site_code == "TUCH"                 ~ "TFH_U",   # Group Tucannon Hatchery to TFH_U
      site_code == "PENAWC"               ~ "PWA_D",   # Group Penawawa Creek to PWA_U
      TRUE ~ node))

# Snake and Columbia River dams configuration 
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

# Downriver Subbasins configuration
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

# merge into a single configuration file
configuration = bind_rows(sr_config,
                          dam_config,
                          downriver_config) ; rm(sr_config, dam_config, downriver_config)

# convert configuration into a sites_sf
crb_sites_sf = configuration %>%
  select(site_code = node, org_site_code = site_code) %>%
  mutate(site_code = str_replace(site_code, "_D$|_U$", "")) %>%
  distinct() %>%
  group_by(site_code) %>%
  summarize(incl_sites = toString(org_site_code)) %>%
  ungroup() %>%
  mutate(incl_sites = str_remove_all(incl_sites, paste0("(^|,\\s*)", site_code, "(\\s*,|$)")),
         incl_sites = str_replace_all(incl_sites, ",{2,}", ","),
         incl_sites = str_replace_all(incl_sites, "^,|,$", ""),
         incl_sites = str_trim(incl_sites)) %>%
  left_join(ptagis_sf, by = "site_code") %>%
  # remove an extra Prosser Dam record
  filter(!(site_code == "PRO" & site_type == "MRR")) %>%
  st_as_sf(crs = default_crs) %>%
  select(site_code,
         site_name,
         site_type,
         site_type_name,
         incl_sites,
         everything())

#----------------------
# download the NHDPlus v2 Flowlines
dwn_flw = T # do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
nhd_list = queryFlowlines(sites_sf = crb_sites_sf,
                          root_site_code = "LGR",
                          min_strm_order = 2,
                          dwnstrm_sites = dwn_flw, 
                          dwn_min_stream_order_diff = 4)

# compile the upstream and downstream flowlines
flowlines = nhd_list$flowlines
if(dwn_flw == T) {
  flowlines %<>%
    rbind(nhd_list$dwn_flowlines)
}

# now cut off areas too far upstream
library(ggmap)
upstrm_loc = "Hells Canyon Dam" # upstream extent of study area 
upstrm_comid = ggmap::geocode(upstrm_loc, output = "latlon") %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) %>%
  nhdplusTools::discover_nhdplus_id()
 
nhd_upstrm_lst = nhdplusTools::plot_nhdplus(outlets = list(upstrm_comid),
                                            streamorder = min(nhd_list$flowlines$StreamOrde),
                                            actually_plot = F) 

# cut off flowlines upstream of Hells Canyon Dam
flowlines %<>%
  anti_join(nhd_upstrm_lst$flowline %>%
              st_drop_geometry() %>%
              select(Hydroseq))

#----------------------
# build parent-child table
parent_child = crb_sites_sf %>%
  # just keep values at and upstream of LGR for the initial parent-child table
  filter(
    grepl("^522\\.[0-9]{3}", rkm) & 
      as.numeric(sub("^522\\.([0-9]{3}).*", "\\1", rkm)) >= 173
  ) %>% 
  filter(!site_code == "GRS") %>%
  buildParentChild(flowlines = flowlines,
                   rm_na_parent = FALSE,
                   add_rkm = FALSE) %>%
  # fix site mapping issues above LGR
  editParentChild(fix_list = 
                    list(# Lemhi River
                         c("KEN", "AGC", "EVU"),
                         c("LLS", "LBS", "LRW"),
                         c("LLS", "LB8", "LRW"),
                         c("LLS", "LCL", "LRW"),
                         c("LBS", "HEC", "LRW"),
                         c("LBS", "BTL", "LRW"),
                         c("LBS", "CAC", "LRW"),
                         # Imnaha River
                         c("IR3", "IML", "IR4"),
                         c("IR4", "IR5", "IML"))) %>%
  filter(!is.na(parent)) %>%
  # now add in parent-child relationships below LGR
  bind_rows(
    tribble(
      ~parent, ~child,
      "LGR", "GRS",    # Lower Granite Spillway
      "GRS", "PWA",    # Penawawa Creek
      "GRS", "GOA",    # Little Goose Dam
      "GOA", "LTR",    # Tucannon River
      "LTR", "MTR",
      "MTR", "UTR",
      "UTR", "TFH",
      "TFH", "TPJ",
      "GOA", "LMA",    # Lower Monumental Dam
      "LMA", "IHR",    # Ice Harbor Dam
      "IHR", "PRA",    # Upper Columbia River
      "IHR", "PRO",    # Yakima River
      "IHR", "WWB",    # Walla Walla River
      "IHR", "MCN",    # McNary Dam
      "MCN", "UMW",    # Umatilla Dam
      "MCN", "JD1",    # John Day River
      "MCN", "JDA",    # John Day Dam
      "JDA", "DRM",    # Deschutes River
      "JDA", "TDA",    # The Dalles Dam
      "TDA", "KLR",    # Klickitat River
      "TDA", "FID",    # Hood River
      "TDA", "RCX",    # White Salmon River
      "TDA", "LWL",    # Little White Salmon River
      "TDA", "WRA",    # Wind River
      "TDA", "BON")) %>% # Bonneville Dam
  select(-parent_hydro,
         -child_hydro)

#----------------------
# population info for snake river dabom sites
sr_site_pops = crb_sites_sf %>%
  # filter for sites at or above LTR
  filter(str_extract(rkm, "^[0-9]{3}") == "522" & as.numeric(str_extract(rkm, "(?<=\\.)[0-9]{3}")) >= 100) %>%
  # remove GOA, LGR, and GRS
  filter(!site_code %in% c("GOA", "LGR", "GRS")) %>%
  select(site_code, site_type, incl_sites) %>%
  # join steelhead, sp/sum chinook, and (made up) coho populations
  st_join(sthd_pops %>%
            select(-sthd_dps)) %>%
  st_join(spsm_pop %>%
            st_transform(crs = default_crs) %>%
            select(chnk_mpg = MPG,
                   chnk_popid = TRT_POPID,
                   chnk_popname = POP_NAME)) %>%
  left_join(read_excel(here("data/coho_populations/coho_populations.xlsx")) %>%
              select(-coho_esu_dps)) %>%
  # move geometry column to the end
  select(-geometry, everything(), geometry) %>%
  mutate(
    chnk_popid   = if_else(site_code %in% c("SC1", "SC2"), "SCUMA", chnk_popid),
    chnk_popname = if_else(site_code %in% c("SC1", "SC2"), "Upper South Fork Clearwater", chnk_popname)
  ) %>%
  mutate(
    chnk_popid   = if_else(site_code %in% c("IR1", "IR2"), "IRMAI/IRBSH", chnk_popid), # We don't necessarily know whether IR1 and IR2 Chinook end up in IRMAI or IRBSH
    chnk_popname = if_else(site_code %in% c("IR1", "IR2"), "Imnaha River mainstem/Big Sheep Creek", chnk_popname)
  ) %>%
  mutate(
    chnk_popid   = if_else(site_code %in% c("SW1", "SW2"), "SEUMA/SEMEA/SEMOO", chnk_popid),
    chnk_popname = if_else(site_code %in% c("SW1", "SW2"), "Upper Selway River/Meadow Creek/Moose Creek", chnk_popname)
  ) %>%
  mutate(
    chnk_popid   = if_else(site_code == "WR1", "GRLOS/GRMIN", chnk_popid),
    chnk_popname = if_else(site_code == "WR1", "Lostine River/Minam River", chnk_popname)
  ) %>%
  mutate(
    chnk_popid   = if_else(site_code == "SFG", "SFSMA/SFSEC/SFEFS", chnk_popid),
    chnk_popname = if_else(site_code == "SFG", "South Fork Salmon River mainstem/Secesh River/East Fork South Fork Salmon River", chnk_popname)
  ) %>%
  mutate(
    chnk_popid   = if_else(site_code %in% c("USE", "USI"), "SRLMA/SRPAH/SREFS/SRYFS/SRVAL/SRUMA", chnk_popid),
    chnk_popname = if_else(site_code %in% c("USE", "USI"), "Salmon River lower mainstem below Redfish Lake/
                           Pahsimeroi River/East Fork Salmon River/Yankee Fork/Valley Creek/Salmon River upper
                           mainstem above Redfish Lake", chnk_popname)
  ) %>%
  mutate(
    sthd_popid   = if_else(site_code %in% c("SC1", "SC2"), "CRSFC-s", sthd_popid),
    sthd_popname = if_else(site_code %in% c("SC1", "SC2"), "South Fork Clearwater River", sthd_popname)
  ) %>%
  mutate(
    sthd_popid   = if_else(site_code %in% c("USE", "USI"), NA, sthd_popid), # We don't necessarily know which population USI, USE steelhead end up in
    sthd_popname = if_else(site_code %in% c("USE", "USI"), NA, sthd_popname)
  )

#----------------------
# write configuration, parent-child table, flowlines, etc.
save(configuration,
     crb_sites_sf,
     sr_site_pops,
     flowlines,
     parent_child,
     file = here("data/configuration_files/site_config_LGR_20240927.rda"))

# write sites_sf and flowlines out to geopackage, if desired
st_write(crb_sites_sf, dsn = "data/spatial/dabom_sites.gpkg", layer = "sites_sf", driver = "GPKG", append = F)
st_write(flowlines, dsn = "data/spatial/dabom_sites.gpkg", layer = "flowlines", driver = "GPKG", append = F)

# END SCRIPT