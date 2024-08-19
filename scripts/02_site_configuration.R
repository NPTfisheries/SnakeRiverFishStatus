# -----------------------
# Author(s): Mike Ackerman, Kevin See, and Ryan Kinzer
# Purpose: Download and configure Snake River IPTDS infrastructure for tag observation
#   processing and the DABOM model. The outputs from this script are used for processing
#   tag observations and visualizing infrastructure.
# 
# Created Date: October 10, 2023
#   Last Modified: August 19, 2024
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

#----------------------
# prep to summarize Snake River sites of interest

# get Snake River steelhead DPS polygon to filter area of interest
load(here("data/spatial/SR_pops.rda")) ; rm(fall_pop, spsm_pop)
sr_sthd_pops = st_as_sf(sth_pop) %>%
  select(sthd_dps = ESU_DPS,
         sthd_mpg = MPG,
         sthd_popid = TRT_POPID,
         sthd_popname = POP_NAME); rm(sth_pop)

# plot sthd populations
sr_sthd_pops %>%
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
  st_join(sr_sthd_pops) %>%
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
  # CONTINUE HERE!!!
  # add species and spawn year
  mutate(file_name = str_replace(file_name, ".*/", ""), 
         species = str_extract(file_name, "(?<=_)[^_]+"),       
         spawn_year = str_extract(file_name, "SY[0-9]{4}")) %>%
  select(-file_name) %>%
  # summarize number of unique tags observed by site
  select(species,
         spawn_year,
         site_code = `Event Site Code Value`,
         tag_code = `Tag Code`) %>%
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
  st_join(sr_sthd_pops) %>%
  filter(!is.na(sthd_DPS)) %>%
  # grab only MRR sites
  filter(site_type == "MRR") %>%
  # join unique tags observed by site
  left_join(tags_by_site %>%
              select(species, site_code, mean) %>%
              pivot_wider(names_from = species,
                          values_from = mean)) %>%
  # first, let's trim down our list of MRR sites significantly by filtering for only sites where, on average, there 
  # are 10 or greater detections for either Chinook salmon or steelhead (among only years where fish were detected
  # i.e., na.rm = T)
  filter(Chinook >= 10 | Steelhead >= 10) %>%
  select(-Chinook, -Steelhead) %>%
  # remove Lower Granite Dam MRR sites
  filter(!str_detect(site_code, "LGR")) %>%
  # remove acclimation ponds 
  filter(!str_detect(site_description, "AcclimationPond")) %>%
  # now, remove some additional sites that pass the above filter, but are not necessarily useful for escapement analyses.
  # this can be because we can't safely assign detections to another INT or MRR locations, or that carcass recoveries maybe
  # occurred prior to an adjacent array and detections muddle interpretation of escapement estimates there
  filter(!site_code %in% c("ASOTIC",  # Asotin Creek
                           "CATHEC",  # Catherine Creek
                           "CLWH",    # Clearwater Hatchery
                           "OXBO",    # Oxbow Hatchery
                           "GRAND2",  # Grande Ronde River - Wallowa River to headwaters (km 131-325)
                           "PANTHC",  # Panther Creek carcass recoveries (some occurred prior to PCA being installed)
                           "FISHC",   # Fish Creek carcass recoveries (occurred prior to LRU being installed)
                           "FISTRP",  # Fish Creek trap (occured prior to LRU being installed)
                           "SALR4",   # Salmon River - Pahsimeroi River to headwaters (km 489-650)
                           "SALRSF",  # South Fork Salmon River
                           "SNAKE4"   # Snake River - Salmon River to Hells Canyon Dam (km 303-397)
                           )) %>%
  # finally, add back in some MRR sites that are useful i.e., that they can easily be grouped to an adjacent INT site or
  # have sufficient, regular detections that they can be used to estimate a detection probability at a downstream site, etc.
  bind_rows(ptagis_sf %>%
              filter(site_code %in% c(# UPPER SALMON
                                      "ALTULC",   # Aluturas Lake Creek, Upper Salmon
                                      "BEAVEC",   # Beaver Creek, Upper Salmon
                                      "HAYDNC",   # Hayden Creek, Lemhi River
                                      "CARMEC",   # Carmen Creek - tributary to Salmon River
                                      "YANKWF",   # West Fork Yankee Fork Salmon River 
                                      # MIDDLE FORK SALMON
                                      "BEAV4C",   # Beaver Creek, Big Creek
                                      "BIG2C",    # Big Creek, Middle Fork Salmon River
                                      # SOUTH FORK SALMON
                                      "BURNLC",   # Burntlog Creek, Johnson Creek
                                      "GROUSC",   # Grouse Creek, Secesh River
                                      "LAKEC",    # Lake Creek, Secesh River
                                      "MCCA",     # McCall Hatchery
                                      "SUMITC",   # Summit Creek, Secesh River
                                      # LITTLE SALMON
                                      "RAPH",     # Rapid River Hatchery
                                      "RPDTRP",   # Rapid River Smolt Trap
                                      # SOUTH FORK CLEARWATER
                                      "KOOS",     # Kooskia National Fish Hatchery
                                      # POTLATCH RIVER
                                      "EFPW",     # East Fork Potlatch River weir
                                      # IMNAHA RIVER
                                      "CAMP4C",   # Camp Creek, tributary to Big Sheep Creek, Imnaha drainage
                                      "FREEZC",   # Freezeout Creek, Imnaha River
                                      "LSHEEF",   # Little Sheep Facility
                                      "MAHOGC",   # Mahogany Creek, Imnaha River
                                      # GRANDE RONDE
                                      "BCANF",    # Big Canyon Facility
                                      "WALH",     # Wallowa Hatchery
                                      # LOWER SNAKE
                                      "ALMOTC",   # Almota Creek - tributary to Snake River
                                      "CHARLC",   # Charley Creek, Asotin Creek watershed
                                      "PENAWC",   # Penawawa Creek - tributary to Snake River
                                      "TENMC2"))) # Tenmile Creek, tributary to Snake River

# MRR sites that have been included in past DABOM models that I've excluded:  
# "DWOR"   - Dworshak National Fish Hatchery
# "LYFE"   - Lyons Ferry Hatchery
# "NPTH"   - Nez Perce Tribal Hatchery

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
      site_code %in% c("SAWT", "BEAVEC", 
                       "ALTULC")          ~ "STL",    # Do I need to add SAWTRP??? Group Sawtooth Hatchery/Ladder Array & Upper Salmon carcass recoveries all to STL
      site_code %in% c("CEY", "YANKFK",
                       "YANKWF")          ~ "YFK_U",  # Group Yankee Fork and Cearley Creek obs to YFK_U
      site_code == "SALREF"               ~ "SALEFT", # Group EF Salmon River obs (e.g., carcass recoveries) w/ trap
      site_code == "HAYDNC"               ~ "HYC_U",  # Group Hayden Creek carcass recoveries to HYC_U
      site_code == "CARMEC"               ~ "CRC_U",  # Carmen Creek weir
      # MIDDLE FORK SALMON
      site_code %in% c("BIG2C", "BEAV4C") ~ "TAY_U",  # Group Big and Beavers creek obs w/ TAY_U
      # SOUTH FORK SALMON
      site_code %in% c("KNOXB","MCCA", 
                       "SALSFW")          ~ "STR",    # South Fork Salmon River weir
      site_code %in% c("SECESR", "GROUSC",
                       "SUMITC", "LAKEC") ~ "ZEN_U",  # Group Secesh River obs (e.g., carcass recoveries) w/ ZEN
      site_code == "BURNLC"               ~ "JOHNSC", # Burntlog Creek to Johnson Creek
      # LITTLE SALMON
      site_code == "RPDTRP"               ~ "RAPH",   # Group Rapid trap with Rapid Hatchery
      # SOUTH FORK CLEARWATER
      site_code == "KOOS"                 ~ "CLC",    # Group Kooskia Hatchery to CLC
      # IMNAHA RIVER
      site_code == "IMNAHW"               ~ "IML_U",   # Group Imnaha Weir w/ trap array
      site_code == "CAMP4C"               ~ "CMP_U",   # Group Camp Creek to CMP_U
      site_code %in% c("FREEZC", "MAHOGC")~ "IR3_U",
      # GRANDE RONDE
      site_code == "GRANDW"               ~ "UGS_U",   # Group GRANDW to UGS_U (based on site configuration)
      site_code == "CATHEW"               ~ "CCW_U",   # Group CATHEW to CCW_U (based on site configuration)
      site_code %in% c("LOOKGC", "LOOH")  ~ "LGW_U",   # Group LGW, LOOKGC, and LOOH into a single node
      # LOWER SNAKE
      site_code == "TUCH"                 ~ "TFH_U",   # Group Tucannon Hatchery to TFH_U
      site_code == "CHARLC"               ~ "CCA_U",   # Group Charley Creek to CCA_U
      site_code == "PENAWC"               ~ "PWA_U",   # Group Penawawa Creek to PWA_U
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
sites_sf = configuration %>%
  select(site_code = node) %>%
  mutate(site_code = str_replace(site_code, "_D$|_U$", "")) %>%
  distinct() %>%
  left_join(ptagis_sf) %>%
  # remove an extra Prosser Dam record
  filter(!(site_code == "PRO" & site_type == "MRR")) %>%
  st_as_sf(crs = 4326) %>%
  st_transform(crs = 5070) # NAD83

#----------------------
# download the NHDPlus v2 Flowlines
dwn_flw = T # do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
nhd_list = queryFlowlines(sites_sf = sites_sf,
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
parent_child = sites_sf %>%
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
                         # Potlatch River
                         c("HLM", "EPR", "EFPW"),
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
      "GRS", "ALMOTC", # Almota Creek
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

# add nodes to parent-child table (I don't think I need to create this here!)
pc_nodes = parent_child %>%
  select(parent, child) %>%
  addParentChildNodes(.,  configuration = configuration)  

# build paths to each node
node_paths = buildNodeOrder(parent_child = pc_nodes,
                            direction = "u")

#----------------------
# write configuration, parent-child table, flowlines, etc.
save(configuration,
     sites_sf,
     flowlines,
     parent_child,
     pc_nodes,
     node_paths,
     file = here("data/configuration_files/site_config_LGR_20240304.rda"))

# write sites_sf and flowlines out to geopackage, if desired
st_write(sites_sf, dsn = "data/spatial/dabom_sites.gpkg", layer = "sites_sf", driver = "GPKG", append = F)
st_write(flowlines, dsn = "data/spatial/dabom_sites.gpkg", layer = "flowlines", driver = "GPKG", append = F)

# END SCRIPT