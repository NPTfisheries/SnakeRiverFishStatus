# -----------------------
# Author(s): Mike Ackerman
# Purpose: Explore Lower Granite Dam tailrace (LGRTAL) detections to determine whether these
#   detections need to be handled differently in complete tag histories.
# 
# Created Date: June 16, 2026
#   Last Modified:
#
#   Notes:

# clear environment
rm(list = ls())

# load necessary libraries
library(tidyverse)
library(janitor)

# read in complete tag histories since SY2010
cth_df = list.files(path = "data/complete_tag_histories/",
                    pattern = "\\.csv$",
                    full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~read_csv(.x, show_col_types = F), .id = "file_name") %>%
  clean_names() %>%
  # add species and spawn year
  mutate(file_name = str_replace(file_name, ".*/", ""), 
         species = str_extract(file_name, "(?<=_)[^_]+"),       
         spawn_year = str_extract(file_name, "SY[0-9]{4}"),
         event_site_code_value = if_else(event_site_code_value == "3BV", "BV3", event_site_code_value)) %>%
  select(-file_name)

# return all records for any tag_code with a LGRTAL event or release site
lgrtal_cth_df = cth_df %>%
  group_by(tag_code) %>%
  filter(any(event_site_code_value == "LGRTAL" | event_release_site_code_code == "LGRTAL")) %>%
  ungroup()

# summarize how LGRTAL observations are recorded for records where LGRTAL appears in either event or release site code
lgrtal_code_summ = lgrtal_cth_df %>%
  filter(
    event_site_code_value == "LGRTAL" |
      event_release_site_code_code == "LGRTAL"
  ) %>%
  count(
    species,
    spawn_year,
    event_site_code_value,
    event_type_name,
    event_release_site_code_code,
    name = "n_records"
  ) %>%
  arrange(species, spawn_year)

# for tag_codes with any LGRTAL record, summarize release-site codes for events with LGR event_site_code_value
lgrtal_release_summ = lgrtal_cth_df %>%
  filter(event_site_code_value == "LGR") %>%
  count(
    species,
    spawn_year,
    event_type_name,
    event_release_site_code_code,
    name = "n_records"
  ) %>%
  arrange(species, spawn_year)

# for all LGR records, what various release site codes are being used
lgr_release_summ = cth_df %>%
  filter(event_site_code_value == "LGR") %>%
  count(
    species,
    spawn_year,
    event_site_code_value,
    event_type_name,
    event_release_site_code_code,
    name = "n_records"
  )

# read in processed observations since 2010
proc_df = list.files(path = "output/PITcleanr/human_reviewed/",
                     pattern = "\\.xlsx$",
                     full.names = T) %>%
  setNames(nm = .) %>%
  map_df(~ read_xlsx(.x), .id = "file_name") %>%
  # add species and spawn year
  mutate(
    file_name = basename(file_name),
    species = str_extract(file_name, "^[^_]+"),
    spawn_year = str_extract(file_name, "SY\\d{4}")
  ) %>%
  select(-file_name)

lgrtal_misfits = cth_df %>%
  group_by(tag_code) %>%
  filter(
    any(
      event_site_code_value == "LGR" & event_release_site_code_code %in% c("LGRRBR", "LGRRRR", "LGRTAL", "LGRRTR", "SNAKE2")
    )
  )

### END SCRIPT
