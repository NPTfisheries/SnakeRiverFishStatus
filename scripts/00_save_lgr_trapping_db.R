# -----------------------
# Author(s): Ryan Kinzer and Mike Ackerman
# Purpose: Open connection to the Lower Granite Dam trapping database (LGTrappingDB)
#   and save table as a .csv
# 
# Created Date: Unknown
#   Last Modified: February 24, 2024
#
# Notes:

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(here)

# Requires a copy of the Lower Granite Dam trapping database and the ODBC driver
if(.Platform$OS.type != "unix") { # includes windows machines
  
  # source LGTrappingDB function
  source(here("R/loadLGTrappingDBase.R"))
  
  # set path to LGTrappingDB Access database
  trap_filepath = here("data/LGTrappingDB/LGTrappingExport_2024-02-21.accdb")
  
  # connect to LGTrappingDB
  con = loadLGTrappingDBase(trapDB_filepath = trap_filepath)
  
  # read master table within LGTrappingDB
  trap_dbase = DBI::dbReadTable(con, 'tblLGDMasterCombineExport')
  
  # and disconnect
  DBI::dbDisconnect(con)
  
} # end if not unix

# write .csv of LGTrappingDB for later use
write_csv(trap_dbase, file = paste0(here("data/LGTrappingDB/LGTrappingDB_/"), Sys.Date(), ".csv"))

# END SCRIPT