# -----------------------
# Author(s): Ryan Kinzer and Mike Ackerman
# Purpose: Open connection to the Lower Granite Dam trapping database (LGTrappingDB)
#   and save master table as a .csv
# 
# Created Date: Unknown
#   Last Modified: May 21, 2025
#
# Notes:

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(here)

# Requires a copy of the Lower Granite Dam trapping database and the ODBC driver
if(.Platform$OS.type != "unix") { # includes windows machines
  
  # set path to LGTrappingDB Access database
  trap_filepath = here("data/LGTrappingDB/LGTrappingExport_2025-05-21.accdb")
  
  # source loadLGTrappingDBase function
  source(here("R/loadLGTrappingDBase.R"))
  
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