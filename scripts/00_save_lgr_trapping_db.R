# -----------------------
# Author(s): Ryan Kinzer and Mike Ackerman
# Purpose: Open connection to the Lower Granite Dam trapping database (LGTrappingDB)
#   and save master table to a .csv
# 
# Created Date: Unknown
#   Last Modified: January 6, 2026
#
# Notes: due to file size, the Lower Granite trapping database is saved to a .gitignore directory i.e., the user needs to download their own local copy

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# requires a copy of the Lower Granite Dam trapping database and the ODBC driver
if(.Platform$OS.type != "unix") {                                       # includes windows machines
  
  trap_filepath = "data/LGTrappingDB/LGTrappingExport_2026-01-06.accdb" # set path to LGTrappingDB Access database
  source("R/loadLGTrappingDBase.R")                                     # source loadLGTrappingDBase function
  con = loadLGTrappingDBase(trapDB_filepath = trap_filepath)            # connect to database
  trap_dbase = DBI::dbReadTable(con, 'tblLGDMasterCombineExport')       # read master table within LGTrappingDB
  DBI::dbDisconnect(con)                                                # disconnect
}

# write database to .csv
write_csv(trap_dbase, file = paste0(here("data/LGTrappingDB/LGTrappingDB_/"), Sys.Date(), ".csv"))

# END SCRIPT