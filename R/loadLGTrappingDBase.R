#------------------------------------------------------------------------------
# Function to connect to LGTrapping Database
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------

loadLGTrappingDBase <- function(trapDB_filepath)  {

  # check if filepath exists before connecting

  if (!file.exists(trapDB_filepath)) {
    stop("Trapping database file does not exist at ", trap_filepath,".")
  }

  # build connection strings

  path_string <- paste0("DBQ=", trapDB_filepath)
  driver_string <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
  connection_string <- paste0(driver_string, path_string)

  con <- DBI::dbConnect(odbc::odbc(),
                      .connection_string = connection_string)
  return(con)
}