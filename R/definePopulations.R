#' @title Define Population Groups - LGR
#'
#' @description Define which main branches or detection sites are summed for TRT populations.
#'
#' @author Ryan Kinzer and Mike Ackerman
#'
#' @param spp Species, either "Chinook" or "Steelhead"
#'
#' @import dplyr stringr
#' @export
#' @return NULL
#' @examples definePopulations()
 
definePopulations = function(spc = c("Chinook", "Steelhead")) {
  
  spc = match.arg(spc)
  
  if(spc == "Chinook") {
    
    report_df = list(
      # Upper Salmon River MPG
      "SRUMA" = c("RFL", "STL"),
      "SRVAL" = "VC2",
      "SRYFS" = "YFK",
      "SREFS" = "SALEFT",
      "SRLMA" = c("USE_bb", "USI_bb"),
      "SRPAH" = "PAHH",
      "SRLEM" = c("CRC", "LLR"),
      "SRNFS" = "NFS",
      "SRPAN" = "PCA",
      # Middle Fork Salmon River MPG
      "MFBEA" = "BRC",
      "MFMAR" = "MAR",
      #"MFSUL" = NA,
      #"MFUMA" = NA,
      #"MFLOO" = NA,
      #"MFLMA" = NA,
      "MFBIG" = "TAY",
      #"MFCAM" = NA,
      #"SRCHA" = NA,
      # South Fork Salmon River MPG
      "SFSMA" = c("SFG_bb", "KRS"),
      "SFEFS" = "ESS",
      "SFSEC" = "ZEN",
      "SRLSR" = c("WB1", "RAPH"),
      # Wet Clearwater MPG
      "SEUMA/SEMEA/SEMOO" = "SW1",
      "CRLOC" = "LRL",
      "CRLOL" = "LC1",
      #"NCUMA" = NA,
      #"NCLMA" = NA,
      # Dry Clearwater MPG
      "SCUMA" = "SC1",
      "SCLAW" = c("SIX", "LAW", "CLC"),  
      "CRLAP" = c("LAP", "JA1"),           # add BED in SY2024
      "CRPOT" = "JUL",
      # Grande Ronde / Imnaha MPG
      "IRMAI" = c("COC", "IR1_bb", "IR2_bb", "IR3"), # COC is not operated for Chinook, but keeping here for record keeping and bc it is all 0s
      "IRBSH" = c("CMP", "LSC", "BSC"), # CMP and LSC are not operated for Chinook, but keeping here for record keeping and because they are all 0s
      "GRUMA" = "UGS",
      "GRCAT" = "CCU",
      "GRLOS" = c("WR1_bb", "WR2"),
      "GRMIN" = "MR1",
      "GRLOS/GRMIN" = "WR1",
      "GRLOO" = "LGW",
      "GRWEN" = "WEN",
      "Joseph"= "JOC",
      # Lower Snake MPG
      "SNASO" = "ACM",
      "SNTUC" = "LTR"
      ) %>%
      stack() %>%
      as_tibble() %>%
      select(TRT_POPID = ind,
             site = values)
  } # end if Chinook
  
  if(spc == "Steelhead") {
    report_df = list(
      # Salmon River MPG
      "SRUMA-s" = c("YFK", "VC2", "RFL", "STL"),
      "SREFS-s" = c("USI_bb", "SALEFT"),
      "SRPAH-s" = c("USE_bb", "PAHH"),
      "SRLEM-s" = c("CRC", "LLR"),
      "SRNFS-s" = "NFS",
      "SRPAN-s" = "PCA",
      #"SRCHA-s" = NA,
      "MFUMA-s" = c("BRC", "MAR"), # BRC is only operated for Chinook, but keeping it here for record-keeping and all estimates to BRC are 0.
      "MFBIG-s" = "TAY",
      "SFMAI-s" = c("SFG_bb", "KRS", "ESS"),
      "SFSEC-s" = "ZEN",
      "SRLSR-s" = c("WB1", "RAPH"),
      # Clearwater River MPG
      "CRSFC-s" = "SC1",
      "CRSEL-s" = "SW1",
      "CRLOC-s" = "LRL",
      "CRLOL-s" = "LC1",
      #"CRNFC-s" = NA,
      "CRLMA-s" = c("JUL", "LAP", "JA1", "SIX", "LAW", "CLC"), # add BED in SY2024
      # Hells Canyon MPG
      #"SNHCT-s" = NA,
      # Imnaha River MPG
      "IRMAI-s" = c("IR1", "COC"),
      # Grande Ronde River MPG
      "GRUMA-s" = c("UGR", "LGW"), 
      "GRWAL-s" = "WR1",
      "GRLMT-s" = "WEN",
      "GRJOS-s" = "JOC",
      # Lower Snake MPG
      "SNASO-s" = c("ACM", "ALMOTC", "ALPOWC", "TENMC2"),
      "SNTUC-s" = c("LTR", "PWA")
      ) %>%
      stack() %>%
      as_tibble() %>%
      select(TRT_POPID = ind,
             site = values)
  } # end if Steelhead
  
  return(report_df)
  
} # end definePopulations()
