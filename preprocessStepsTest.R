library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(optparse)

ldaModel <- readRDS("data/ldaModel.rds")

filetoparse <- "data/exampleData.csv"
normalisationPanel <- "data/standardisationPanel.tsv"

#Take the name of the dataset specifically
name <- gsub(".csv", "", basename(filetoparse))

#Import input datafiles
subselectedDataset <- tibble(fread(filetoparse))

#Read-in the panel
panel <- fread("data/standardisationPanel.tsv")

#Define function to scale a continuous variable according to mean and SD values indicated in a reference panel
#x is a numeric vector
#ID is an identifier for the row against which to scale
#panel defines the reference panel, with at least the columns: "ID",  "mean",  "sd"
scaleToPanel <- function(x, ID, panel) {
  rowID <- panel$ID == ID        #define the panel row to scale against
  x <- x - panel$mean[rowID]    #center on mean
  x <- x / panel$sd[rowID]      #Scale to reference variance
  return(x)
}

if ("Age_at_onset_years" %in% colnames(subselectedDataset)) {
  subselectedDataset$AGEONS_nml <-
    scaleToPanel(subselectedDataset$Age_at_onset_years, "Age_at_onset_years", panel)
}
if ("Diagnostic_delay_years" %in% colnames(subselectedDataset)) {
  subselectedDataset$DELAY_nml <-
    scaleToPanel(subselectedDataset$Diagnostic_delay_years,
                 "Diagnostic_delay_years",
                 panel)
}

if ("Time_to_death_or_last_followup_years" %in% colnames(subselectedDataset)) {
  subselectedDataset$SRV_nml <-
    scaleToPanel(
      subselectedDataset$Time_to_death_or_last_followup_years,
      "Time_to_death_or_last_followup_years",
      panel
    )
}

columnNames <- colnames(subselectedDataset)
columnNames[columnNames %in% c("Phenotype_1", "Phenotype_2", "Phenotype_3")] <- c("PHE_1", "PHE_2", "PHE_3")

colnames(subselectedDataset) <- columnNames

predicted <- predict(ldaModel, subselectedDataset)
