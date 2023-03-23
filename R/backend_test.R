# libraries
library(shiny)
library(ggplot2)
library(shinydashboard)
library(knitr)
library(rmarkdown)
library(ggthemes)
library(plotly)
library(DT)
library(data.table)
library(dplyr)
library(shinythemes)
library(RColorBrewer)
library(markdown)
library(shinycssloaders)
library(MASS)
library(stringr)
library(ggpubr)
library(rlang)
library(shinybusy)
library(caret)
library(xgboost)


clinicalData <- fread("data/exampleData.csv") %>% as.data.frame()

head(clinicalData)

mlModel <- readRDS("data/data_pheno/XGBoost_AUC_maybeFinal/results/XGBfinalTune.Rds")

clinicalData$class <- predict(mlModel, newdata = clinicalData)

head(clinicalData)
