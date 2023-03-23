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
library(MASS)
library(ggplot2)
library(tidyr)
library(dplyr)
library(moments) #For checking skewness
library(data.table)
library(writexl)
library(reshape2)
library(psych)
library(optparse)

opt <- NULL
opt$dataset2 <- NULL
keepcols <- c("C", "CTYDEL", "AGEONS_nml", "Site_of_Onset", "Sex_at_birth", "Phenotype", "origin")
predictorcols <- c("CTYDEL","AGEONS_nml","Site_of_Onset","Sex_at_birth","Phenotype")
opt$intPredict <- TRUE
path <- "data"

data <- data.table::fread(file="data/saved_mplusfit.tsv"#,select=keepcols, header=TRUE
                          )
data$C <- as.factor(data$C)
data <- na.omit(data)

#Dummy code Phenotype
if("Phenotype" %in% names(data) && "Phenotype" %in% predictorcols){
  data <- data %>%
    mutate(row=row_number()) %>% #Row is a unique identifier
    pivot_wider(names_from = Phenotype, values_from = Phenotype, names_prefix= "PHE_",values_fill=0,names_sort=TRUE) %>%
    mutate(across(starts_with("PHE_"), ~ if_else(.>=2,1, 0))) %>%
    dplyr::select(-c(PHE_1,row))

  predictorcols <- c(predictorcols[!grepl("Phenotype",predictorcols)], "PHE_2","PHE_3")
}

if("origin" %in% names(data)){
  data$origin <- factor(data$origin,
                        levels=c("BE","CH","ES","FR","GB","IE","IL","IT","NL","PT","SE","TR","US"),
                        labels= c("Belgium","Switzerland","Spain","France",
                                  "United Kingdom","Ireland","Israel","Italy","Netherlands",
                                  "Portugal","Sweden","Turkey","USA"))
}


accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

runLDA <- function(data,predictorcols,testdata=NULL,path=NULL){

  f_rh <- reformulate(predictorcols)
  lda_fit  <- lda(update.formula(f_rh,C~.), data=data,CV=FALSE)

  #If the main feature of LD1 (diag delay) is negative coef, then flip the valence of LD1
  if("CTYDEL" %in% rownames(lda_fit$scaling)){
    if(lda_fit$scaling["CTYDEL",1]<0){
      lda_fit$scaling[,1] <- lda_fit$scaling[,1]*-1
    }
  }

  #If the main feature of LD2 (disease duration) is negative coef, then flip the valence of LD2
  if("SRV_nml" %in% rownames(lda_fit$scaling)){
    if(lda_fit$scaling["SRV_nml",2]<0){
      lda_fit$scaling[,2] <- lda_fit$scaling[,2]*-1
    }
  }

  #Save eigenvalues
  eigen <- lda_fit$svd

  #If testdata argument is provided, predict in the new dataset. Otherwise predict in the original
  if(!is.null(testdata)){
    predicted <- predict(lda_fit,testdata) #Predict with 'newdata'
    whichdata <- "LDA_externalaxisvals_" #Indicate part of file name for writing the prediction data
    testdata$C <- droplevels(testdata$C)
    data <- testdata
  } else {
    predicted <- predict(lda_fit)
    whichdata <- "LDA_axisvals_"
  }

  #create confusion matrix, test accuracy and reformat for writing to file
  tabulate <- table(predicted$class,data$C)
  perc_correct <- accuracy(tabulate)
  tabulate <- matrix(tabulate,nrow=nrow(tabulate))
  rownames(tabulate) <- paste0("Pred_",1:nrow(tabulate))

  # #Adjust formatting of posterior probs so that they can be summarised
  post <- data.frame(predicted$posterior)
  names(post) <- paste0("C",dimnames(predicted$posterior)[[2]])

  post<- cbind(C=data$C, post)

  # #Generate table of average posterior probabilities
  post_prob <- post %>%
    group_by(C) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

  #Save axis loadings and predicted vs actual
  axis<- data.frame(C=data$C,predictedC=predicted$class, predicted$x)
  data.table::fwrite(axis,file=paste0(path,whichdata,"predictions.csv"))

  #Define data for plotting, combining data and predicted axis values
  plot_dat <- data.table(cbind(data, predicted$x))

  #Manual palette: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  cbbPalette <- c("gold", "violetred", "#56B4E9", "#000000", "limegreen", "#D55E00", "#CC79A7", "#F0E442")

  #Syntax to obtain centroids
  #https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot
  ##centroids <- aggregate(cbind(LD1,LD2)~C,plot_dat,mean)

  #Construct plot across the first two linear discriminant axes
  lda_plot <- ggplot(plot_dat, aes(LD1, LD2)) +
    scale_colour_manual(values=cbbPalette,name = "Class",labels=c(paste0(levels(plot_dat$C), " [n=",unname(table(plot_dat$C)),"]")))+
    theme_bw()+
    theme(legend.position = 'right')+
    xlab("Linear Discriminant 1 (LD1)")+
    ylab("Linear Discriminant 2 (LD2)")

  #Depreciated option to differentiate phenotypes by shape
  # if("Phenotype" %in% names(plot_dat) && length(unique(plot_dat$Phenotype))>1){
  #   lda_plot <- lda_plot +
  #     geom_point(aes(color = C,shape=Phenotype),size=0.7)+
  #     scale_shape_manual(name="Phenotype",values=c(20,6,4))
  # } else {
  lda_plot <- lda_plot +
    geom_point(aes(color = C),size=0.7)
  #}

  #Save the figure
  ggsave(plot=lda_plot, filename = paste0(path,"lda_plot_",Sys.Date(),".pdf"),
         units="mm",width=175,height=150)

  #If origin column is in the dataset, produce additional plot with facets by origin
  if("origin" %in% names(plot_dat) && length(unique(plot_dat$origin))>1){
    lda_plot_origin <- lda_plot+
      facet_wrap(~origin,ncol=3)+
      guides(color=guide_legend(nrow=2, byrow=TRUE))+
      theme(legend.position = c(1, 0),
            legend.justification = c(1.1,-0.2),
      )
    theme(legend.position = 'top')

    #Save the figure
    ggsave(plot=lda_plot_origin,
           filename = paste0(path,"lda_plot_origin_",Sys.Date(),".pdf"),
           units="mm",width=150,height=200)

  }

  ### Run pooled correlations, then compare correlations between variables and LD axes
  ld_corr <- psych::statsBy(plot_dat[,c("C",predictorcols,grep("LD",colnames(plot_dat),value=TRUE)),with=FALSE],group="C",cors=TRUE)
  LDs <- grep("LD",colnames(ld_corr$rwg))
  pooled_cor<- ld_corr$rwg[-LDs,LDs]
  if(length(LDs)==1) {
    cors<- data.frame(variable=names(pooled_cor),pooled_cor)
  } else {
    cors<- data.frame(variable=rownames(pooled_cor),pooled_cor)
  }

  #Write list summary into xlsx file
  summary <- list(coefs = data.frame(variable=rownames(lda_fit$scaling),lda_fit$scaling),
                  eigen = data.frame(eigenvalues=lda_fit$svd,proptrace=prop.table(lda_fit$svd^2)),
                  prediction = as.data.frame(perc_correct),
                  predic_tab = data.frame(cbind(rows=rownames(tabulate),tabulate)),
                  post_prob = post_prob,
                  corrs = cors
  )

  writexl::write_xlsx(x=summary,
                      path=paste0(path,"LDA_results_",Sys.Date(),".xlsx"),
                      col_names = TRUE)

  saveRDS(lda_fit, file = paste0(path,"LDA_results_",Sys.Date(),".rds"))
}



if(is.null(opt$dataset2)){

  #Run LDA on the full training data
  runLDA(data = data,predictorcols = predictorcols, path = path)

  #If using internal validation, split the dataset into train and test samples
  if(opt$intPredict==TRUE){
    #Obtain directory tree minus model name. This is where output files will be saved
    INTpredpath <- paste0(path,"LDApredictINT/")
    if(!dir.exists(INTpredpath)){dir.create(INTpredpath,recursive = TRUE)}

    #Assign test and train values in 8:2 proportion
    INTdata <- data

    #Randomly assign to training and sample dataset
    set.seed(84)
    assigngroup <- sample(c("TRAIN", "TEST"), nrow(INTdata), replace=TRUE, prob=c(0.8,0.2))
    internaltraindata<- INTdata[assigngroup=="TRAIN",]
    internaltestdata <- INTdata[assigngroup=="TEST",]

    runLDA(data=internaltraindata,predictorcols = predictorcols,testdata=internaltestdata,path=INTpredpath)

  }
}
