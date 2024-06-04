library('Hmisc')
library('GGally')
library("forestmodel")
library("survival")
library(plyr)
library("dplyr")
library("ggpubr")
library("lattice")
library("Formula")
library(ggplot2)
library(magrittr)
require("survival")
library("survminer")
library(rms)
library(Gmisc)
library(Greg)
library(splines)
library(smoothHR)
# library("svglite")
# library("systemfonts")
library(Kendall)
# library(ltm)
library(devtools)



# --------------------- Random Forest ---------------------- #
random_forest <-function(){
  model <- coxph( Surv(time, event) ~ .,
                  data = dataset )
  ggforest(model)
}
# --------------------- Linear Regression ---------------------- #

linear_regression <-function(dataset, feature_one_list, feature_two_list,x_label, y_label,file_name="",save_plot=FALSE){
  model <- lm(feature_one_list~feature_two_list, data=dataset)
  print(summary(model))


  #create plot to visualize fitted linear regression model
  p <- ggplot(dataset,aes(feature_one_list, feature_two_list)) + geom_point() + geom_smooth(method='lm')+labs(x = x_label, y =y_label)
  p

  if(save_plot){
      pdf(file = file_name)
      print(p,newpage = FALSE)
      dev.off()
  }

}


# --------------------- Box Plot ---------------------- #
box_plot<-function(dataset, feature_one, feature_two,x_label, y_label,file_name="",save_plot=FALSE){


  p <- boxplot(feature_one~feature_two,# here you write the feature column names in the csv file
  data=dataset,
  xlab=x_label,
  ylab=y_label,
  col="gray",
  border="black"
  )
  p

  if(save_plot){
    pdf(file = file_name)
    print(p,newpage = FALSE)
    dev.off()
  }
}
# --------------------- Correlation ---------------------- #
correlation_method <-function(dataset, feature_one, feature_two,x_label, y_label, corr_method="spearman",file_name = '',save_plot=FALSE){


  p <- ggscatter(dataset, y = feature_two, x = feature_one,
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE, cor.method = corr_method,cor.coef.size = 6)#, xlim=c(5, 75))
  p <- p+theme(axis.text = element_text(size = 14, color = "black", face = "bold"), axis.title.x = element_text(size = 14, color = "black", face = "bold"),legend.text = element_text(size = 14, color = "black", face = "bold"),
                                                     axis.title.y = element_text(size = 14, color = "black", face = "bold"))
  p <- p+labs(x = x_label, y =y_label)

  print(p)


}
#--------------------- Binary Correlation --------------
binary_correlation<-function(feature_one, feature_two){
    cor.test(x=feature_one,y=feature_two, method="kendall")
    biserial.cor(x=feature_one,y=feature_two)
    summary(Kendall(feature_one,feature_two))
}
#--------------------- Forest Plot ---------------------
forest_plot<-function(dataset,save_plot=FALSE,file_name=""){
  #to factor the feature with category as numbers
  dataset <- within(dataset, {
      tumour_size_cat <- factor(tumour_size_cat)
      grade <- factor(grade)
      axillary_nodes <- factor(axillary_nodes)
    #differ <- factor(differ, labels = c("well", "moderate", "poor"))
    #extent <- factor(extent, labels = c("submuc.", "muscle", "serosa", "contig."))
  })


  model <- coxph( Surv(time, event) ~Digi_sTILs + TAS_scores_max + tumour_size_cat+ axillary_nodes+age,
                  data = dataset )
  ggforest(model)
  p <- ggforest(model)
  print(p)
  if (save_plot){
    forest_plot(model)
    ggsave(file_name, plot=p, width=7.0, height=3.0)}

}
# --------------------- Correlation (Loop over features or marker) ---------------------- #
multi_correlation<-function(save_csv=FALSE,show_plot=FALSE,method='spearman',file_name='',target_col=1,start_cols=1, end_cols=1){

  correlation_markers <- data.frame(matrix(ncol = 3, nrow = 0))

  for (n in seq(start_cols,end_cols, by=1)){

    col_one <- dataset[[n]]

    col_two <- dataset[[target_col]]

    try({
      res <- cor.test(col_one,col_two,
                         method = method,exact=FALSE)

      if (res$p.value < 0.05){

        print(n)
        print(colnames(dataset)[n])

        correlation_markers[nrow(correlation_markers) + 1,] <-c(colnames(dataset)[n],res$estimate,(res$p.value))
      }
      if (show_plot){

        feature_one <- colnames(dataset)[n]
        feature_two <- colnames(dataset)[target_col]

        ggscatter(dataset, x = feature_one, y = feature_two,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = TRUE, cor.method = "spearman",
                   xlab = col_one, ylab = col_two)
    }
    })
  }

  print(correlation_markers)
  print(colnames(dataset)[target_col])

  if(save_csv){
    write.csv(correlation_markers,paste(file_name,'.csv'))
  }
}
# --------------------- Correlation (Loop over features or marker with adjusted p-values) ---------------------- #
multi_correlation_adjusted_p_values<-function(save_csv=FALSE,adjusted_method ='bonferroni',method='spearman',file_name='',target_col=1,start_cols=1, end_cols=1){

  pvaluesList <- list()
  featureList <- list()
  rvaluesList <- list()
  for (n in seq(start_cols,end_cols, by=1)){

    col_one <- dataset[[n]]

    col_two <- dataset[[target_col]]

    featureList <- append(featureList,colnames(dataset)[n])

    try({
      res <- cor.test(col_one,col_two,
                         method = method,exact=FALSE)
      pvaluesList <-append(pvaluesList,res$p.value)

      rvaluesList <- append(rvaluesList,res$estimate)

    })
  }
  psAdjusted <- p.adjust(pvaluesList,method=adjusted_method)

  print(colnames(dataset)[target_col])
  mpFeature<-do.call(rbind, Map(data.frame, A=featureList, B=rvaluesList, C=psAdjusted))

  if(save_csv){
    write.csv(mpFeature,paste(file_name,'.csv'))
  }
}

rcs_plotHR<-function(dataset,save_plot=FALSE, file_name=""){


  fit.cph <- cph(Surv(time, event)  ~ rcs(Digi_sTILs,3)+score, data=dataset)
  anova(fit.cph)

  m <- cph(Surv(time, event)  ~ ., data=dataset)

  p.val <- 1- pchisq( (fit.cph$loglik[2]- m$loglik[2]), 2 )
  print(p.val)
  #summary(fit.cph)
  p <- plotHR(m, term=1, plot.bty="l", xlim=c(0, 0.1),col.term = "#DF76EE",col.se = "#0F5129",xlab = "Digi-sTILs + Age + Lymph Nodes")
  print(p)

  if (save_plot){
    pdf(file = file_name)   # The directory you want to save the file in
    print(p,newpage = FALSE)
    dev.off()
  }
}
#--------------------- Scattered Plots ---------------------#
scattered_plot<-function(dataset,feature_one,feature_two,feature_three,x_lab, y_lab){
  cols <- c("blue","orange","green")
  p<-ggplot(dataset,aes(feature_one,feature_two,colour=feature_three))+labs(x = x_lab, y = y_lab, color ="RCB Score")+geom_point()
  p<-p+scale_color_manual(values = cols)
  p
}