# functions for sensitivity, specificity, ppv, npv, accuracy and auc----

hmean<-function(x) {1/mean(1/x)}
sd.hmean<-function(x) {sqrt((mean(1/x))^(-4)*var(1/x)/length(x))}


fx_sensitivity <- function(tp, fn){
  sens <- as_tibble(binconf(tp,
                            (tp + fn),
                            alpha = 0.05,
                            method = "wilson"))
  colnames(sens) <- c("value", "lower", "upper")
  return(sens)
}


fx_specificity <- function(tn, fp){
  spec <- as_tibble(binconf(tn,
                            (tn + fp),
                            alpha = 0.05,
                            method = "wilson"))
  colnames(spec) <- c("value", "lower", "upper")
  return(spec)
}

fx_ppv <- function(tp, fp){
  ppv <- as_tibble(binconf(tp,
                           (tp + fp),
                           alpha = 0.05,
                           method = "wilson"))
  colnames(ppv) <- c("value", "lower", "upper")
  return(ppv)
}

fx_npv <- function(tn, fn){
  npv <- as_tibble(binconf(tn,
                           (tn + fn),
                           alpha = 0.05,
                           method = "wilson"))
  colnames(npv) <- c("value", "lower", "upper")
  return(npv)
}

fx_accuracy <- function(tp, fp, tn, fn){
  accuracy <- as_tibble(binconf((tp + tn),
                                (tp + fp + tn + fn),
                                alpha = 0.05,
                                method = "wilson"))
  colnames(accuracy) <- c("value", "lower", "upper")
  return(accuracy)
}


fx_ci.auc <- function(boot.roc){
  auc <- tibble(value = boot.roc$AUC)
  auc$lower <- boot_ci(boot.roc, AUC, in_bag = TRUE, alpha = 0.05)[1,][["values"]]
  auc$upper <- boot_ci(boot.roc, AUC, in_bag = TRUE, alpha = 0.05)[2,][["values"]]
  return(auc)
}




# function to calculate the roc curve, the cut point and a summary table with cutpointer ----
fx_roc_analysis <- function(data = NULL,
                   predictor.vars = NULL,
                   outcome.vars = NULL,
                   method = maximize_metric, # method to determine the cutpoints
                   metric = youden, # function for computing a metric when maximize_metric or minimize_metric is used
                   pos_class = NULL, # value indicatin the positive class
                   neg_class = NULL,
                   direction = NULL,
                   boot_runs = 1000, # number of runs to assess the variability
                   boot_stratify = FALSE,
                   use_midpoints = FALSE,
                   break_ties = mean,
                   na.rm = FALSE,
                   allowParallel = TRUE,
                   silent = FALSE,
                   tol_metric = 1e-06,
                   ...
){
  
  # first run
  
  i = predictor.vars[1]
  j = outcome.vars[1]
  
  data1 <- data[,c(i,j)]
  
  names(data1) <- c("predictor.var", "outcome.var")
  
  roc_temp <- cutpointr(data = data1, # Data
                        x = predictor.var, # predictor
                        class = outcome.var, # class membership
                        method = method, # method to determine the cutpoints
                        metric = metric, # function for computing a metric when maximize_metric or minimize_metric is used
                        pos_class = pos_class, # value indicatin the positive class
                        neg_class = neg_class,
                        direction = direction,
                        boot_runs = boot_runs, # number of runs to assess the variability
                        boot_stratify = boot_stratify,
                        use_midpoints = use_midpoints,
                        break_ties = break_ties,
                        na.rm = na.rm,
                        allowParallel = allowParallel,
                        silent = silent,
                        tol_metric = tol_metric
  )
  rm(data1)
  
  sum_temp <- summary(roc_temp)
  
  tp <- sum_temp[["confusion_matrix"]][[1]][["tp"]]
  fn <- sum_temp[["confusion_matrix"]][[1]][["fn"]]
  fp <- sum_temp[["confusion_matrix"]][[1]][["fp"]]
  tn <- sum_temp[["confusion_matrix"]][[1]][["tn"]]
  
  temp <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                tibble(outcome = j, predictor = i, metric = "auc", fx_ci.auc(roc_temp))
  )
  rm(list = c("tp","fp","tn","fn"))
  
  metrics <- NULL
  metrics <- rbind(metrics, temp)
  rm(temp)
  
  roc_temp$metrics<-list(metrics)
  rm(metrics)
  
  roc_temp[["outcome"]] <- j
  roc_temp[["predictor"]] <- i
  
  roc_temp[["data"]][[1]][["outcome"]] <- j
  roc_temp[["data"]][[1]][["predictor"]] <- i
  
  roc_temp[["roc_curve"]][[1]][["outcome"]] <- j 
  roc_temp[["roc_curve"]][[1]][["predictor"]] <- i 
  
  names(roc_temp[["data"]])<-paste(j,i,sep = "_")
  names(roc_temp[["roc_curve"]])<-paste(j,i,sep = "_")
  names(roc_temp[["boot"]])<-paste(j,i,sep = "_")
  names(roc_temp[["metrics"]])<-paste(j,i,sep = "_")
  
  roc_analysis<-roc_temp
  rm(list = c("roc_temp","i","j"))
  
  if (length(predictor.vars) > 1){
    for (j in outcome.vars){
      for (i in predictor.vars[-1]){
        
        data1 <- data[,c(i,j)]
        
        names(data1) <- c("predictor.var", "outcome.var")
        
        roc_temp <- cutpointr(data = data1, # Data
                              x = predictor.var, # predictor
                              class = outcome.var, # class membership
                              method = method, # method to determine the cutpoints
                              metric = metric, # function for computing a metric when maximize_metric or minimize_metric is used
                              pos_class = pos_class, # value indicatin the positive class
                              neg_class = neg_class,
                              direction = direction,
                              boot_runs = boot_runs, # number of runs to assess the variability
                              boot_stratify = boot_stratify,
                              use_midpoints = use_midpoints,
                              break_ties = break_ties,
                              na.rm = na.rm,
                              allowParallel = allowParallel,
                              silent = silent,
                              tol_metric = tol_metric)
        rm(data1)
        
        sum_temp <- summary(roc_temp)
        
        tp <- sum_temp[["confusion_matrix"]][[1]][["tp"]]
        fn <- sum_temp[["confusion_matrix"]][[1]][["fn"]]
        fp <- sum_temp[["confusion_matrix"]][[1]][["fp"]]
        tn <- sum_temp[["confusion_matrix"]][[1]][["tn"]]
        
        temp <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                      tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                      tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                      tibble(outcome = j, predictor = i, metric = "auc", fx_ci.auc(roc_temp))
        )
        rm(list = c("tp","fp","tn","fn"))
        
        metrics <- NULL
        metrics <- rbind(metrics, temp)
        rm(temp)
        
        roc_temp$metrics<-list(metrics)
        rm(metrics)
        
        roc_temp[["outcome"]] <- j
        roc_temp[["predictor"]] <- i
        
        roc_temp[["data"]][[1]][["outcome"]] <- j
        roc_temp[["data"]][[1]][["predictor"]] <- i
        
        roc_temp[["roc_curve"]][[1]][["outcome"]] <- j 
        roc_temp[["roc_curve"]][[1]][["predictor"]] <- i 
        
        names(roc_temp[["data"]])<-paste(j,i,sep = "_")
        names(roc_temp[["roc_curve"]])<-paste(j,i,sep = "_")
        names(roc_temp[["boot"]])<-paste(j,i,sep = "_")
        names(roc_temp[["metrics"]])<-paste(j,i,sep = "_")
        
        roc_analysis <- rbind(roc_analysis,roc_temp)
        
      }
    }
    rm(list = c("roc_temp","i","j"))
  }
  
  metrics_comb <- NULL
  
  for (i in names(roc_analysis[["metrics"]])){
    temp<-roc_analysis[["metrics"]][[i]]
    metrics_comb<-rbind(metrics_comb,temp)
    rm(temp)
  }
  
  data_comb <- NULL
  
  for (i in names(roc_analysis[["data"]])){
    temp<-roc_analysis[["data"]][[i]]
    data_comb<-rbind(data_comb,temp)
    rm(temp)
  }
  
  
  rocc_comb <- NULL
  
  for (i in names(roc_analysis[["roc_curve"]])){
    temp<-roc_analysis[["roc_curve"]][[i]]
    rocc_comb<-rbind(rocc_comb,temp)
    rm(temp)
  }
  
  sens_points <- c(roc_analysis[["sensitivity"]], rep(NA, times = (nrow(rocc_comb) - nrow(roc_analysis))))
  spec_points <- c((1 - roc_analysis[["specificity"]]), rep(NA, times = (nrow(rocc_comb) - nrow(roc_analysis))))
  
  roc_plot <- ggplot(rocc_comb, aes(x=fpr, y=tpr, group=predictor)) +
    geom_line(aes(linetype = predictor, color=predictor), size = 1.25) +
    geom_point(aes(x = spec_points,
                   y = sens_points), size = 2) +
    labs(title = "ROC Curve",
         y = "Sensitivity",
         x = "1 - Specificity") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(axis.text = element_text(size = 11)) +
    theme(legend.position="bottom")
  
  metrics_plot <- ggplot(metrics_comb, aes(x = predictor, y = value, ymin = lower, ymax = upper)) +
    geom_pointrange() + 
    geom_hline(yintercept = 0.8, lty = 2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Predictor") + ylab("Mean (95% CI)") +
    theme_bw() +  # use a white background
    facet_wrap( ~ metric, ncol = 3)
  
  roc_analysis <- tibble(roc_analysis = list(roc_analysis),
                         data = list(data_comb),
                         roc_curve = list(rocc_comb),
                         metrics = list(metrics_comb),
                         roc_plot = list(roc_plot),
                         metrics_plot = list(metrics_plot)
                         )
  
  return(roc_analysis)
  
}


# function to calculate the roc curve and statisticals differences in AUC, the cut point (threshold) and a summary table with pROC ----

fx_roc_analysis2 <- function(response = NULL,
                             predictors = NULL,
                             data = NULL,
                             ci.method = "bootstrap",
                             boot.n = 10000,
                             boot.stratified = T,
                             best.method = "youden",
                             conf.level = .95,
                             x = "best",
                             input = NULL
){
  # Constructing ROC curve
  roc_analysis <- NULL
  temp2 <- NULL
  
  for(i in response){
    temp1 <- pROC::roc(as.formula(paste0(i, sep = "~", paste(predictors, collapse = "+"))),
                       data = data,
                       ci = T,
                       ci.method = ci.method,
                       boot.stratified = boot.stratified,
                       boot.n = boot.n)
    temp2[i] <- list(temp1)
    roc_analysis <- c(roc_analysis, temp2)
  }
  
  rm(list=c("temp1", "temp2"))
  # Calculating metrics pluss ci: Sens, Spec, ppv, npv, Auccuracy, threshold 
  metrics <-NULL
  for(j in levels(as.factor(names(roc_analysis)))){
    for(i in names(roc_analysis[[j]])){
      temp1 <- coords(roc_analysis[[j]][[i]],
                      x = x,
                      input = input,
                      best.method = best.method,
                      best.policy = "random",
                      ret = c("tp", "tn", "fp", "fn"),
                      transpose = F)
      
      tp<-temp1[["tp"]]
      fp<-temp1[["fp"]]
      tn<-temp1[["tn"]]
      fn<-temp1[["fn"]]
      
      temp2 <- ci.coords(roc_analysis[[j]][[i]],
                         x = "best",
                         ret = "threshold",
                         best.method = best.method,
                         best.weights = c(1, 0.5),
                         best.policy = "random",
                         boot.n = boot.n,
                         conf.level = conf.level,
                         transpose = T)
      
      temp3 <- rbind(tibble(outcome = j, predictor = i, metric = "sensitivity", fx_sensitivity(tp = tp, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "specificity", fx_specificity(tn = tn, fp = fp)),
                     tibble(outcome = j, predictor = i, metric = "ppv", fx_ppv(tp = tp, fp = fp)),
                     tibble(outcome = j, predictor = i, metric = "npv", fx_npv(tn = tn, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "accuracy", fx_accuracy(tp = tp, fp = fp, tn = tn, fn = fn)),
                     tibble(outcome = j, predictor = i, metric = "auc", value = roc_analysis[[j]][[i]][["ci"]][[2]],
                            lower = roc_analysis[[j]][[i]][["ci"]][[1]],
                            upper = roc_analysis[[j]][[i]][["ci"]][[3]]),
                     tibble(outcome = j, predictor = i, metric = "threshold", value = temp2[["threshold"]][[2]],
                            lower = temp2[["threshold"]][[1]],
                            upper = temp2[["threshold"]][[3]])
      )  
      
      metrics<-rbind(metrics, temp3)
      rm(list = c("temp1", "temp2", "temp3", "tp", "fp", "tn", "fn"))
    }}
  
  
  
  auc.test <- tibble(response = as.character(NULL),
                     comparison = as.character(NULL),
                     roc1 = as.character(NULL),
                     roc2 = as.character(NULL),
                     auc1 = as.double(NULL),
                     auc2 = as.double(NULL),
                     difference = as.double(NULL),
                     p.value = as.double(NULL),
                     power = as.double(NULL))
  
  for (h in levels(as.factor(names(roc_analysis)))){
    for(i in names(roc_analysis[[h]])){
      for(j in names(roc_analysis[[h]])){
        if (i != j & !(paste(j, i, sep = "-") %in% auc.test[auc.test[["response"]] == h,][["comparison"]])){
          temp <- tibble(response = h, comparison = paste(i, j, sep = "-"), roc1 = i, roc2 = j, auc1 = pROC::auc(roc_analysis[[h]][[i]]), auc2 = pROC::auc(roc_analysis[[h]][[j]]))
          temp$difference <- pROC::auc(roc_analysis[[h]][[i]]) - pROC::auc(roc_analysis[[h]][[j]])
          temp1 <- roc.test(roc_analysis[[h]][[i]],roc_analysis[[h]][[j]], method = ci.method)
          temp$p.value <- temp1$p.value
          temp2 <- power.roc.test(roc_analysis[[h]][[i]], roc_analysis[[h]][[j]], method = "delong")
          temp$power <- temp2$power
          auc.test <- rbind(auc.test, temp) 
        }
      }
    }  
  }
  
  
  rm(list = c("i", "j", "temp", "temp1", "temp2"))  
  
  
  roc_curves <- NULL
  
  for (h in levels(as.factor(names(roc_analysis)))){
    for (i in names(roc_analysis[[h]])){
      temp <- tibble(response = h,
                     predictor = i,
                     sensitivities = roc_analysis[[h]][[i]][["sensitivities"]],
                     specificities = roc_analysis[[h]][[i]][["specificities"]],
                     thresholds = roc_analysis[[h]][[i]][["thresholds"]])
      roc_curves <- rbind(roc_curves,temp)
      rm(temp)
    }  
  }
  
  roc_analysis <- list(roc_analysis = roc_analysis,
                       metrics = metrics,
                       roc_curves = roc_curves,
                       auc.test = auc.test)
  
  return(roc_analysis)
  
}


# function to perfomre the heirarchical clustering method ----

fx_HClust.Analysis <- function(data = NULL,
                               names.variables = NULL,
                               names.rows = NULL,
                               distance = "auto",
                               method = "auto",
                               NbC.index ="alllong"){
  
  # Step One: Data preparation
  data_df<-as.data.frame(data[names.variables])
  row.names(data_df)<-data[[names.rows]]
  data_matrix <- as.matrix(scale(data_df))
  
  # step two: Measure pairwise distances
  # Define the method: auto, euclidean, manhattan,
  # canberra, binary, minkowski
  
  # Step three: define the clustering algorithm (two possibilities)
  ## Agglomrative clustering (hclust)
  # and the method (hc.method) to compute dissimilarities
  #between clusters:
  #hclust
  # "ward.D", "ward.D2", "single",
  #"complete", "average" (= UPGMA),
  #"mcquitty" (= WPGMA),
  #"median" (= WPGMC) or "centroid" (= UPGMC).
  
  # Step four: define the number of clusters
  
  if(distance == "auto" & method == "auto"){
    d <- c("euclidean", "manhattan",
           "canberra", "binary", "minkowski")
    m <- c("ward.D", "ward.D2", "single",
           "complete", "average",
           "mcquitty")
    ac <- NULL
    for (i in d){
      temp1 <- get_dist(data_matrix, method = i)
      for(j in m){
        temp2 <- hclust(temp1, method = j)
        temp3 <- tibble(distance = i, method = j, ac = coef.hclust(temp2))
        ac <- rbind(ac, temp3)
      }
    }
    dist.matrix <- get_dist(data_matrix, method = ac[which.max(ac$ac),][["distance"]])
    h.clust <- hclust(dist.matrix, method = ac[which.max(ac$ac),][["method"]])
    Nb.Clust<-NbClust(data_matrix,
                      distance = ac[which.max(ac$ac),][["distance"]],
                      method = ac[which.max(ac$ac),][["method"]],
                      index = NbC.index)
    rm(list = c("i", "j", "temp1", "temp2", "temp3"))
  }else if(distance == "auto" & method != "auto"){
    d <- c("euclidean", "manhattan",
           "canberra", "binary", "minkowski")
    ac <- NULL
    for (i in d){
      temp1<-get_dist(data_matrix, method = i)
      temp2 <- hclust(temp1, method = method)
      temp3 <- tibble(distance = i, method = method, ac = coef.hclust(temp2))
      ac <- rbind(ac, temp3)
    }
    dist.matrix <- get_dist(data_matrix, method = ac[which.max(ac$ac),][["distance"]])
    h.clust <- hclust(dist.matrix, method = method)
    Nb.Clust<-NbClust(data_matrix,
                      distance = ac[which.max(ac$ac),][["distance"]],
                      method = method,
                      index = NbC.index)
    rm(list = c("i", "temp1", "temp2", "temp3"))
  }else if(distance != "auto" & method == "auto"){
    m <- c("ward.D", "ward.D2", "single",
           "complete", "average",
           "mcquitty")
    ac <- NULL
    temp1<-get_dist(data_matrix, method = distance)
    for(j in m){
      temp2 <- hclust(temp1, method = j)
      temp3 <- tibble(distance = distance, method = j, ac = coef.hclust(temp2))
      ac <- rbind(ac, temp3)
    }
    dist.matrix <- get_dist(data_matrix, method = distance)
    h.clust <- hclust(dist.matrix, method = ac[which.max(ac$ac),][["method"]])
    Nb.Clust<-NbClust(data_matrix,
                      distance = distance,
                      method = ac[which.max(ac$ac),][["method"]],
                      index = NbC.index)
    rm(list = c("j", "temp1", "temp2", "temp3"))
  }else{
    dist.matrix <- get_dist(data_matrix, method = distance)
    h.clust <- hclust(dist.matrix, method = method)
    ac <- tibble(distance = distance, method = method, ac = coef.hclust(h.clust))
    Nb.Clust<-NbClust(data_matrix,
                      distance = distance,
                      method = method,
                      index = NbC.index)
  }
  
  HClust.Analysis <- list(dist.matrix = dist.matrix,
                          h.clust = h.clust,
                          ac = ac,
                          Nb.Clust = Nb.Clust)
  return(HClust.Analysis)
}


# function to create a clusterde heatmap plot ----

fx_HClust.Heatmap <-function(data = NULL,
                             names.variables = NULL,
                             names.rows = NULL,
                             dendrogram = NULL,
                             clusters = NULL
){# Step One: Data preparation
  data_df<-as.data.frame(data[names.variables])
  row.names(data_df)<-data[[names.rows]]
  
  #dendogram data
  dendro <- dendro_data(dendrogram)
  # add cluster info to dendro object
  dendro[["labels"]][["clusters"]]<-NA
  for(i in names(clusters)){
    dendro[["labels"]][dendro[["labels"]][["label"]]==i,"clusters"] <- clusters[i]
  }
  
  # define cluster dimensions 
  clusters_dim <- NULL
  
  for(i in levels(as.factor(dendro[["labels"]][["clusters"]]))){
    temp <- tibble(cluster = i, n = nrow(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,]))
    temp$x_min <- min(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,"x"])
    temp$x_max <- max(dendro[["labels"]][dendro[["labels"]][["clusters"]]==i,"x"])
    if(i == "1"){
      temp$y_max <- as.numeric(levels(as.factor(dendro$segments$y)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
      temp$yend_max <- as.numeric(levels(as.factor(dendro$segments$yend)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
    }else{
      temp$y_max <- as.numeric(levels(as.factor(dendro$segments$y)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
      temp$yend_max <- as.numeric(levels(as.factor(dendro$segments$yend)))[nrow(dendro[["labels"]]) - as.integer(as.numeric(i)/2)-1]
    }
    clusters_dim<-rbind(clusters_dim, temp)
  }
  
  # create a function to plot the dendrogram in ggplot2 format
  ggdend <- function(df) {
    ggplot() +
      geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend, color = cluster)) +
      labs(x = "", y = "") + theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none")
  }
  
  # add clusters to segments
  dendro$segments$cluster <- NA
  for(i in clusters_dim[["cluster"]]){
    dendro[["segments"]][dendro$segments$x >= clusters_dim[clusters_dim[["cluster"]]==i,][["x_min"]] &
                           dendro$segments$xend >= clusters_dim[clusters_dim[["cluster"]]==i,][["x_min"]] &
                           dendro$segments$x <= clusters_dim[clusters_dim[["cluster"]]==i,][["x_max"]] &
                           dendro$segments$xend <= clusters_dim[clusters_dim[["cluster"]]==i,][["x_max"]] &
                           (dendro$segments$y <= clusters_dim[clusters_dim[["cluster"]]==i,][["y_max"]] |
                              dendro$segments$yend <= clusters_dim[clusters_dim[["cluster"]]==i,][["yend_max"]]),][["cluster"]]<-i
    
  }
  
  # x/y dendograms
  dendro_plot <- ggdend(dendro$segments) + coord_flip()
  # heatmap
  col.ord <- order.dendrogram(dendrogram)
  data_matrix <- as.matrix(scale(data_df[col.ord,]))
  
  data_names <- attr(data_matrix, "dimnames")
  data_df2 <- as.data.frame(data_matrix)
  colnames(data_df2) <- data_names[[2]]
  data_df2$id <- data_names[[1]]
  data_df2$id <- with(data_df2, factor(id, levels=id, ordered=TRUE))
  melt_df <- reshape2::melt(data_df2, id.vars="id")
  Heatmap_plot <- ggplot(melt_df, aes(x = variable, y = id)) +
    geom_tile(aes(fill = value))+
    theme_bw()+
    scale_fill_gradient(low="white", high="red")+
    theme(axis.text.y= element_text(size=8))+
    theme(axis.text.x= element_text(angle = 0, 
                                    vjust = 1, 
                                    size = 10.5, 
                                    hjust = .5))+
    labs(y="eaid",
         x="Metrics",
         fill = "Scale Value")+
    theme(plot.title = element_text(hjust=0.5, size=12))+
    theme(axis.title.x = element_text(size=10))+
    theme(axis.title.y = element_text(size=10))+
    theme(legend.title = element_text(size=10))+
    theme(legend.position="bottom")
  
  Clust.Heatmap.Plots <- list(dendro = dendro,
                              clusters_dim = clusters_dim,
                              dendro_plot = dendro_plot,
                              Heatmap_plot = Heatmap_plot)
  
  return(Clust.Heatmap.Plots)
}


