# Load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(cocor)
library(pROC)
library(ggpubr)
library(cowplot)

# Set colors
BAR_COL <- c("#FB8072", "#80B1D3", "#FDB462", "#CCCCFF", 
             "#F19C79", "#006D77", "#6A3569")
GROUP_COL <- c("#A5170E", "#1965B0", "#F4A736", "#4EB265", 
              "#E6BEFF", "#9A6324", "#FFFAC8", "#800000")
ROC_COL <- c("#E41A1C","#377EB8")

# Function 1: test for R2 or AUC
perform.test <- function(pheno = NULL, 
                         pred = NULL,
                         type = "gaussian"
){
  
  if (type == "gaussian") {
    
    perform_test <- cor.test(pheno, pred)
    perform <- c(perform_test$estimate, perform_test$conf.int)^2
  } else {
    
    perform_test <- roc(pheno, pred, levels = c(0, 1))
    perform <- ci.auc(perform_test)[c(2, 1, 3)]
  }
  return(list("test" = perform_test,
              "perform" = perform))
}

# Function 2: density plot
dens.plt <- function(pheno = NULL,
                     pred = NULL,
                     trait_label = NULL,
                     type = "gaussian", 
                     watermark = FALSE
){
  
  # R square
  r2 <- perform.test(pheno = pheno, 
                     pred = pred, 
                     type = type)[["perform"]][1] %>% round(4)
  datt <- data.frame(Value = c(pheno, pred),
                     Group = factor(c(rep(trait_label, length(pheno)),
                                      rep("PGS", length(pred))), 
                                    levels = c(trait_label, "PGS")))
  anno_lab <- bquote(italic(R)^2~"="~.(r2))
  # plot
  plt <- ggplot(datt) +
    geom_density(aes(x = Value, fill = Group, color = Group),
                 alpha = 0.2, size = 0.6) +
    scale_fill_manual(values = ROC_COL) + 
    scale_color_manual(values = ROC_COL) +
    ylab("Density") + xlab("") +
    theme_bw() +
    annotate("text", colour="grey10", label = anno_lab,
             x = 0, y = 0.5, size = 5) + 
    theme(legend.position = "inside",
          legend.position.inside=c(.8, .9), 
          panel.grid = element_blank(),
          legend.title = element_text(size = 12, face = "bold", color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          title = element_text(size = 15, face = "bold", color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"))
  if (watermark)
    plt <- plt + 
    draw_label("Created by PGSFusion\nUKBB application 144904",size = 30, angle = 50, 
               color ="grey", alpha = 0.7)
  return(plt)
}

# Function 3: ROC plot
roc.plt <- function(pheno = NULL,
                    pred = NULL,
                    trait_label = NULL, 
                    type = "b", 
                    watermark = FALSE
){
  
  # AUC
  perform_test <- perform.test(pheno = pheno, 
                               pred = pred, 
                               type = type)
  auc_vec <- perform_test[["perform"]] %>% round(4)
  auc <- paste0("AUC of ",
                trait_label, ": \n",
                auc_vec[1] , "(",
                auc_vec[2], "~", 
                auc_vec[3], ")")
  roc <- perform_test[["test"]]
  datt <- data.frame(Model = "ROC",
                     rev_specificity = 1 - roc$specificities,
                     sensitivity = roc$sensitivities)
  # plot
  plt <- ggplot(datt) +
    geom_line(aes(x = rev_specificity, y = sensitivity), 
              color = "#4DAF4A", linetype = 2, size = 1.3) +
    geom_abline(slope = 1, intercept = 0, 
                color = "grey10", linetype = 2)  + 
    xlab("1-Specificity") + ylab("Sensitivity") +
    theme_bw() + 
    annotate("text", colour="grey10", label = auc,
               x = 0.7, y = 0.2, size = 5) +
    theme(axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"))
  if (watermark)
    plt <- plt  + 
    draw_label("Created by PGSFusion\nUKBB application 144904",size = 30, angle = 50, 
               color ="grey", alpha = 0.7)
  return(plt)
}

# Function 4: PGS trend plot 
PGS.trend.plt <- function(datt = NULL,
                          n = NULL,
                          type = "gaussian", 
                          watermark = FALSE
){
  
  ## fit model
  datt <- subset(datt, rowSums(is.na(datt)) == 0)
  datt$PGS_g <- as.factor(datt$PGS_g)
  fomu_cov_PGSg <- paste0("y ~ .")
  glm_cov_PGSg <- glm(fomu_cov_PGSg, data = datt, family = type)
  
  ## format estimation
  summ_coef_PGSg <- summary(glm_cov_PGSg)$coefficients[paste0("PGS_g", 2:n), ]
  est_PGS_g <- data.frame(est = c(0, summ_coef_PGSg[,1]),
                          lowCI = c(0, (summ_coef_PGSg[,1] - 1.96*summ_coef_PGSg[,2])),
                          highCI = c(0, (summ_coef_PGSg[,1] + 1.96*summ_coef_PGSg[,2])))
  if (type == "binomial") 
    est_PGS_g <- exp(est_PGS_g) %>% as.data.frame()
  high <- seq(0, 100, length = n + 1)[-1] %>% round()
  low <- c(0, high[-n] + 1)
  est_PGS_g$Group <- seq(0, 100, length = n + 1)[-1] %>% round()

  x_lab <- "PGS (%)"
  y_lab <- ifelse(type == "binomial", "OR (95% CI)", "Beta (95% CI)")
  ref_y <- ifelse(type == "binomial", 1, 0)
  title_lab <- paste0("PGS by ", n, " Groups")
  
  ## p trend for OR/Beta
  p_trend <- mk.test(est_PGS_g$est, continuity = TRUE)$p.value
  p_trend <- format(p_trend, scientific = T, digits = 4) %>%
    gsub("e", "E", .)
  anno_lab <- bquote(italic(P)[trend]~"="~.(p_trend))
  
  ## plot
  plt <- ggplot(est_PGS_g, aes(x = Group, y = est, color = Group)) +
    geom_hline(yintercept = ref_y, 
               color = "grey50", 
               linetype = 2, 
               size = 1) +
    geom_point(show.legend = F) + 
    geom_pointrange(aes(ymin = lowCI, ymax = highCI),
                    stat='identity',
                    position = position_dodge(0.4),
                    show.legend = F,
                    alpha = 1, 
                    size = 1) +
    annotate("text",
             x = 85, y = min(est_PGS_g$lowCI) - 0.1,
             size = 5, colour="grey10",
             label = anno_lab) +
    scale_x_continuous(breaks = est_PGS_g$Group, 
                       labels = paste0(low,  "~", high)) + 
    scale_color_viridis_c() +
    stat_smooth(method = "lm",
                se = F, linetype = 2, color = "red") +
    xlab(x_lab) + ylab(y_lab) +
    ggtitle(title_lab) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"),
          title = element_text(size = 15, face = "bold", color = "black"))
  if (watermark)
    
    plt <- plt + 
      annotate(geom = "text",
               x = median(est_PGS_g$Group), 
               y = mean(c(min(est_PGS_g$lowCI), max(est_PGS_g$highCI))), 
               label = 'Created by PGSFusion\nUKBB application 144904', 
               color = 'grey', angle = 45, 
               size = 10, alpha = 0.5)
   
  return(plt)
}

# Function 5: PGS in strata plot
PGS.strata.plt <- function(datt = NULL,
                           strata_var = NULL, 
                           watermark = FALSE
){
  
  # plot parameters
  x_lab <- ""
  y_lab <- "PGS"
  title_lab <- paste0("PGS by ", strata_var)
  n_group <- length(unique(datt$Group))
  # p value
  datt <- subset(datt, rowSums(is.na(datt)) == 0)
  if (n_group == 1) {
    anno_lab <- NULL
  } else {
    
    if (n_group == 2) {
      p_diff <- wilcox.test(PGS~Group, datt)$p.value
    } else {
      p_diff <- kruskal.test(PGS~Group, datt)$p.value
    }
    p_diff <- p_diff %>% 
      format(., scientific = T, digits = 4) %>%
      gsub("e", "E", .)
    anno_lab <- bquote(italic(P)~"="~.(p_diff))
  }
  # plot
  plt <- ggplot(datt) + 
    geom_boxplot(aes(x = Group, y = PGS, color = Group),
                 width = 0.3, size = 1.2, fill = NA) + 
    scale_color_manual(values = BAR_COL) + 
    annotate("text", 
             x = n_group * 0.5, y = max(datt$PGS)*1.1, 
             size = 5, colour="grey10",
             label = anno_lab) + 
    xlab(x_lab) + ylab(y_lab) +
    ggtitle(title_lab) +
    theme_bw() + 
    theme(legend.position = "none",
          title = element_text(size = 15, face = "bold", color = "black"),
          axis.text.x = element_text(size = 12, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"))
  if(watermark)
    plt <- plt  +
      annotate(geom = "text",
               x = mean(as.numeric(datt$Group)), y = mean(datt$PGS),
               label = 'Created by PGSFusion\nUKBB application 144904',
               color = 'grey', angle = 45,
               size = 8, alpha = 0.5)
  
  return(plt)
}

# Function 6: compare prediction performance
perform.compare <- function(model1 = NULL,
                            model2 = NULL,
                            type = "gaussian"
){
  
  if (type == "gaussian") {
    perform_compare <- cocor.indep.groups(r1.jk = model1$estimate,
                                          r2.hm = model2$estimate,
                                          n1 = model1$parameter + 2,
                                          n2 = model2$parameter + 2,
                                          alternative = "two.sided",
                                          test = "fisher1925",
                                          alpha = 0.05,
                                          conf.level = 0.95,
                                          return.htest = F)
    return(perform_compare@fisher1925$p.value)
  } else {
    
    perform_compare <- roc.test(roc1 = model1, 
                                roc2 = model2, 
                                method = "bootstrap", 
                                alternative = "two.sided",
                                boot.n = 1000, 
                                progress = "none",
                                conf.level = 0.95)
    return(perform_compare$p.value)
  }
}

# Function 7: PGS in strata
perform.strata.test <- function(pheno = NULL,
                                pred = NULL,
                                type = "gaussian",
                                subgroup = NULL
){
  
  ## format input
  datt <- data.frame(pheno = pheno,
                     pred = pred,
                     Group = subgroup)
  datt <- subset(datt, rowSums(is.na(datt)) == 0)
  datt_group_list <- split(datt, 
                           f = ~ datt$Group)
  
  ## calculate performance for each group
  perform_test_strata <- lapply(datt_group_list, function(datt_group_listx){
    
    performx <- tryCatch({
      performx <- perform.test(pheno = datt_group_listx$pheno, 
                               pred = datt_group_listx$pred,
                               type = type)
    }, error = function(e){
      performx <- list("test" = NA,
                       "perform" = rep(NA, 3))
    })
    
    return(performx)
  })
  ## format prediction performance
  perform_strata <- lapply(perform_test_strata, function(x){
    x[["perform"]]
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  colnames(perform_strata) <- c("perform", "lowCI", "highCI")
  perform_strata$Group <- factor(names(datt_group_list), levels = names(datt_group_list))
  ## compare performance in strata
  all_comp_pair <- combn(names(perform_test_strata), 2)
  n_pair <- ncol(all_comp_pair)
  perform_strata_test <- lapply(1:ncol(all_comp_pair), function(comb_x){
    
    pairx <- all_comp_pair[,comb_x, drop = T]
    perform_test_x1 <- perform_test_strata[[pairx[1]]]
    perform_test_x2 <- perform_test_strata[[pairx[2]]]
    
    if (any(is.na(perform_test_x1[["perform"]] + 
                  perform_test_x2[["perform"]]))) {
      
      p_performx <- NA
    } else {
      
      p_performx <- perform.compare(model1 = perform_test_x1[["test"]],
                                    model2 = perform_test_x2[["test"]],
                                    type = type)
    }
    return(matrix(c(pairx, p_performx), nrow = 1))
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  ## format test data frame
  colnames(perform_strata_test) <- c("group1", "group2", "p")
  perform_strata_test$p.adj <- p.adjust(perform_strata_test$p, method = "BH") %>% 
    round(3)
  perform_strata_test$p.adj[perform_strata_test$p.adj < 0.001] <- "< 0.001"
  ref_perform <- ifelse(type == "gaussian", 0, 0.5)
  range_perform <- max(perform_strata$highCI) - min(c(ref_perform, perform_strata$lowCI))
  perform_strata_test$y.position <- max(perform_strata$highCI) + 
    (c(n_pair:1) - 1) * range_perform * 0.08
  
  return(list("perform_strata" = perform_strata,
              "test" = perform_strata_test))
}

# Function 8: performance strata plot
perform.strata.plt <- function(datt = NULL,
                               test = NULL,
                               type = "gaussian",
                               strata_var = NULL, 
                               watermark = FALSE
){
  
  # plot parameters
  x_lab <- ""
  y_lab <- ifelse(type == "gaussian", 
                  as.expression(bquote(italic(R)^2~" (95% CI)")),
                  "AUC (95% CI)")
  ref_y <- ifelse(type == "gaussian", 0, 0.5)
  title_lab <- paste0("Performance by ", strata_var)
  
  # plot 
  datt <- subset(datt, rowSums(is.na(datt)) == 0)
  plt <- ggplot(datt, aes(x = Group, y = perform, color = Group)) + 
    geom_hline(yintercept = ref_y, 
               color = "grey50", 
               linetype = 2, 
               size = 1) +
    geom_pointrange(aes(ymin = lowCI, ymax = highCI),
                    stat='identity',
                    alpha = 1, 
                    show.legend = F, 
                    position = position_dodge(0.6),
                    size = 1) +
    scale_x_discrete(breaks = datt$Group) +
    xlab(x_lab) + ylab(y_lab) +
    ggtitle(title_lab) + 
    scale_color_manual(values = GROUP_COL) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"),
          title = element_text(size = 15, face = "bold", color = "black"))
  ## add test
  if (!is.null(test)) 
    plt <- plt + stat_pvalue_manual(test, label = "p.adj", tip.length = 0.01, 
                                    label.size = 4, bracket.size = 0.6)
  if(watermark)
    plt <- plt  + annotate(geom="text",
                         x = mean(as.numeric(datt$Group)), 
                         y = mean(c(ref_y, max(datt$highCI))),
                         label='Created by PGSFusion\nUKBB application 144904', 
                         color='grey', angle = 45, 
                         size = 8, alpha=0.5)
  return(plt)
}

# Function 9: effect strata data
effect.strata.plt <- function(datt = NULL,
                              strata_var = NULL,
			                        watermark = FALSE,
                              model = "gaussian"
){
  
  ## model parameters
  datt <- subset(datt, rowSums(is.na(datt)) == 0)
  colnames(datt)[which(str_detect(colnames(datt),'NA') | is.na(colnames(datt)))]<-paste0('noname',which(str_detect(colnames(datt),'NA') | is.na(colnames(datt))))
  if (strata_var == "Sex")
    datt$V2 <- NULL
  if (strata_var == "Age_65")
    datt$V1 <- NULL
  use_cov <- setdiff(colnames(datt), c("pheno", "Group", "PGS"))
  n_group <- length(unique(datt$Group))
  ## p value of interaction 
  datt_inter <- datt
  datt_inter$Group <- as.numeric(datt_inter$Group)
  fomu_cov_PGS_inter <- paste0("pheno ~ ", 
                               paste(c(use_cov, "PGS*Group"), collapse = " + "))
  glm_cov_PGS_inter <- glm(as.formula(fomu_cov_PGS_inter), data = datt_inter, family = model)
  p_all <- summary(glm_cov_PGS_inter)$coefficients[, 4]
  p_plt <- p_all[grep("PGS:Group", names(p_all))] %>% 
    format(., scientific = T, digits = 4) %>%
    gsub("e", "E", .)
  summ_out <- summary(glm_cov_PGS_inter)$coefficients[, c(1, 2, 4)] %>% as.data.frame()
  summ_out <- summ_out[grep("PGS|Group", row.names(summ_out)), ]  
  colnames(summ_out) <- c("effect", "sd", "P")
  summ_ci <- confint(glm_cov_PGS_inter, parm = c("PGS", "Group", "PGS:Group"), 
                     type = "Wald")
  summ_out$CI_low <- summ_ci[, 1] %>% round(.,4)
  summ_out$CI_high <- summ_ci[, 2] %>% round(.,4)
  summ_out_i <- data.frame(effect = summ_out$effect %>% round(.,4), 
                           sd = summ_out$sd %>% round(.,4), 
                           CI_low = summ_out$CI_low, 
                           CI_high = summ_out$CI_high, 
                           P = summ_out$P %>% format(., scientific = T, digits = 4)) 
  summ_out_i <- cbind.data.frame(Variable = c("PGS", strata_var, "Interaction"), summ_out_i)
  anno_lab <- bquote(italic(P)[interaction]~"="~.(p_plt))
  
  ## fit model in each strata
  fomu_cov_PGS <- paste0("pheno ~ ", 
                         paste(use_cov, collapse = " + "), " + ",
                         "PGS")
  est_PGS <- lapply(levels(datt$Group), function(gx){

    dattx <- subset(datt, datt$Group == gx)
    glm_cov_PGSx <- glm(fomu_cov_PGS, data = dattx, family = model)
    est_PGSx <- c(summary(glm_cov_PGSx)$coefficients["PGS", 1], 
                  confint(glm_cov_PGSx, parm = "PGS", type = "Wald"))

    return(est_PGSx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  colnames(est_PGS) <- c("est", "lowCI", "highCI")
  if (model == "binomial") 
    est_PGS <- exp(est_PGS) %>% as.data.frame()
  est_PGS$Group <- levels(datt$Group) %>% factor(., levels = levels(datt$Group))

  ## other parameters for plot
  x_lab <- ""
  y_lab <- ifelse(model == "binomial", "OR (95% CI)", "Beta (95% CI)")
  title_lab <- paste0("Strata by ", strata_var)
  ## plot
  plt <- ggplot(est_PGS, aes(x = Group, y = est))+
    geom_bar(aes(color = Group),
             stat = "identity", width = 0.5, size = 1,
             fill = NA) +
    scale_color_manual(values = BAR_COL) + 
    geom_errorbar(aes(ymin = lowCI, ymax = highCI),
                  size = 0.8, width = 0.2, color = "grey10") +
    annotate("text", 
             x = n_group * 0.6, y = max(est_PGS$highCI)*1.1, 
             size = 5, colour="grey10",
             label = anno_lab) + 
    xlab(x_lab) + ylab(y_lab) +
    ggtitle(title_lab) +
    theme_bw() + 
    theme(legend.position = "none",
          title = element_text(size = 15, face = "bold", color = "black"),
          axis.text.x = element_text(size = 12, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", color = "black"))
  if(watermark)
    plt <- plt  + annotate(geom="text",
                            x = mean(as.numeric(est_PGS$Group)), 
                            y = mean(c(0, max(est_PGS$highCI))),
                            label='Created by PGSFusion\nUKBB application 144904', 
                            color='grey', angle = 45, 
                            size = 6, alpha=0.5)
  
  return(list(plt, summ_out_i))
}

# Function 10: subgroup plot
subgroup.plt <- function(datt = NULL,
                         type = NULL,
                         sub_var = NULL,
                         eff_var = NULL, 
                         watermark = FALSE
){
  
  model <- ifelse(type == "gaussian", "gaussian", "binomial")
  # sub eff_var by sub_var
  var_use_sub <- setdiff(colnames(datt), c("y", sub_var, "Group"))
  fomu_joint <- paste0("y ~ ",
                       paste(var_use_sub, collapse = " + "))
  est_sub <- lapply(levels(datt[[sub_var]]), function(subx){
    
    datt_sub <- subset(datt, datt[[sub_var]] == subx)
    glm_subx <- glm(fomu_joint, data = datt_sub, family = model)
    summ_coef_allx <- summary(glm_subx)$coefficients %>% as.data.frame()
    summ_coef_subx <- summ_coef_allx[grep(eff_var, rownames(summ_coef_allx)),]
    rownames(summ_coef_subx) <- gsub(paste0("^", eff_var), "", 
                                     rownames(summ_coef_subx))
    est_subx<- matrix(NA, ncol = 4, nrow = length(levels(datt_sub[[eff_var]]))) %>% 
      as.data.frame()
    est_subx[1,] <- c(rep(0, 3), 1)
    dimnames(est_subx) <- list(levels(datt_sub[[eff_var]]),
                               c("est", "lowCI", "highCI", "P"))
    est_subx[rownames(summ_coef_subx),] <- 
      data.frame(est = summ_coef_subx[,1],
                 lowCI = summ_coef_subx[,1] - 1.96*summ_coef_subx[,2],
                 highCI = summ_coef_subx[,1] + 1.96*summ_coef_subx[,2],
                 P = summ_coef_subx[,4])
    #
    if (model == "binomial") {
      est_subx[, c("est", "lowCI", "highCI")] <- 
        exp(est_subx[, c("est", "lowCI", "highCI")]) %>% as.data.frame()
    }
    # add group
    est_subx <- tibble::rownames_to_column(est_subx, "eff_var")
    est_subx$sub_var <- subx
    
    return(est_subx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  est_sub$eff_var <- factor(est_sub$eff_var, levels = levels(datt[[eff_var]]))
  est_sub$sub_var <- factor(est_sub$sub_var, levels = levels(datt[[sub_var]]))
  ##
  x_lab <- ""
  y_lab <- ifelse(model == "binomial", "OR (95% CI)", "Beta (95% CI)")
  ref_y <- ifelse(model == "binomial", 1, 0)
  eff_var_title <- ifelse(eff_var == "PGS_g", "PGS", eff_var)
  eff_var_legend <- ifelse(eff_var == "PGS_g", "PGS groups", eff_var)
  sub_var_title <- ifelse(sub_var == "PGS_g", "PGS", sub_var)
  title_lab <- paste0("Effects of ", eff_var_title, " in different ", sub_var_title, " groups")
  ## plot
  plt <- ggplot(est_sub, aes(x = sub_var, y = est, color = eff_var)) +
    geom_hline(yintercept = ref_y, 
               color = "grey50", 
               linetype = 2, 
               size = 1) +
    geom_pointrange(aes(ymin = lowCI, ymax = highCI),
                    stat='identity',
                    position = position_dodge(0.4),
                    show.legend = T,
                    alpha = 1, 
                    size = 1) +
    guides(color = guide_legend(title = eff_var_legend))+
    xlab(x_lab) + ylab(y_lab) +
    ggtitle(title_lab) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"),
          title = element_text(size = 15, face = "bold", color = "black"))
  if(watermark)
    plt <- plt  + annotate(geom="text",
                           x = mean(as.numeric(est_sub$sub_var)), 
                           y = mean(c(min(est_sub$lowCI), max(est_sub$highCI))),
                           label='Created by PGSFusion\nUKBB application 144904', 
                           color='grey', angle = 45, 
                           size = 6, alpha=0.5)
  ## add color
  if (eff_var == "PGS_g") {
    
    plt <- plt + scale_color_viridis_d()
  } else {
    
    plt <- plt + scale_color_manual(values = GROUP_COL)
  }
  est_sub_eff <- est_sub[which(est_sub$eff_var != levels(datt[[eff_var]])[1]),]
  est_sub_eff$eff_var_comb <- paste0(est_sub_eff$eff_var, 
                                     " VS ", 
                                     levels(datt[[eff_var]])[1])
  sub_tab <- data.frame(Effect_variable = est_sub_eff$eff_var_comb,
                        Strata_variable = est_sub_eff$sub_var,
                        Estimation = paste0(round(est_sub_eff$est, 3),
                                            "(",
                                            round(est_sub_eff$lowCI, 3),
                                            "~",
                                            round(est_sub_eff$highCI, 3),
                                            ")"),
                        P = formatC(est_sub_eff$P, format = "E", digits = 3))
  return(list("plt" = plt,
              "tab" = sub_tab))
}
