library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsignif)
library(viridis)
library(cowplot)
library(readxl)
library(ggnewscale)
library(ggpubr)
library(rstatix)

raw_data2 <- read_xlsx("data_impact_LST_50k_FINAL.xlsx")
df <- as_tibble(raw_data2)
df$GROUP <- recode(df$GROUP, "BRCA1_2_IDH1_2" = "BRCA1/2",
                   "BRCA1_2_IDH1_2_Other_DDR" = "DDR, Other",
                   "BRCA1_2_Only" = "BRCA1/2",
                   "BRCA1_2_Other_DDR" = "BRCA1/2",
                   "IDH1_2_Only" = "IDH1/2",
                   "IDH1_2_Other_DDR" = "DDR, Other",
                   "None" = "WT",
                   "Other_DDR_Only" = "DDR, Other")
df$GROUP <- factor(df$GROUP, levels = c("BRCA1/2", "DDR, Other", "IDH1/2", "WT"))
#df$GROUP <- recode(df$GROUP, "ATM" = "DDR, Other")
df <- df %>% mutate(is_glioma = CANCER_TYPE %in% c("Glioma", "CNS Cancer"))
df$GROUP <- df$GROUP %>% recode("WT" = "Not DDR/IDH")
df <- df %>% filter(FACETS_QC == "TRUE")

GI <- c("Colorectal Cancer", "Pancreatic Cancer",
        "Esophagogastric Cancer", "Small Bowel Cancer", "Appendiceal Cancer",
        "Tubular Adenoma of the Colon", "Ampullary Cancer")
Sarcoma <- c("Soft Tissue Cancer", "Breast Sarcoma", "Soft Tissue Sarcoma", "Uterine Sarcoma",
             "Bone Cancer")
SCC <- c("Cervical Cancer", "Head and Neck Cancer", "Anal Cancer")
Hematopoietic <- c("Lymphatic Cancer, NOS", "Myeloproliferative Neoplasms", "Leukemia",
                   "B-Lymphoblastic Leukemia/Lymphoma", "T-Lymphoblastic Leukemia/Lymphoma",
                   "Myelodysplastic Syndromes", "Non-Hodgkin Lymphoma", "Mature T and NK Neoplasms",
                   "Mature B-Cell Neoplasms")
Brain <- c("Glioma", "Miscellaneous Brain Tumor", "CNS Cancer")
Testicular <- c("Germ Cell Tumor", "Sex Cord Stromal Tumor")

df <- df %>% mutate(tumor_group = ifelse(CANCER_TYPE %in% Sarcoma, "Sarcoma", ifelse(
  CANCER_TYPE %in% SCC, "SCC", ifelse(
    CANCER_TYPE %in% Brain, "Glioma", ifelse(
      CANCER_TYPE %in% Hematopoietic, "Hematopoietic", ifelse(
        CANCER_TYPE == "Breast Cancer", "Breast", ifelse(
          CANCER_TYPE == "Non-Small Cell Lung Cancer", "NSCLC", ifelse(
            CANCER_TYPE == "Bladder Cancer", "Bladder", ifelse(
              CANCER_TYPE == "Hepatobiliary Cancer", "Hepatobiliary", ifelse(
                CANCER_TYPE == "Prostate Cancer", "Prostate", ifelse(
                  CANCER_TYPE == "Melanoma", "Melanoma", ifelse(
                    CANCER_TYPE == "Skin Cancer, Non-Melanoma", "Skin Non-Melanoma", ifelse(
                      CANCER_TYPE == "Cancer of Unknown Primary", "CUP", ifelse(
                        CANCER_TYPE == "Pancreatic Cancer", "Pancreas", ifelse(
                          CANCER_TYPE == "Colorectal Cancer", "CRC", '')))))))))))))))

wes_palette("Darjeeling2")
Darjeeling2_manual <- c("#899da4", "#046C9A", "#D69C4E", "#ECCBAE")
Darjeeling2_manual_dark <- c("#6B7A80", "#014867", "#997037", "#998472")
Rushmore_manual <- c("#899da4", "#35274a", "#0c775f", "#e1be6d")

#plot for all tumor types
plot1 <- function(dat, yvar) {
  require("dplyr")
  #yvar <- as.name(yvar)
  upper_limit = max(dat[,yvar], na.rm = TRUE) * 1
  stat_box_data <- function(x) {
    return(data.frame(y = 1.38 * upper_limit,
                      label = paste('N =', length(x))))
  }
  #compare_means(DDR_status ~ yvar, data = dat, ref.group = "DDR_wt", method = "t.test")
  ggplot(dat, mapping = aes(x = GROUP, y = .data[[yvar]], fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0.9, color = "grey50", size = 3, show.legend = FALSE) +
    stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2")), size = 3, tip.length = 0.01) +
    labs(title = "LST by Mutation - All Tumors - IMPACT Data") +
    ylab("Large Scale Transition Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(filename = paste("20200721_IMPACT_", yvar, ".pdf", sep = ""), height = 8, width = 12)
}

plot1b <- function(dat, yvar) {
  require("dplyr")
  dat <- as.data.frame(dat)
  dat$var1 <- dat[,yvar]
  #yvar <- as.name(yvar)
  upper_limit = max(dat$var1, na.rm = TRUE) * 1
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length(x))))
  }
  dunn_stat <- dat %>% dunn_test(var1 ~ GROUP, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "GROUP", step.increase = 0.05)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = TRUE, digits = 3)
  print(dunn_stat)
  kruskal_test <- compare_means(var1 ~ GROUP, data = dat, method = "kruskal")
  print(kruskal_test)
  ggplot(dat, mapping = aes(x = GROUP, y = var1, fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0.9, color = "grey50", size = 3, show.legend = FALSE) +
    stat_pvalue_manual(dunn_stat, label = "p.adj.sci", hide.ns = TRUE, inherit.aes = FALSE,
                       size = 3, tip.length = 0.01) +
    #stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2")), size = 3, tip.length = 0.01) +
    labs(title = "LST by Mutation - All Tumors - IMPACT Data") +
    ylab("Large Scale Transition Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(filename = paste("20200721_IMPACT_", yvar, ".pdf", sep = ""), height = 8, width = 12)
}
plot1b(df, "LST")

#plot for select tumor types
plot2 <- function(tumor_group, na.rm = TRUE, ...) {
  data_f <- df %>% filter(.$tumor_group == {{tumor_group}})
  df_filtered_WT <- data_f %>% filter(.$GROUP == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[["LST"]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  data_f %>% 
    ggplot(mapping = aes(x = GROUP, y = .data[["LST"]], fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    # scale_color_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    # scale_fill_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5, scale = "width") +
    # geom_violin(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", alpha = 1/5, scale = "width") +
    # geom_boxplot(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    #facet_wrap(~condition_combined, nrow = 1) +
    #stat_compare_means(comparisons = list(c("0.5% O2", "21% O2")), tip.length = 0.02) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2")), size = 3, tip.length = 0.01) +
    ylim(0,47) +
    scale_x_discrete(labels = c("BRCA1/2" = "BRCA", "DDR, Other" = "DDR", "IDH1/2" = "IDH", "Not DDR/IDH" = "unalt")) +
    labs(title = paste0("LST - ", {{tumor_group}})) +
    ylab("LST Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(file = paste0("20220104 IMPACT_", tumor_short, ".png"), width = 7, height = 5)
}

plot2b <- function(tumor_group, na.rm = TRUE, ...) {
  data_f <- df %>% as.data.frame %>% filter(.$tumor_group == {{tumor_group}})
  df_filtered_WT <- data_f %>% filter(.$GROUP == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[["LST"]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.4 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  dunn_stat <- data_f %>% dunn_test(LST ~ GROUP, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "GROUP", step.increase = 0.05)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = 2)
  print(dunn_stat)
  kruskal_test <- compare_means(LST ~ GROUP, data = data_f, method = "kruskal")
  print(kruskal_test)
  data_f %>% 
    ggplot(mapping = aes(x = GROUP, y = .data[["LST"]], fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    # scale_color_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    # scale_fill_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5, scale = "width") +
    # geom_violin(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", alpha = 1/5, scale = "width") +
    # geom_boxplot(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    #facet_wrap(~condition_combined, nrow = 1) +
    #stat_compare_means(comparisons = list(c("0.5% O2", "21% O2")), tip.length = 0.02) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    # stat_compare_means(method = "kruskal.test") +
    stat_pvalue_manual(dunn_stat, label = "{scales::pvalue(p.adj, accuracy = 0.0001)}", hide.ns = FALSE, inherit.aes = FALSE,
                       size = 2.5, tip.length = 0.01) +
    scale_x_discrete(labels = c("BRCA1/2" = "BRCA", "DDR, Other" = "DDR", "IDH1/2" = "IDH", "Not DDR/IDH" = "unalt")) +
    labs(title = paste0("LST - ", {{tumor_group}})) +
    ylab("LST Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(file = paste0("20220104 IMPACT_", tumor_short, ".png"), width = 7, height = 5)
}

plot2b("Breast")

plot2c <- function(tumor_group, na.rm = TRUE, ...) {
  data_f <- df %>% as.data.frame %>% filter(.$tumor_group == {{tumor_group}})
  df_filtered_WT <- data_f %>% filter(.$GROUP == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[["LST"]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  dunn_stat <- data_f %>% dunn_test(LST ~ GROUP, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "GROUP", step.increase = 0.05)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = 2)
  print(dunn_stat)
  kruskal_test <- compare_means(LST ~ GROUP, data = data_f, method = "kruskal")
  print(kruskal_test)
  data_f %>% 
    ggplot(mapping = aes(x = GROUP, y = .data[["LST"]], fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    # scale_color_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    # scale_fill_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5, scale = "width") +
    # geom_violin(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", alpha = 1/5, scale = "width") +
    # geom_boxplot(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    #facet_wrap(~condition_combined, nrow = 1) +
    #stat_compare_means(comparisons = list(c("0.5% O2", "21% O2")), tip.length = 0.02) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    # stat_compare_means(method = "kruskal.test") +
    stat_pvalue_manual(dunn_stat, label = "{scales::pvalue(p.adj, accuracy = 0.0001)}", hide.ns = FALSE, inherit.aes = FALSE,
                       size = 2.5, tip.length = 0.01) +
    scale_x_discrete(labels = c("BRCA1/2" = "BRCA", "DDR, Other" = "DDR", "IDH1/2" = "IDH", "Not DDR/IDH" = "unalt")) +
    labs(title = paste0("LST - ", {{tumor_group}})) +
    ylab("LST Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(file = paste0("20220104 IMPACT_", tumor_short, ".png"), width = 7, height = 5)
}

custom_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        #         text = element_text(size=16), 
        axis.title = element_text(face = "bold"),
        #         axis.text.x = element_text(angle = 90, hjust = 1),
        #         plot.caption = element_text(face = "bold", size = "20"), 
        legend.position='none')

#plot for not tumor types e.g. not glioma
plot3 <- function(tumor_group, na.rm = TRUE, ...) {
  data_f <- df %>% filter(.$tumor_group != {{tumor_group}})
  df_filtered_WT <- data_f %>% filter(.$GROUP == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[["LST"]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  data_f %>% 
    ggplot(mapping = aes(x = GROUP, y = .data[["LST"]], color = GROUP, fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    # scale_color_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    # scale_fill_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5) +
    # geom_violin(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    #facet_wrap(~condition_combined, nrow = 1) +
    #stat_compare_means(comparisons = list(c("0.5% O2", "21% O2")), tip.length = 0.02) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2")), size = 3, tip.length = 0.01) +
    ylim(0,47) +
    scale_x_discrete(labels = c("BRCA1/2" = "BRCA", "DDR, Other" = "DDR", "IDH1/2" = "IDH", "Not DDR/IDH" = "unalt")) +
    labs(title = paste0("LST - Not ", {{tumor_group}})) +
    ylab("LST Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(file = paste0("20220104 IMPACT_", tumor_short, ".png"), width = 7, height = 5)
}

plot3b <- function(tumor_group, na.rm = TRUE, ...) {
  data_f <- df %>% as.data.frame %>% filter(.$tumor_group != {{tumor_group}})
  df_filtered_WT <- data_f %>% filter(.$GROUP == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[["LST"]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.4 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  dunn_stat <- data_f %>% dunn_test(LST ~ GROUP, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "GROUP", step.increase = 0.05)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = TRUE, digits = 3)
  print(dunn_stat)
  kruskal_test <- compare_means(LST ~ GROUP, data = data_f, method = "kruskal")
  print(kruskal_test)
  data_f %>% 
    ggplot(mapping = aes(x = GROUP, y = .data[["LST"]], fill = GROUP)) +
    geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
    # scale_color_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    # scale_fill_manual(values = wes_palette(5, name = "Rushmore1", type = "continuous")) +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = GROUP), alpha = 1/5, scale = "width") +
    # geom_violin(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", alpha = 1/5, scale = "width") +
    # geom_boxplot(data = df_filtered_WT, aes(x = GROUP, y = .data[["LST"]]), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    #facet_wrap(~condition_combined, nrow = 1) +
    #stat_compare_means(comparisons = list(c("0.5% O2", "21% O2")), tip.length = 0.02) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    # stat_compare_means(method = "kruskal.test") +
    stat_pvalue_manual(dunn_stat, label = "p.adj.sci", hide.ns = FALSE, inherit.aes = FALSE,
                       size = 2.5, tip.length = 0.01) +
    scale_x_discrete(labels = c("BRCA1/2" = "BRCA", "DDR, Other" = "DDR", "IDH1/2" = "IDH", "Not DDR/IDH" = "unalt")) +
    labs(title = paste0("LST - Not ", {{tumor_group}})) +
    ylab("LST Score") +
    xlab(element_blank()) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(lineend = "square"),
          #         text = element_text(size=16), 
          axis.title = element_text(face = "bold"),
          #         axis.text.x = element_text(angle = 90, hjust = 1),
          #         plot.caption = element_text(face = "bold", size = "20"), 
          legend.position='none')
  #         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
  #         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
  #ggsave(file = paste0("20220104 IMPACT_", tumor_short, ".png"), width = 7, height = 5)
}

plot2b("CRC")

#prepare data for outlines of Not DDR/IDH group
df_WT <- filter(df, GROUP == "Not DDR/IDH")

#figure1_panelA all tumor types
all_tumors <- plot1(df, "LST")
all_tumors

#All tumors with Kruskal-Wallis test stats
upper_limit = max(df$LST, na.rm = TRUE) * 1
stat_box_data <- function(x) {
  return(data.frame(y = 1.3 * upper_limit,
                    label = paste('N =', length(x))))
}
compare_means(LST ~ GROUP, df, method = "kruskal")
dunn_stats <- df %>% dunn_test(LST ~ GROUP, p.adjust.method = "hochberg")
dunn_stats <- dunn_stats %>% add_xy_position(x = "GROUP", step.increase = 0.05)
dunn_stats$p.adj.sci <- format(dunn_stats$p.adj, scientific = 2)
dunn_stats <- dunn_stats %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
print(dunn_stats)
all_tumors2 <- ggplot(df, mapping = aes(x = GROUP, y = LST, fill = GROUP)) +
  geom_point(aes(color = GROUP), size = 2.5, alpha = 1/10, position = "jitter") +
  scale_color_manual(values = rev(Darjeeling2_manual)) +
  new_scale_color() +
  geom_boxplot(aes(color = GROUP), outlier.shape = NA, alpha = 1/5) +
  geom_violin(aes(color = GROUP), alpha = 1/5) +
  # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
  # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
  scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
  scale_fill_manual(values = rev(Darjeeling2_manual)) +
  stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0.9, color = "grey50", size = 3, show.legend = FALSE) +
  #stat_compare_means(method = "kruskal.test") +
  stat_pvalue_manual(dunn_stats, hide.ns = TRUE, label = "p.adj.sci", inherit.aes = FALSE, tip.length = 0.01, size = 3, step_increase = 0.05) +
  #stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2")), size = 3, tip.length = 0.01) +
  labs(title = "LST by Mutation - All Tumors - IMPACT Data") +
  ylab("Large Scale Transition Score") +
  xlab(element_blank()) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        #         text = element_text(size=16), 
        axis.title = element_text(face = "bold"),
        #         axis.text.x = element_text(angle = 90, hjust = 1),
        #         plot.caption = element_text(face = "bold", size = "20"), 
        legend.position='none')
#         strip.background = element_rect(color = "#899da4", fill="#899da4", size=1, linetype="solid"),
#         strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"))
#ggsave(filename = paste("20200721_IMPACT_", yvar, ".pdf", sep = ""), height = 8, width = 12)

all_tumors3 <- plot1b(df, "LST")
plot3b("Glioma")

#figure1_panelsB:E tumor subtypes

tumor_group <- unique(df$tumor_group)
tumor_types_IMPACT <- tibble(tumor_group)
list_by_tumor_type <- tumor_types_IMPACT %>% pmap(plot2b)
names(list_by_tumor_type) <- tumor_types_IMPACT %>% as_vector()

list_by_tumor_type2 <- tumor_types_IMPACT %>% pmap(plot2c)
names(list_by_tumor_type2) <- tumor_types_IMPACT %>% as_vector()

##not glioma
not_glioma_list <- tibble(tumor_group = list("Glioma"))
not_glioma <- not_glioma_list %>% pmap(plot3b)






IMPACT_LST_all_tumors <- all_tumors
IMPACT_LST_all_tumors2 <- all_tumors2
IMPACT_LST_all_tumors3 <- all_tumors3
IMPACT_LST_list_by_tumor_type <- list_by_tumor_type
IMPACT_LST_list_by_tumor_type2 <- list_by_tumor_type2
IMPACT_LST_not_glioma <- not_glioma

### save files for figure
save(IMPACT_LST_all_tumors,
     IMPACT_LST_all_tumors2,
     IMPACT_LST_all_tumors3,
     IMPACT_LST_list_by_tumor_type,
     IMPACT_LST_list_by_tumor_type2,
     IMPACT_LST_not_glioma,
     file = "20220408 IMPACT LST all plots.rdata")