library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsignif)
library(viridis)
library(cowplot)
library(readxl)
library(ggnewscale)
library(ggbreak)

scores <- read.table("TCGA_DDR_Data_Resources/Scores.tsv")
DDRscores <- read.table("TCGA_DDR_Data_Resources/DDRscores.tsv")
samples <- read.table("TCGA_DDR_Data_Resources/Samples.tsv")
genes <- read.table("TCGA_DDR_Data_Resources/Genes.tsv")
gene_mutations <- read.table("TCGA_DDR_Data_Resources/GeneMutations.tsv")
gene_deletions <- read.table("TCGA_DDR_Data_Resources/GeneDeletions.tsv")
gene_alterations <- read.table("TCGA_DDR_Data_Resources/GeneAlterations.tsv")
TCGA_LST <- read_xlsx("/Users/juanma/GDrive/LabData/cbioportal/TCGA_HRD_Scores.xlsx", skip = 1, col_names = c("sample", "NtAI", "LST", "HRD_LOH", "Nmut", "FLOH", "wGII", "Ploidy")) %>% 
  mutate(sample = paste(sample, "-01", sep = ""))
TCGA_IDH2 <- read_xlsx("/Users/juanma/GDrive/LabData/R Projects/IDH_TMB/new_TCGA_from_Henry/TCGA_HRD_Scores_with_IDH2.xlsx",
                       skip = 1,
                       col_names = c("patient", "sample", "tumor_type", "subtype", "HRD_LST", "IDH2")) %>% select(sample, IDH2)
TCGA_cosmic <- read_xlsx("/Users/juanma/GDrive/LabData/R Projects/IDH_TMB/new_TCGA_from_Henry/TCGA_WES_sigProfiler_SBS_signatures_in_samples.xlsx") %>% 
  select(-c(tumor_type, Accuracy))

scores_v <- as.vector(scores$V1)
colnames(DDRscores) <- scores_v
DDRscores$sample <- samples$V1
DDRscores$tumor_type <- samples$V2

gene_mutations_t <- t(gene_mutations)
colnames(gene_mutations_t) <- genes$V1
selected_genes <- gene_mutations_t[,c("IDH1", "BRCA1", "BRCA2", "MLH1", "MLH3", "MSH2", "MSH6", "FANCA",
                                      "FANCD2", "LIG4", "XRCC4", "XRCC6", "POLE", "POLB", "REV3L", "ATM",
                                      "ATR", "CHEK2", "TP53BP1", "BRIP1", "RAD51", "TOP3A", "ALKBH3", "MGMT")]
all_selected_data <- cbind(DDRscores, selected_genes)

all_selected_data <- inner_join(all_selected_data, TCGA_IDH2, by = "sample")

modified_data <- all_selected_data %>% 
  mutate(mut_class = ifelse(str_detect(IDH1, "1"), "IDH1/2", ifelse(
    str_detect(IDH2, "1"), "IDH1/2", ifelse(
      str_detect(BRCA1, "1"), "BRCA1/2", ifelse(
        str_detect(BRCA2, "1"), "BRCA1/2", ifelse(
          str_detect(MLH1, "1"), "DDR, Other", ifelse(
            str_detect(MLH3, "1"), "DDR, Other", ifelse(
              str_detect(MSH2, "1"), "DDR, Other", ifelse(
                str_detect(MSH6, "1"), "DDR, Other", ifelse(
                  str_detect(FANCA, "1"), "DDR, Other", ifelse(
                    str_detect(FANCD2, "1"), "DDR, Other", ifelse(
                      str_detect(LIG4, "1"), "DDR, Other", ifelse(
                        str_detect(XRCC4, "1"), "DDR, Other", ifelse(
                          str_detect(XRCC6, "1"), "DDR, Other", ifelse(
                            str_detect(POLE, "1"), "DDR, Other", ifelse(
                              str_detect(POLB, "1"), "DDR, Other", ifelse(
                                str_detect(REV3L, "1"), "DDR, Other", ifelse(
                                  str_detect(ATM, "1"), "DDR, Other", ifelse(
                                    str_detect(ATR, "1"), "DDR, Other", ifelse(
                                      str_detect(CHEK2, "1"), "DDR, Other", ifelse(
                                        str_detect(TP53BP1, "1"), "DDR, Other", ifelse(
                                          str_detect(BRIP1, "1"), "DDR, Other", ifelse(
                                            str_detect(RAD51, "1"), "DDR, Other", ifelse(
                                              str_detect(TOP3A, "1"), "DDR, Other", ifelse(
                                                str_detect(ALKBH3, "1"), "DDR, Other", ifelse(
                                                  str_detect(MGMT, "1"), "DDR, Other", "Not DDR/IDH"))))))))))))))))))))))))))

modified_data$mut_class <- factor(modified_data$mut_class, levels = c("BRCA1/2", "DDR, Other", "IDH1/2", "Not DDR/IDH"))
df_TCGA <- modified_data %>% as_tibble() %>% select(1:46, 72) %>%
  mutate(disease_long = recode(tumor_type,
                               `ACC` = "Adrenocortical carcinoma",
                               `BLCA` = "Bladder Uro. Carcinoma",
                               `BRCA` = "Breast invasive carcinoma",
                               `CESC` = "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
                               `CHOL` = "Cholangiocarcinoma",
                               `COAD` = "Colon adenocarcinoma",
                               `DLBC` = "DLBCL",
                               `ESCA` = "Esophageal carcinoma",
                               `GBM` = "Glioblastoma multiforme",
                               `LGG` = "Brain Lower Grade Glioma",
                               `HNSC` = "HNSCC",
                               `KICH` = "Kidney Chromophobe",
                               `KIRC` = "Kidney RCC",
                               `KIRP` = "Kidney PCC",
                               `LAML` = "Acute Myeloid Leukemia",
                               `LIHC` = "HCC",
                               `LUAD` = "Lung adenocarcinoma",
                               `LUSC` = "Lung squa. cell carcinoma",
                               `MESO` = "Mesothelioma",
                               `OV` = "Ovarian serous cystadenocarcinoma",
                               `PAAD` = "Pancreatic adenocarcinoma",
                               `PCPG` = "Paraganglioma/Pheochr.",
                               `PRAD` = "Prostate adenocarcinoma",
                               `READ` = "Rectum adenocarcinoma",
                               `SARC` = "Sarcoma",
                               `SKCM` = "Skin Cutaneous Melanoma",
                               `STAD` = "Stomach adenocarcinoma",
                               `THCA` = "Thyroid carcinoma",
                               `THYM` = "Thymoma",
                               `TGCT` = "Testicular Germ Cell Tumors",
                               `UCEC` = "Endometrial Carcinoma",
                               `UCS` = "Uterine Carcinosarcoma",
                               `UVM` = "Uveal Melanoma"))

### NB need to shorten sample name in TCGA_cosmic
TCGA_cosmic$sample <- TCGA_cosmic$sample %>% str_sub(1, 15)
df_TCGA <- left_join(df_TCGA, TCGA_cosmic, by = "sample")



modified_data_LST <- inner_join(modified_data, TCGA_LST[c(1,3)], by = "sample") %>% select(sample, mut_class, LST)

Darjeeling2_manual_dark <- c("#6B7A80", "#014867", "#997037", "#998472")

plot2_tcga <- function(var, var_title, na.rm = TRUE, ...) {
  data_f <- df_TCGA
  data_WT <- df_TCGA %>% filter(.$mut_class == "Not DDR/IDH")
  upper_limit <- data_f %>% summarise(max = max(.data[[var]], na.rm = TRUE)) %>% as.double()
  #upper_limit = max(data_f${{var}}, na.rm = TRUE) * 1
  print(upper_limit)
  stat_box_data <- function(x) {
    return(data.frame(y = 1.38 * upper_limit,
                      label = paste('N =', length({{x}}))))
  }
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  ggplot(data_f, mapping = aes(x = mut_class, y = .data[[var]], fill = mut_class)) +
    geom_point(aes(color = mut_class), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = mut_class), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = mut_class), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2"), tip.length = 0.01), size = 3) +
    labs(title = paste0("TCGA - ", var_title)) +
    ylab(var_title) +
    xlab(element_blank()) +
    #facet_wrap( ~ disease_long) +
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
  #ggsave(file = paste0("20211227 TCGA ", var_title, ".png"), width = 7, height = 5)
}

### test plot function
plot2_tcga("HRD_LST", "HRD_LST")

plot2a_tcga <- function(var, var_title, na.rm = TRUE, ...) {
  data_f <- as.data.frame(df_TCGA)
  data_f$yvar <- data_f[,var]
  data_WT <- df_TCGA %>% filter(.$mut_class == "Not DDR/IDH")
  upper_limit = max(data_f$yvar, na.rm = TRUE) * 1
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length(x))))
  }
  dunn_stat <- data_f %>% dunn_test(yvar ~ mut_class, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "mut_class", step.increase = 0.1)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = TRUE, digits = 3)
  print(dunn_stat)
  kruskal_test <- compare_means(yvar ~ mut_class, data = data_f, method = "kruskal")
  print(kruskal_test)
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  ggplot(data_f, mapping = aes(x = mut_class, y = yvar, fill = mut_class)) +
    geom_point(aes(color = mut_class), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = mut_class), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = mut_class), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_pvalue_manual(dunn_stat, label = "p.adj.sci", hide.ns = TRUE, inherit.aes = FALSE,
                       size = 3, tip.length = 0.01) +
    #stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2"), tip.length = 0.01), size = 3) +
    labs(title = paste0("TCGA - ", var_title)) +
    ylab(var_title) +
    xlab(element_blank()) +
    #facet_wrap( ~ disease_long) +
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
  #ggsave(file = paste0("20211227 TCGA ", var_title, ".png"), width = 7, height = 5)
}

### test plot function
plot2a_tcga("HRD_LST", "HRD_LST")

plot2b_tcga <- function(var, var_title, na.rm = TRUE, ...) {
  data_f <- as.data.frame(df_TCGA)
  data_f$yvar <- data_f[,var]
  data_WT <- df_TCGA %>% filter(.$mut_class == "Not DDR/IDH")
  upper_limit = max(data_f$yvar, na.rm = TRUE) * 1
  stat_box_data <- function(x) {
    return(data.frame(y = 1.3 * upper_limit,
                      label = paste('N =', length(x))))
  }
  dunn_stat <- data_f %>% dunn_test(yvar ~ mut_class, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "mut_class", step.increase = 0.075)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = TRUE, digits = 3)
  print(dunn_stat)
  kruskal_test <- compare_means(yvar ~ mut_class, data = data_f, method = "kruskal")
  print(kruskal_test)
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  ggplot(data_f, mapping = aes(x = mut_class, y = yvar, fill = mut_class)) +
    geom_point(aes(color = mut_class), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = mut_class), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = mut_class), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_pvalue_manual(dunn_stat, label = "p.adj.sci", hide.ns = TRUE, inherit.aes = FALSE,
                       size = 2.5, tip.length = 0.01) +
    #stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2"), tip.length = 0.01), size = 3) +
    labs(title = paste0("TCGA - ", var_title)) +
    ylab(var_title) +
    xlab(element_blank()) +
    #facet_wrap( ~ disease_long) +
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
  #ggsave(file = paste0("20211227 TCGA ", var_title, ".png"), width = 7, height = 5)
}

### test plot function
plot2b_tcga("HRD_LST", "HRD_LST")

plot2c_tcga <- function(var, var_title, na.rm = TRUE, ...) {
  data_f <- as.data.frame(df_TCGA)
  data_f$yvar <- data_f[,var]
  data_WT <- df_TCGA %>% filter(.$mut_class == "Not DDR/IDH")
  upper_limit = max(data_f$yvar, na.rm = TRUE) * 1
  stat_box_data <- function(x) {
    return(data.frame(y = 1.4 * upper_limit,
                      label = paste('N =', length(x))))
  }
  dunn_stat <- data_f %>% dunn_test(yvar ~ mut_class, p.adjust.method = "hochberg")
  dunn_stat <- dunn_stat %>% filter(group1 == "IDH1/2" | group2 == "IDH1/2")
  dunn_stat <- dunn_stat %>% add_xy_position(x = "mut_class", step.increase = 0.05)
  dunn_stat$p.adj.sci <- format(dunn_stat$p.adj, scientific = TRUE, digits = 3)
  print(dunn_stat)
  kruskal_test <- compare_means(yvar ~ mut_class, data = data_f, method = "kruskal")
  print(kruskal_test)
  #compare_means(HRD_Score ~ mut_class, data = data_f, ref.group = "IDH1", method = "t.test")
  ggplot(data_f, mapping = aes(x = mut_class, y = yvar, fill = mut_class)) +
    geom_point(aes(color = mut_class), size = 2.5, alpha = 1/10, position = "jitter") +
    scale_color_manual(values = rev(Darjeeling2_manual)) +
    new_scale_color() +
    geom_boxplot(aes(color = mut_class), outlier.shape = NA, alpha = 1/5) +
    geom_violin(aes(color = mut_class), alpha = 1/5) +
    # geom_violin(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", alpha = 1/5) +
    # geom_boxplot(data = df_WT, aes(x = GROUP, y = LST), color = "#3f5151", outlier.shape = NA, alpha = 1/5) +
    scale_color_manual(values = rev(Darjeeling2_manual_dark)) +
    scale_fill_manual(values = rev(Darjeeling2_manual)) +
    stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9), hjust = 0.5, vjust = 0, color = "grey50", size = 3, show.legend = FALSE) +
    stat_pvalue_manual(dunn_stat, label = "p.adj.sci", hide.ns = TRUE, inherit.aes = FALSE,
                       size = 2.5, tip.length = 0.01) +
    #stat_compare_means(comparisons = list(c("BRCA1/2", "IDH1/2"), c("DDR, Other", "IDH1/2"), c("Not DDR/IDH", "IDH1/2"), tip.length = 0.01), size = 3) +
    labs(title = paste0("TCGA - ", var_title)) +
    ylab(var_title) +
    xlab(element_blank()) +
    #facet_wrap( ~ disease_long) +
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
  #ggsave(file = paste0("20211227 TCGA ", var_title, ".png"), width = 7, height = 5)
}

### generate plot list
variables_tcga <- tibble(var = list("HRD_LST", "mutSig1", "mutSig6", "mutSig10", "mutLoad_silent", "mutLoad_nonsilent"),
                         var_title = list("HRD LST Score", "MS1 - AC/T>AN", "MS6 - CpG", "MS10 - Toxin", "Silent Mut. Load", "Non-silent Mut. Load"))
plots_tcga <- variables_tcga %>% pmap(plot2a_tcga)
names(plots_tcga) <- variables_tcga$var %>% as.character()
plots_tcga2 <- variables_tcga %>% pmap(plot2b_tcga)
names(plots_tcga2) <- variables_tcga$var %>% as.character()
plots_tcga3 <- variables_tcga %>% pmap(plot2c_tcga)
names(plots_tcga3) <- variables_tcga$var %>% as.character()

save(plots_tcga,
     plots_tcga2,
     plots_tcga3,
     file = "20220111 TCGA LST + mutsigcv.rdata")


select_COSMIC <- c("SBS1", "SBS2", "SBS3", "SBS11", "SBS13", "SBS18")
plots_wCOSMIC_select <- variables_all_tb %>% filter(var %in% select_COSMIC) %>% pmap(plot2b_tcga)
names(plots_wCOSMIC_select) <- select_COSMIC
plots_wCOSMIC_select$SBS1 + scale_y_break(breaks = c(2500,2500), scales = 0.6)

plots_wCOSMIC_select2 <- variables_all_tb %>% filter(var %in% select_COSMIC) %>% pmap(plot2c_tcga)
names(plots_wCOSMIC_select2) <- select_COSMIC

save(plots_wCOSMIC_select,
     plots_wCOSMIC_select2,
     file = "TCGA_plots_wCOSMIC_select.rdata")