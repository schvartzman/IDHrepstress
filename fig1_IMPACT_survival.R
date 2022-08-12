library(survminer)
library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsignif)
library(viridis)
library(cowplot)
library(lemon)
library(survival)
library(vtree)
library(ggvenn)
library(survival)

BRCA1 <- read.delim("20220105_IMPACT_data/brca1.tsv") %>% mutate(gene = ifelse(str_detect(BRCA1, "no alteration"), "", 'BRCA1')) %>% as_tibble() %>%
  filter(!duplicated(Patient.ID))
names(BRCA1)[names(BRCA1) == "Altered"] <- "altered_brca1"

BRCA2 <- read.delim("20220105_IMPACT_data/BRCA2.tsv") %>% mutate(gene = ifelse(str_detect(BRCA2, "no alteration"), "", 'BRCA2')) %>% as_tibble() %>%
  filter(!duplicated(Patient.ID))
names(BRCA2)[names(BRCA2) == "Altered"] <- "altered_brca2"

IDH1 <- read.delim("20220105_IMPACT_data/IDH1.tsv") %>% mutate(gene = ifelse(str_detect(IDH1, "no alteration"), "", 'IDH1')) %>% as_tibble() %>%
  filter(!duplicated(Patient.ID))
names(IDH1)[names(IDH1) == "Altered"] <- "altered_idh1"

IDH2 <- read.delim("20220105_IMPACT_data/IDH2.tsv") %>% mutate(gene = ifelse(str_detect(IDH2, "no alteration"), "", 'IDH2')) %>% as_tibble() %>%
  filter(!duplicated(Patient.ID))
names(IDH2)[names(IDH2) == "Altered"] <- "altered_idh2"

BRCA1.2 <- left_join(BRCA1[c(3,4,5,10)], BRCA2[c(3,4,5,10)], by = "Patient.ID")
IDH1.2 <- left_join(IDH1[c(3,4,5,10)], IDH2[c(3,4,5,10)], by = "Patient.ID")
all_genes <- left_join(BRCA1.2, IDH1.2, by = "Patient.ID")

all_genes <- all_genes %>% mutate(mut_sum = altered_brca1 + altered_brca2 + altered_idh1 + altered_idh2) %>% 
  mutate(gene = ifelse(altered_brca1 + altered_brca2 == "2" & mut_sum == 2, "BRCA1+2", ifelse(
    altered_idh1 + altered_idh2 == "2" & mut_sum == 2, "IDH1+IDH2", ifelse(
      altered_brca1 == "1" & mut_sum == 1, "BRCA1", ifelse(
        altered_brca2 == "1" & mut_sum == 1, "BRCA2", ifelse(
          altered_idh1 == "1" & mut_sum == 1, "IDH1", ifelse(
            altered_idh2 == "1" & mut_sum == 1, "IDH2", ifelse(
              mut_sum == 0, "unaltered", "mixed")))))))) %>%
  mutate(gene_group = ifelse(gene %in% c("BRCA1", "BRCA2", "BRCA1+2"), "BRCA1/2", ifelse(
    gene %in% c("IDH1", "IDH2", "IDH1+IDH2"), "IDH1/2", ifelse(
      gene == "mixed", "mixed", "unaltered")))) %>%
  select(-4, -7, -10, -13)

all_genes %>% group_by(gene) %>% count()
all_genes %>% group_by(gene_group) %>% count()

data_clin <- read_tsv("20220105_IMPACT_data/mskimpact_clinical_data.tsv") %>% as_tibble() %>% select(2, 5, 6, 8, 9, 41, 43, 57)
data_clin_colnames <- c("Patient.ID", "age_at_seq", "age_current", "cancer_type", "cancer_type_det", "OS_months", "OS_status", "sex")
colnames(data_clin) <-  data_clin_colnames
data_clin <- data_clin %>% filter(!duplicated(Patient.ID))

data_all <- left_join(data_clin, all_genes, by = "Patient.ID") %>% filter(!is.na(OS_status))

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

tumor_types <- unique(data_all$cancer_type)
data_all$OS_status <- data_all$OS_status %>% recode("1:DECEASED" = 1, "0:LIVING" = 0)
data_all <- data_all %>% mutate(tumor_group = ifelse(cancer_type %in% GI, "GI", ifelse(
  cancer_type %in% Sarcoma, "Sarcoma", ifelse(
        cancer_type %in% SCC, "SCC", ifelse(
            cancer_type %in% Brain, "Brain", ifelse(
              cancer_type %in% Hematopoietic, "Hematopoietic", ''))))))
data_all$cancer_type <- data_all$cancer_type %>% recode("Gastrointestinal Neuroendocrine Tumors of the Esophagus/Stomach" = "GI NETs of the Esophagus/Stomach",
                                            "Myelodysplastic/Myeloproliferative Neoplasms" = "MDS/MPNs",
                                            "Gastrointestinal Neuroendocrine Tumor" = "GI NETs",
                                            "B-Lymphoblastic Leukemia/Lymphoma" = "B-cell Leukemia/Lymphoma",
                                            "T-Lymphoblastic Leukemia/Lymphoma" = "T-cell Leukemia/Lymphoma")

brca_types <- data_all %>% filter(gene_group %in% c("BRCA1/2")) %>% group_by(cancer_type) %>% summarise(count = n())
idh_types <- data_all %>% filter(gene_group %in% c("IDH1/2")) %>% group_by(cancer_type) %>% summarise(count = n())
selected_types <- inner_join(brca_types, idh_types, by = "cancer_type") %>% filter(cancer_type != "NA")

selected_types_short <- c("ACC", "Appendiceal", "B-ALL/Lymphoma", "Bladder", "Bone", "Breast", "CUP", "CNS", "Colorectal", "Endometrial",
                          "Germ Cell", "Glioma", "HNSCC", "HPB", "Histiocytosis", "Leukemia", "Mastocytosis", "B-ALL/Lymphoma", "T/NK Neoplasm", "Melanoma",
                          "MDS", "MDS", "MPN", "MPN", "Neopastic/Reactive", "NSCLC", "PDAC", "Prostate", "RCC", "Non-Melanoma Skin", "Small Bowel",
                          "STS", "Thyroid", "Colorectal")
names(selected_types_short) <- selected_types$cancer_type

data_all_sub <- data_all %>% filter(cancer_type %in% selected_types$cancer_type) %>% mutate(cancer_type = recode(cancer_type, !!!selected_types_short))

colors <- c("#49cf8a", "#785698", "#49cf8a", "#785698")
colsIDH_2old <- c("#ECCBAE", "#046C9A")
colsIDH_2 <- c("#e1be6d", "#35274a")
colsIDH <- c("#899da4", "#482576FF", "#43BF71FF")
colsIDH_4 <- viridis_pal(option = "D")(4)
colsIDH3 <- c("#e1be6d", "#35274a", "#899da4")
colsIDH3old <- c("#ECCBAE", "#046C9A", "#899da4")
c("#899da4", "#046C9A", "#D69C4E", "#ECCBAE")

tumor_frequencies_plot <- data_all_sub %>% filter(gene_group %in% c("BRCA1/2", "IDH1/2") & cancer_type != "NA") %>%
  ggplot(mapping = aes(x = gene_group, y = fct_rev(fct_infreq(cancer_type)), color = gene_group)) +
  scale_color_manual(values = colsIDH_2old) +
  scale_fill_manual(values = colsIDH_2old) +
  geom_count() + #aes(size = after_stat(prop), group = gene)) +
  scale_size_area(max_size = 6) +
  scale_y_discrete(expand = expansion(add = c(1,1))) +
  labs(color = "Mutation", size = "Number", x = element_blank(), y = element_blank(),
       title = "Sample Frequencies by\nTumor Type and Mutation Group") +
  guides(color = "none") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey90"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        #legend.position = c(0.1, 0.15)
        )
ggsave(filename = "20220105 all tumor types brca vs IDH.png", height = 6, width = 3.5)
save(tumor_frequencies_plot, file = "20220106 tumor frequencies plot brca v idh.rdata")

theme_survminer2 <- function (base_size = 10, base_family = "",
                              font.main = c(14, "bold", "black"),
                              font.submain = c(10, "plain", "black"), 
                              font.x = c(11, "bold", "black"),
                              font.y = c(11, "bold", "black"),
                              font.caption = c(11, "plain", "black"),
                              font.tickslab = c(9, "plain", "black"),
                              legend = c("top", "bottom", "left", "right", "none"),
                              font.legend = c(10, "plain", "black"),
                              base_line_size = base_size/22,
                              ...) 
{
  font.main <- ggpubr:::.parse_font(font.main)
  font.x <- ggpubr:::.parse_font(font.x)
  font.y <- ggpubr:::.parse_font(font.y)
  font.submain <- ggpubr:::.parse_font(font.submain)
  font.caption <- ggpubr:::.parse_font(font.caption)
  font.tickslab <- ggpubr:::.parse_font(font.tickslab)
  font.legend <- ggpubr:::.parse_font(font.legend)
  if (!is(legend, "numeric")) 
    legend <- match.arg(legend)
  tickslab <- element_text(size = font.tickslab$size, face = font.tickslab$face, 
                           colour = font.tickslab$color, angle = 0)
  legend.text <- element_text(size = font.legend$size, face = font.legend$face, 
                              colour = font.legend$color)
  result <- theme_classic(base_size = base_size, base_family = base_family, base_line_size = base_line_size) + 
    theme(plot.title = element_text(size = font.main$size, 
                                    lineheight = 1, face = font.main$face, colour = font.main$color), 
          plot.subtitle = element_text(size = font.submain$size, 
                                       lineheight = 1, face = font.submain$face, colour = font.submain$color), 
          axis.title.x = element_text(size = font.x$size, face = font.x$face, 
                                      colour = font.x$color), axis.title.y = element_text(angle = 90, 
                                                                                          size = font.y$size, face = font.y$face, colour = font.y$color), 
          axis.line = element_line(size = rel(1), lineend = "square"),
          plot.caption = element_text(size = font.caption$size, 
                                      lineheight = 1, face = font.caption$face, colour = font.caption$color), 
          axis.text.x = tickslab, axis.text.y = tickslab, legend.position = legend, 
          legend.text = legend.text, legend.title = legend.text)
  class(result) <- "theme"
  result
}


#survival plot all -liquid
require("survival")
data_select <- data_all %>% filter(gene_group != "mixed" & tumor_group != "Hematopoietic")
fit_all <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
plot_all <- ggsurvplot(fit_all, data = data_select,
                       risk.table = TRUE,
                       #tables.theme = theme_cleantable(),
                       tables.y.text = FALSE,
                       pval = "Log-rank mIDH1/2 vs BRCA1/2:\np = 0.0014",
                       pval.method = TRUE,
                       palette = colsIDH3old,
                       axes.offset = TRUE,
                       conf.int = TRUE,
                       conf.int.alpha = 1/5,
                       conf.int.style = "ribbon",
                       censor.size = 3,
                       censor.shape = 124,
                       surv.median.line = c("v"),
                       #size = 1.5,
                       pval.size = 4,
                       pval.coord = c(0.1, 0.1),
                       pval.method.size = 4,
                       pval.method.coord = c(0.1, 0.165),
                       fontsize = 3.5,
                       legend.labs = c("BRCA1/2", "mIDH1/2", "BRCA/IDH WT"),
                       legend.title = "Alteration",
                       legend = c(0.9, 0.85),
                       xlab = "Time (months)",
                       title = "Survival for All Solid Tumors",
                       font.title = "bold",
                       ggtheme = theme_survminer2())
save(plot_all, file = "20220106 BRCAvIDHvUnaltered survival plot all.rdata")

plot_grid(plotlist = list(plot_all[[1]], plot_all[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
ggsave(filename = "20220105 BRCAvIDHvUNALTERED survival plot all.png", height = 6, width = 8)

require("survival")
data_select <- data_all %>% filter(gene_group %in% c("BRCA1/2", "IDH1/2") & tumor_group != "Hematopoietic")
fit_all <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
plot_all2 <- ggsurvplot(fit_all, data = data_select,
                       risk.table = TRUE,
                       #tables.theme = theme_cleantable(),
                       tables.y.text = FALSE,
                       pval = TRUE,
                       pval.method = TRUE,
                       palette = colsIDH3old,
                       axes.offset = TRUE,
                       conf.int = TRUE,
                       conf.int.alpha = 1/5,
                       conf.int.style = "ribbon",
                       censor.size = 3,
                       censor.shape = 124,
                       surv.median.line = c("v"),
                       #size = 1.5,
                       pval.size = 4,
                       pval.coord = c(0.1, 0.1),
                       pval.method.size = 4,
                       pval.method.coord = c(0.1, 0.165),
                       fontsize = 3.5,
                       legend.labs = c("BRCA1/2", "mIDH1/2"),
                       legend.title = "Alteration",
                       legend = c(0.9, 0.85),
                       xlab = "Time (months)",
                       title = "Survival for All Solid Tumors",
                       font.title = "bold",
                       ggtheme = theme_survminer2())

#survival plot all -liquid -brain
require("survival")
data_select_nobrain <- data_all %>% filter(gene_group != "mixed" & tumor_group != "Hematopoietic" & tumor_group != "Brain")
fit_all_nobrain <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select_nobrain)
plot_all_nobrain <- ggsurvplot(fit_all_nobrain, data = data_select_nobrain,
                       risk.table = TRUE,
                       pval = TRUE,
                       pval.method = TRUE,
                       palette = colsIDH3old,
                       tables.y.text = FALSE,
                       axes.offset = TRUE,
                       conf.int = TRUE,
                       conf.int.alpha = 1/5,
                       conf.int.style = "ribbon",
                       censor.size = 3,
                       censor.shape = 124,
                       surv.median.line = c("v"),
                       pval.size = 4,
                       pval.coord = c(0.1, 0.1),
                       pval.method.size = 4,
                       pval.method.coord = c(0.1, 0.25),
                       fontsize = 3.5,
                       legend.labs = c("BRCA1/2", "mIDH1/2", "unaltered"),
                       legend.title = "Alteration",
                       legend = c(0.9, 0.9),
                       xlab = "Time (months)", #,
                       title = "All Tumors Except Glioma",
                       font.title = "bold",
                       ggtheme = theme_survminer2())

save(plot_all_nobrain, file = "20220106 BRCAvIDHvUnaltered survival plot all no brain.rdata")

require("survival")
data_select_nobrain <- data_all %>% filter(gene_group %in% c("BRCA1/2", "IDH1/2") & tumor_group != "Hematopoietic" & tumor_group != "Brain")
fit_all_nobrain <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select_nobrain)
plot_all_nobrain2 <- ggsurvplot(fit_all_nobrain, data = data_select_nobrain,
                               risk.table = TRUE,
                               pval = TRUE,
                               pval.method = TRUE,
                               palette = colsIDH3old,
                               tables.y.text = FALSE,
                               axes.offset = TRUE,
                               conf.int = TRUE,
                               conf.int.alpha = 1/5,
                               conf.int.style = "ribbon",
                               censor.size = 3,
                               censor.shape = 124,
                               surv.median.line = c("v"),
                               pval.size = 4,
                               pval.coord = c(0.1, 0.1),
                               pval.method.size = 4,
                               pval.method.coord = c(0.1, 0.25),
                               fontsize = 3.5,
                               legend.labs = c("BRCA1/2", "mIDH1/2"),
                               legend.title = "Alteration",
                               legend = c(0.9, 0.9),
                               xlab = "Time (months)",
                               ylab = "Survival Prob.",
                               title = "All Tumors Except Glioma",
                               font.title = "bold",
                               ggtheme = theme_survminer2())





plot_grid(plotlist = list(plot_all_nobrain[[1]], plot_all_nobrain[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
ggsave(filename = "20220105 BRCAvIDHvUNALTERED survival plot all no brain.png", height = 6, width = 8)

survival_plot <- function(cancer_type, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group != "mixed" & .$cancer_type == {{cancer_type}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = TRUE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2", "unaltered"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            xlab = "Time (months)",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
# plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
# ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}

survival_plot2 <- function(cancer_type, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group %in% c("BRCA1/2", "IDH1/2") & .$cancer_type == {{cancer_type}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = TRUE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            xlab = "Time (months)",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
  # plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
  # ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}


survival_plot_small <- function(cancer_type, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group != "mixed" & .$cancer_type == {{cancer_type}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = FALSE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2", "unaltered"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            xlab = "Time (months)",
                            ylab = "Survival Prob.",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
  # plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
  # ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}

survival_plot_small2 <- function(cancer_type, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group %in% c("BRCA1/2", "IDH1/2") & .$cancer_type == {{cancer_type}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = FALSE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            xlab = "Time (months)",
                            ylab = "Survival Prob.",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
  # plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
  # ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}

survival_plot_tumorgroup_small <- function(tumor_group, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group != "mixed" & .$tumor_group == {{tumor_group}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = FALSE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2", "unaltered"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            #legend = "none",
                            xlab = "Time (months)",
                            ylab = "Survival Prob.",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
  # plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
  # ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}

survival_plot_tumorgroup_small2 <- function(tumor_group, title) {
  require("survival")
  data_select <- data_all %>% filter(.$gene_group %in% c("BRCA1/2", "IDH1/2") & .$tumor_group == {{tumor_group}})
  fit_select <- survfit(Surv(OS_months, OS_status) ~ gene_group, data = data_select)
  plot_select <- ggsurvplot(fit_select, data = data_select,
                            risk.table = FALSE,
                            risk.table.y.text = FALSE,
                            pval = TRUE,
                            pval.method = TRUE,
                            palette = colsIDH3old,
                            axes.offset = TRUE,
                            conf.int = TRUE,
                            conf.int.alpha = 1/5,
                            conf.int.style = "ribbon",
                            censor.size = 3,
                            censor.shape = 124,
                            surv.median.line = c("v"),
                            pval.size = 4,
                            pval.coord = c(0.1, 0.1),
                            pval.method.size = 4,
                            pval.method.coord = c(0.1, 0.25),
                            fontsize = 3.5,
                            legend.labs = c("BRCA1/2", "mIDH1/2"),
                            legend.title = "Alteration",
                            legend = c(0.9, 0.9),
                            #legend = "none",
                            xlab = "Time (months)",
                            ylab = "Survival Prob.",
                            title = paste0({{title}}),
                            font.title = "bold",
                            ggtheme = theme_survminer2())
  # plot_grid(plotlist = list(plot_select[[1]], plot_select[[2]]), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")
  # ggsave(filename = paste0("20220106 BRCAvIDH survival plot ", {{title}}, ".png"), height = 6, width = 8)
}

map2(c("GI"), c("GI Tumors"), survival_plot2)

test_list <- list("Glioma", "Glioma")
test_plots <- test_list %>% pmap(survival_plot_small2)
plot_grid(plotlist = list(test_plots[[1]]$plot, test_plots[[1]]$table), ncol = 1, rel_heights = c(0.75, 0.25), align = "v", axis = "l")

cancer_type_idh <- data_all %>% filter(gene_group == "IDH1/2") %>% count(cancer_type, sort = TRUE) %>% filter(n > 0) %>% select(-n)
cancer_type_brca <- data_all %>% filter(gene_group == "BRCA1/2") %>% count(cancer_type, sort = TRUE) %>% filter(n > 0) %>% select(-n)
cancer_type_common <- inner_join(cancer_type_idh, cancer_type_brca, by = "cancer_type")
cancer_type_common <- cancer_type_common %>% filter(cancer_type != "NA")
variables <- tibble(cancer_type = unlist(cancer_type_common), title = unlist(cancer_type_common))

cancer_type_selected <- pmap(variables[c(1:28),], survival_plot)
names(cancer_type_selected) <- variables[c(1:28),1] %>% as_vector()

cancer_type_selected2 <- pmap(variables[c(1:28),], survival_plot2)
names(cancer_type_selected2) <- variables[c(1:28),1] %>% as_vector()

cancer_type_selected_small <- pmap(variables[c(1:28),], survival_plot_small)
names(cancer_type_selected_small) <- variables[c(1:28),1] %>% as_vector()

cancer_type_selected_small2 <- pmap(variables[c(1:28),], survival_plot_small2)
names(cancer_type_selected_small2) <- variables[c(1:28),1] %>% as_vector()

cancer_type_selected_mp <- cancer_type_selected[c(4,12)]

tumor_group <- data_all$tumor_group %>% unique()
tumor_group <- tumor_group[-1]
title <- c("Gastrointestinal Cancers", "Sarcomas", "Squamous Cell Carcinomas", "Primary Brain", "Hematopoietic Malignancies")
variables_tumor_groups <- tibble(tumor_group = tumor_group, title = title)

tumor_group_plots <- pmap(variables_tumor_groups, survival_plot_tumorgroup_small)
names(tumor_group_plots) <- variables_tumor_groups[c(1:5),1] %>% as_vector()

tumor_group_plots2 <- pmap(variables_tumor_groups, survival_plot_tumorgroup_small2)
names(tumor_group_plots2) <- variables_tumor_groups[c(1:5),1] %>% as_vector()



survival_all <- plot_all
survival_all2 <- plot_all2
survival_all_nobrain <- plot_all_nobrain
survival_all_nobrain2 <- plot_all_nobrain2
survival_cancer_type_selected <- cancer_type_selected
survival_cancer_type_selected2 <- cancer_type_selected2
survival_cancer_type_selected_small <- cancer_type_selected_small
survival_cancer_type_selected_small2 <- cancer_type_selected_small2
survival_tumor_group_plots <- tumor_group_plots
survival_tumor_group_plots2 <- tumor_group_plots2

write_csv(data_all, "data_all_MSK-IMPACT.csv")

save(tumor_frequencies_plot,
     survival_all,
     survival_all2,
     survival_all_nobrain,
     survival_all_nobrain2,
     survival_cancer_type_selected,
     survival_cancer_type_selected2,
     survival_cancer_type_selected_small,
     survival_cancer_type_selected_small2,
     survival_tumor_group_plots,
     survival_tumor_group_plots2,
     file = "20220106 IMPACT survival plots all.rdata")










