library(survminer)
library(survival)
library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsignif)
library(viridis)
library(cowplot)
library(lemon)
library(ggpubr)

load("/Users/schvartj/Google Drive/LabData/R Projects/IDH_TMB/20220107 IMPACT LST all plots.rdata")
load("/Users/schvartj/Google Drive/LabData/R Projects/IDH_TMB/20220107 IMPACT TMB all plots.rdata")
load("/Users/schvartj/Google Drive/LabData/R Projects/IDH_TMB/20220111 TCGA LST + mutsigcv.rdata")
load("/Users/schvartj/Google Drive/LabData/cbioportal/20220106 tumor frequencies plot brca v idh.rdata")
load("/Users/schvartj/Google Drive/LabData/cbioportal/20220106 IMPACT tumor_group survival plots.rdata")
load("/Users/schvartj/Google Drive/LabData/cbioportal/20220106 BRCAvIDHvUnaltered survival plot all.rdata")
load("/Users/schvartj/Google Drive/LabData/cbioportal/20220106 BRCAvIDHvUnaltered survival plot all no brain.rdata")

plots_tcga_LST <- plots_tcga$HRD_LST + labs(title = "LST by Mutation - All Tumors - TCGA Data", y = "Large Scale Transition Score")
fig1_panel_AB <- plot_grid(IMPACT_LST_all_tumors, plots_tcga_LST,
                           labels = c("A", "B"),
                           label_size = 16,
                           nrow = 1)
ggsave(plot = fig1_panel_AB, "fig1_panels_AB.png", height = 5, width = 10)

fig1_panel_C <- tumor_frequencies_plot + scale_size_area(max_size = 7) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(size = 0.5, color = "grey90"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        axis.title = element_text(face = "bold"),
  )

fig1_panel_D <- plot_grid(plot_all[[1]], plot_all[[2]],
                          rel_heights = c(0.75, 0.25),
                          ncol = 1,
                          align = "hv",
                          axis = "l")

fig1_panel_CD <- plot_grid(fig1_panel_C, fig1_panel_D,
                           nrow = 1,
                           align = "h",
                           axis = "b",
                           rel_widths = c(0.3, 0.7),
                           labels = c("C", "D"),
                           label_size = 16)

fig1_final <- plot_grid(fig1_panel_AB,
                        fig1_panel_CD,
                        rel_heights = c(0.35, 0.65),
                        ncol = 1)


##### Supplementary figure 1
names(not_glioma) <- c("not_glioma")
list_by_tumor_type <- append(list_by_tumor_type, not_glioma)
sFig1_panels_A_H_list <- list_by_tumor_type[c("not_glioma", "GLIOMA", "MELANOMA", "BREAST CANCER", "NON-SMALL CELL LUNG CANCER", "HEPATOBILIARY CANCER", "PROSTATE CANCER", "COLORECTAL CANCER")]
sFig1_panels_A_H <- plot_grid(plotlist = sFig1_panels_A_H_list,
                             nrow = 2,
                             label_size = 16,
                             labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
ggsave(plot = sFig1_panels_A_H, "sFig1_panels_A_H.png", height = 6, width = 12)











