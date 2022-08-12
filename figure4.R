library(CytoExploreR)
library(CytoExploreRData)
library(tidyverse)
library(viridis)
library(wesanderson)
library(ggtext)
library(cowplot)
library(ggsignif)
library(ggpubr)
library(scales)
library(ggcyto)
library(lemon)
library(ggridges)

load("/Users/schvartj/Google Drive/LabData/H2AX/20200127 10T EdU dT block release/20200127_10T_EdU_async_for_figure.rdata")
load("/Users/schvartj/Google Drive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20200127_U2OS_EdU_async_for_figure.rdata")
load("/Users/schvartj/Google Drive/LabData/H2AX/20211123 U2OS sgEHMT1 sgCBX1,5 pH2AX/20211123 U2OS IDH2 by guide.rdata")
load("/Users/schvartj/Google Drive/LabData/Flow Cytometry/20211222 U2OS sgCBX,2f,EHMT_002/20211222 U2OS olaparib + guides.rdata")
load("/Users/schvartj/Google Drive/LabData/H2AX/20210917 IF U2OS CldU for ssDNA/20210917 pH2AX IF U2OS olaparib.rdata")
load("/Users/schvartj/Google Drive/LabData/Flow Cytometry/20211222 U2OS sgCBX,2f,EHMT_002/20211222 U2OS olaparib + guides.rdata")
load("/Users/schvartj/Google Drive/LabData/H2AX/20211123 U2OS sgEHMT1 sgCBX1,5 pH2AX/20211123 U2OS PAR IF.rdata")

f4a <- c10T_H2AX_IF
f4b_top <- c10T_H2AX_IF_by_DNA
f4b_bot <- c10T_DNA
f4d <- U2OS_H2AX_IF
f4e_top <- U2OS_H2AX_IF_by_DNA
f4e_topb <- U2OS_H2AX_IF_by_DNAb
f4e_bot <- U2OS_DNA

IF_c10T_WT <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20200127 10T EdU dT block release/IDH2-WT (RGB) example.png") +
  cowplot::draw_label("IDH2-WT", x = 0.01, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_c10T_R172K <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20200127 10T EdU dT block release/IDH2-R172K (RGB) example.png") +
  cowplot::draw_label("IDH2-R172K", x = 0.01, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
  
f4b <- plot_grid(f4b_top, f4b_bot,
                 ncol = 1,
                 rel_heights = c(0.75, 0.25),
                 align = "hv",
                 axis = "tblr")

fig4_panels_ab <- plot_grid(f4a, f4b,
                            nrow = 1,
                            labels = c("A", "B"),
                            label_size = 16)

fig4_panel_c <- plot_grid(NULL, IF_c10T_WT, NULL, IF_c10T_R172K,
                          nrow = 1,
                          rel_widths = c(0.05, 0.47, 0.02, 0.47),
                          labels = "C",
                          label_size = 16)

fig4_panels_abc <- plot_grid(fig4_panels_ab, fig4_panel_c,
                             nrow = 1,
                             rel_widths = c(0.5, 0.5))

IF_U2OS_WT <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/Annotation 8 (RGB).png") +
  cowplot::draw_label("IDH2-WT", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_U2OS_R172K <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/Annotation 4 (RGB).png") +
  cowplot::draw_label("IDH2-R172K", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")

IF_U2OS_WTz <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20210901 slide3 WT 1-4-Airyscan Processing-10-Scene-01 zoom.png") +
  cowplot::draw_label("IDH2-WT", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_U2OS_R172Kz <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/schvartj/Google Drive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20210901 slide3 R172K 5-8-Airyscan Processing-04-Scene-01 zoom.png") +
  cowplot::draw_label("IDH2-R172K", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")

f4e <- plot_grid(f4e_topb, f4e_bot,
                 ncol = 1,
                 rel_heights = c(0.75, 0.25),
                 align = "hv",
                 axis = "tblr")

fig4_panels_de <- plot_grid(f4d, f4e,
                            nrow = 1,
                            labels = c("D", "E"),
                            label_size = 16)

fig4_panel_f <- plot_grid(NULL, IF_U2OS_WT, NULL, IF_U2OS_R172K,
                          nrow = 1,
                          rel_widths = c(0.05, 0.47, 0.02, 0.47),
                          labels = "F",
                          label_size = 16)

fig4_panel_fz <- plot_grid(NULL, IF_U2OS_WTz, NULL, IF_U2OS_R172Kz,
                          nrow = 1,
                          rel_widths = c(0.05, 0.47, 0.02, 0.47),
                          labels = "F",
                          label_size = 16)

fig4_panels_def <- plot_grid(fig4_panels_de, fig4_panel_f,
                             nrow = 1,
                             rel_widths = c(0.5, 0.5))

fig4_panels_defz <- plot_grid(fig4_panels_de, fig4_panel_fz,
                             nrow = 1,
                             rel_widths = c(0.5, 0.5))

fig4_panels_g <- plot_grid(U2OS_IDH2_by_guide_bis,
                           labels = "G",
                           label_size = 16)

fig4_panels_ag <- plot_grid(fig4_panels_abc, NULL, fig4_panels_def, NULL, fig4_panels_g,
                            ncol = 1,
                            rel_heights = c(0.32, 0.01, 0.32, 0.01, 0.32))

fig4_panels_gh <- plot_grid(fig4_panels_g, NULL, U2OS_H2AX_IF_olaparib_by_IDH2_bis,
                            nrow = 1,
                            rel_widths = c(0.49, 0.02, 0.49),
                            labels = c("", "", "H"),
                            label_size = 16)

fig4_panels_ijk <- plot_grid(U2OS_PAR_IF, NULL, U2OS_RPA70_IF_by_condition, NULL, U2OS_olaparib_guides_R26.1_only,
                             nrow = 1,
                             rel_widths = c(0.33, 0.01, 0.33, 0.01, 0.33),
                             labels = c("I", "", "J", "", "K"),
                             label_size = 16)

fig4_panels_ak <- plot_grid(fig4_panels_abc, NULL,
                            fig4_panels_defz, NULL,
                            fig4_panels_gh, NULL,
                            fig4_panels_ijk,
                            ncol = 1,
                            rel_heights = c(0.20, 0.01, 0.20, 0.01, 0.26, 0.01, 0.26))

ggsave(plot = fig4_panels_ak, "figure4.png", height = 14, width = 12)






