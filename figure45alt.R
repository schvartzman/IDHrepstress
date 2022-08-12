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

load("/Users/juanma/GDrive/LabData/H2AX/20200127 10T EdU dT block release/20200127_10T_EdU_async_for_figure.rdata")
load("/Users/juanma/GDrive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20200127_U2OS_EdU_async_for_figure.rdata")
load("/Users/juanma/GDrive/LabData/H2AX/20211123 U2OS sgEHMT1 sgCBX1,5 pH2AX/20211123 U2OS IDH2 by guide.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20211222 U2OS sgCBX,2f,EHMT_002/20211222 U2OS olaparib + guides.rdata")
load("/Users/juanma/GDrive/LabData/H2AX/20210917 IF U2OS CldU for ssDNA/20210917 pH2AX IF U2OS olaparib.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20211222 U2OS sgCBX,2f,EHMT_002/20211222 U2OS olaparib + guides.rdata")
load("/Users/juanma/GDrive/LabData/H2AX/20211123 U2OS sgEHMT1 sgCBX1,5 pH2AX/20211123 U2OS PAR IF.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20220308 U2OS guides +- ola/20220304 U2OS olaparib.rdata")

f5a <- c10T_H2AX_IF
f5b_top <- c10T_H2AX_IF_by_DNA
f5b_bot <- c10T_DNA
f5d <- U2OS_H2AX_IF
f5e_top <- U2OS_H2AX_IF_by_DNA
f5e_topb <- U2OS_H2AX_IF_by_DNAb
f5e_bot <- U2OS_DNA

IF_c10T_WT <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/juanma/GDrive/LabData/H2AX/20200127 10T EdU dT block release/IDH2-WT (RGB) example.png") +
  cowplot::draw_label("IDH2-WT", x = 0.01, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_c10T_R172K <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/juanma/GDrive/LabData/H2AX/20200127 10T EdU dT block release/IDH2-R172K (RGB) example.png") +
  cowplot::draw_label("IDH2-R172K", x = 0.01, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")

f5b <- plot_grid(f5b_top, f5b_bot,
                 ncol = 1,
                 rel_heights = c(0.75, 0.25),
                 align = "hv",
                 axis = "tblr")

fig5_panels_ab <- plot_grid(f5a, f5b,
                            nrow = 1,
                            labels = c("A", "B"),
                            label_size = 16)

####for slide
fig5_panel_abc_nolabs <- plot_grid(f5a, f5b,
                                   nrow = 1)
ggsave(plot = fig5_panel_abc_nolabs, "fig5_panel_abc_nolabs.png", height = 4, width = 6)



fig5_panel_c <- plot_grid(NULL, IF_c10T_WT, NULL, IF_c10T_R172K,
                          nrow = 1,
                          rel_widths = c(0.05, 0.47, 0.02, 0.47),
                          labels = "C",
                          label_size = 16)

fig5_panels_abc <- plot_grid(fig5_panels_ab, fig5_panel_c,
                             ncol = 1,
                             rel_heights = c(0.5, 0.5))
ggsave(plot = fig5_panels_abc, "fig5_panels_abc.png", height = 6, width = 6)

IF_U2OS_WT <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      " ") +
  cowplot::draw_label("IDH2-WT", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_U2OS_R172K <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/juanma/GDrive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/Annotation 4 (RGB).png") +
  cowplot::draw_label("IDH2-R172K", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")

IF_U2OS_WTz <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/juanma/GDrive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20210901 slide3 WT 1-4-Airyscan Processing-10-Scene-01 zoom.png") +
  cowplot::draw_label("IDH2-WT", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")
IF_U2OS_R172Kz <- cowplot::ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  cowplot::draw_image(halign = 0, valign = 0,
                      "/Users/juanma/GDrive/LabData/H2AX/20210907 U2OS pTIDH2 WT|172/20210901 slide3 R172K 5-8-Airyscan Processing-04-Scene-01 zoom.png") +
  cowplot::draw_label("IDH2-R172K", x = 0.02, y = 0.02, hjust = 0, vjust = 0, size = 14, fontface = "bold", color = "white")

f5e <- plot_grid(f5e_topb, f5e_bot,
                 ncol = 1,
                 rel_heights = c(0.75, 0.25),
                 align = "hv",
                 axis = "tblr")

fig5_panels_de <- plot_grid(f5d, f5e,
                            nrow = 1,
                            labels = c("D", "E"),
                            label_size = 16)

###for slide
fig5_panel_def_nolabs <- plot_grid(f5d, f5e,
                                   nrow = 1)
ggsave(plot = fig5_panel_def_nolabs, "fig5_panel_def_nolabs.png", height = 4, width = 6)




fig5_panel_f <- plot_grid(NULL, IF_U2OS_WT, NULL, IF_U2OS_R172K,
                          nrow = 1,
                          rel_widths = c(0.05, 0.47, 0.02, 0.47),
                          labels = "F",
                          label_size = 16)

fig5_panel_fz <- plot_grid(NULL, IF_U2OS_WTz, NULL, IF_U2OS_R172Kz,
                           nrow = 1,
                           rel_widths = c(0.05, 0.47, 0.02, 0.47),
                           labels = "F",
                           label_size = 16)

fig5_panels_def <- plot_grid(fig5_panels_de, fig5_panel_f,
                             ncol = 1,
                             rel_heights = c(0.5, 0.5))
ggsave(plot = fig5_panels_def, "fig5_panels_def.png", height = 6, width = 6)


fig5_panels_defz <- plot_grid(fig5_panels_de, fig5_panel_fz,
                              ncol = 1,
                              rel_heights = c(0.5, 0.5))
ggsave(plot = fig5_panels_defz, "fig5_panels_defz.png", height = 6, width = 6)


fig5_panels_a_f <- plot_grid(fig5_panels_abc, fig5_panels_defz,
                              nrow = 1)

fig5_panels_g <- plot_grid(U2OS_IDH2_by_guide_bis,
                           labels = "G",
                           label_size = 16)

### for slide
fig5_panel_G_nolabs <- plot_grid(U2OS_IDH2_by_guide_bis)
ggsave(plot = fig5_panel_G_nolabs, "fig5_panel_G_nolabs.png", height = 5, width = 8)

fig5_panels_HIJ_nolabs <- plot_grid(U2OS_H2AX_IF_olaparib_by_IDH2_bis, U2OS_PAR_IF, U2OS_RPA70_IF_by_condition,
                               nrow = 1,
                               rel_widths = c(0.5, 0.25, 0.25))
ggsave(plot = fig5_panels_HIJ_nolabs, "fig5_panels_HIJ_nolabs.png", height = 5, width = 12)



fig5_panels_gh <- plot_grid(fig5_panels_g, NULL, U2OS_H2AX_IF_olaparib_by_IDH2_bis,
                            nrow = 1,
                            rel_widths = c(0.49, 0.02, 0.49),
                            labels = c("", "", "H"),
                            label_size = 16)

fig5_panels_ijk <- plot_grid(U2OS_PAR_IF, NULL, U2OS_RPA70_IF_by_condition, NULL, U2OS_olap2_EdU,
                             nrow = 1,
                             rel_widths = c(0.33, 0.01, 0.33, 0.01, 0.33),
                             labels = c("I", "", "J", "", "K"),
                             label_size = 16)

### for slide
U2OS_olap_EdU <- plot_grid(U2OS_olap2_EdU)
ggsave(plot = U2OS_olap_EdU, "U2OS_olap_EdU.png", height = 5, width = 5)


fig5_panels_gk <- plot_grid(fig5_panels_gh, NULL, fig5_panels_ijk,
                            ncol = 1,
                            rel_heights = c(0.49, 0.02, 0.49))
ggsave(plot = fig5_panels_gk, "fig5_panels_gk.png", height = 8, width = 12)


fig5_panels_ak_bis <- plot_grid(fig5_panels_a_f, NULL,
                            fig5_panels_gh, NULL,
                            fig5_panels_ijk,
                            ncol = 1,
                            rel_heights = c(0.4, 0.01, 0.3, 0.01, 0.3))

ggsave(plot = fig5_panels_ak_bis, "figure5_bis.png", height = 14, width = 12)

###for poster figure
figure5forposter <- plot_grid(fig4_panels_C, f5d, f5e, U2OS_IDH2_by_guide_bis, U2OS_olap_EdU,
          nrow = 1,
          rel_widths = c(0.13, 0.13, 0.135, 0.25, 0.13),
          labels = c("A", "B", "C", "D", "E"),
          label_size = 16)

ggsave(plot = figure5forposter, "figure5_forposter.png", height = 4, width = 17)

##### supplementary figure 5

load("/Users/juanma/GDrive/LabData/H2AX/20201005 U2OS and RKO IDH2-WTvR172K/20201005 RKO H2AX.rdata")
load("/Users/juanma/GDrive/LabData/H2AX/20200127 10T EdU dT block release/20200127_10T_EdU_async_for_figure_int.rdata")

sf5b_top <- c10T_H2AX_IF_by_DNA_int
sf5b_bot <- c10T_DNA_int
sf5a <- c10T_H2AX_IF_int

sf5d_top <- RKO_H2AX_IF_by_DNA
sf5d_bot <- RKO_DNA
sf5c <- RKO_H2AX_IF

sf5b <- plot_grid(sf5b_top, sf5b_bot,
                 ncol = 1,
                 rel_heights = c(0.75, 0.25),
                 align = "hv",
                 axis = "tblr")

sf5d <- plot_grid(sf5d_top, sf5d_bot,
                  ncol = 1,
                  rel_heights = c(0.75, 0.25),
                  align = "hv",
                  axis = "tblr")

sfig5_panels_ab <- plot_grid(sf5a, sf5b,
                            nrow = 1,
                            labels = c("A", "B"),
                            label_size = 16)

sfig5_panels_cd <- plot_grid(sf5c, sf5d,
                             nrow = 1,
                             labels = c("C", "D"),
                             label_size = 16)

sfig5_panels_a_e <- plot_grid(sfig5_panels_ab,
                               sfig5_panels_cd,
                                U2OS_IDH2_by_guide,
                                ncol = 1,
                                labels = c("", "", "E"),
                                label_size = 16)

ggsave(plot = sfig5_panels_a_e, "sfigure5.png", height = 9, width = 7)














