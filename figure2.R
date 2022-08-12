library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsignif)
library(viridis)
library(cowplot)
library(lemon)
library(rstatix)
library(gridExtra)
library(ggtext)


file_names_GCMS <- as.list(dir(path = "/Users/juanma/GDrive/LabData/GCMS/20201202 RKO U2OS pTurn/", pattern = "*.rdata"))
file_names_GCMS <- lapply(file_names_GCMS, function(x) paste0("/Users/juanma/GDrive/LabData/GCMS/20201202 RKO U2OS pTurn/", x))
lapply(file_names_GCMS, load, environment())

file_names_u2os_growth <- as.list(dir(path = "/Users/juanma/GDrive/LabData/Growth Curves/20210205 U2OS pTIDH1,2 +dox d1-d5/", pattern = "*.rdata"))
file_names_u2os_growth <- lapply(file_names_u2os_growth, function(x) paste0("/Users/juanma/GDrive/LabData/Growth Curves/20210205 U2OS pTIDH1,2 +dox d1-d5/", x))
lapply(file_names_u2os_growth, load, environment())

load("/Users/juanma/GDrive/LabData/Microscopy/20220113 U2OS FUCCI.4/20220421 U2OS FUCCI14.rdata")
load("/Users/juanma/GDrive/LabData/Microscopy/20220113 U2OS FUCCI.4/20220421 U2OS FUCCI14 individual tracks.rdata")
load("/Users/juanma/GDrive/LabData/Microscopy/20220113 U2OS FUCCI.4/U2OS_fucci_AS_3phases.rdata")

WB_U2OS_IDH1 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_u2os_idh1_wb.png")
WB_U2OS_IDH2 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_u2os_idh2_wb.png")
#WB_RKO_IDH1 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_rko_idh1_wb.png")
WB_RKO_IDH2 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_rko_idh2_wb.png")

fig2_panels_A_C <- plot_grid(WB_U2OS_IDH1, NULL, n2HG_U2OS_IDH1, NULL, U2OS_IDH1_growth, ncol = 1,
                             rel_heights = c(0.2, 0.01, 0.4, 0.01, 0.4), align = "v", axis = "l",
                             label_size = 16,
                             labels = c("A", "", "B", "", "C"))
fig2_panels_D_F <- plot_grid(WB_U2OS_IDH2, NULL, n2HG_U2OS_IDH2, NULL, U2OS_IDH2_growth, ncol = 1,
                             rel_heights = c(0.2, 0.01, 0.4, 0.01, 0.4), align = "v", axis = "l",
                             label_size = 16,
                             labels = c("D", "", "E", "", "F"))
fig2_panels_A_F <- plot_grid(fig2_panels_A_C, NULL, fig2_panels_D_F, ncol = 3,
                             rel_widths = c(0.49, 0.02, 0.49))
ggsave("fig2_panels_A_F.png", height = 8, width = 6)

####for slides
fig2_panels_A_C_nolabs <- plot_grid(WB_U2OS_IDH1, NULL, n2HG_U2OS_IDH1, NULL, U2OS_IDH1_growth, nrow = 1,
                                    rel_widths = c(0.2, 0.03, 0.3, 0.03, 0.3))
ggsave(plot = fig2_panels_A_C_nolabs, "fig2_panels_A_C_nolabs.png", height = 3, width = 9)

fig2_panels_D_F_nolabs <- plot_grid(WB_U2OS_IDH2, NULL, n2HG_U2OS_IDH2, NULL, U2OS_IDH2_growth, nrow = 1,
                                    rel_widths = c(0.3, 0.01, 0.3, 0.01, 0.3))
ggsave(plot = fig2_panels_D_F_nolabs, "fig2_panels_D_F_nolabs.png", height = 3, width = 9)



load(file = "20220120 combined growth IDH2.rdata")
load(file = "20220120 combined growth IDH1.rdata")
load(file = "20220120 combined growth.rdata")

combined_growth_IDH1
combined_growth_IDH2


title_fig2_panel_G <- ggdraw() + draw_label("Relative Cell Number", fontface = "bold", hjust = 0, size = 13, x = 0.05)
fig2_panel_Gnotitle <- combined_growth_IDH1_2
fig2_panel_G <- plot_grid(NULL, title_fig2_panel_G, NULL, fig2_panel_Gnotitle,
                          rel_heights = c(0.04, 0.04, 0.02, 0.90), ncol = 1,
                          label_size = 16, labels = "G")
fig2_panel_G_nolabel <- plot_grid(fig2_panel_Gnotitle, ncol = 1)
ggsave(plot = fig2_panel_G_nolabel, "fig2_panel_G_nolabel.png", height = 2.3, width = 7)



fig2_panelHt <- cowplot::ggdraw(xlim = c(-0.01, 1.01)) +
  cowplot::draw_image("/Users/juanma/GDrive/LabData/Microscopy/20220113 U2OS FUCCI.4/U2OS Fucci example cells.png")

fucci_schematic <- cowplot::ggdraw(xlim = c(-0.01, 1.01)) +
  cowplot::draw_image("/Users/juanma/GDrive/LabData/IDH_repstress/figure2/fucci4_schematic.png")

ggsave(plot = fucci_schematic, "fucci_schematic.png", height = 4, width = 6)

U2OS_fucci_phase7 <- U2OS_fucci_phase7 + theme(plot.title = element_text(size = 11))
U2OS_fucci_G2M <- U2OS_fucci_G2M + theme(plot.title = element_text(size = 11))

U2OS_fucci_AS_Sphase <- U2OS_fucci_AS_Sphase + theme(plot.title = element_text(size = 11), legend.position = "none")
U2OS_fucci_AS_G2Mphase <- U2OS_fucci_AS_G2Mphase + theme(plot.title = element_text(size = 11), legend.position = "none")


fig2_panelHb <- plot_grid(fucci_schematic, NULL, U2OS_fucci_AS_Sphase, U2OS_fucci_AS_G2Mphase,
                          nrow = 1, rel_widths = c(0.6, 0.01, 0.23, 0.23),
                          labels = c("I", "", "J", "K"),
                          label_size = 16)

fig2_panelH <- plot_grid(NULL, fig2_panelHt, fig2_panelHb, ncol = 1, rel_heights = c(0.02, 0.42, 0.46),
                         label_size = 16, labels = "H")

fig2_panels_G_H <- plot_grid(fig2_panel_G, fig2_panelH, rel_heights = c(0.3, 0.7), ncol = 1)
ggsave(plot = fig2_panels_G_H, "fig2_panels_G_H.png", height = 8, width = 6)

fig2_panels_A_H <- plot_grid(fig2_panels_A_F, NULL, fig2_panels_G_H, nrow = 1, rel_widths = c(0.42, 0.02, 0.56))
ggsave("fig2_panels_A_H.png", height = 8, width = 14)

###for slide
fucci_slide <- plot_grid(U2OS_fucci_AS_Sphase, NULL, U2OS_fucci_AS_G2Mphase,
                         nrow = 1,
                         rel_widths = c(0.49, 0.02, 0.49))
ggsave(plot = fucci_slide, "fig2_panels_J_K_nolab.png", height = 4, width = 6)
ggsave(plot = pn_bis, "sfig2_E_nolab.png", height = 5, width = 7)


##### Supplementary Figure 2

load("/Users/juanma/GDrive/LabData/Microscopy/20200114 FUCCI.1 analysis/20200114 c10T FUCCI.rdata")
load(file = "20200117 FUCCI individual tracks.rdata")
load(file = "20200117_RKO_FUCCI4_plot.rdata")

RKO_fucci_SG2 <- RKO_fucci_SG2 + labs(title = "RKO S-G2/M Phase")
c10T_fucci <- c10T_fucci + labs(title = "c10T1/2 S-G2/M Phase")
FUCCI_ind_tracks <- FUCCI_ind_tracks + labs(title = "U2OS Individual Cell Tracks")

WB_RKO_IDH2 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_RKO_idh2_wb.png", scale = 0.95)
WB_c10T_IDH2 <- cowplot::ggdraw(xlim = c(-0.01, 1.01), ylim = c(-0.05, 1.05)) + cowplot::draw_image("fig2_c10T_idh2_wb.png", scale = 0.95)

sfig2_panels_BCDE <- plot_grid(WB_RKO_IDH2, WB_c10T_IDH2, RKO_fucci_SG2, c10T_fucci,
                              ncol = 2,
                              #rel_widths = c(0.3, 0.05, 0.3, 0.05, 0.3),
                              labels = c("B", "C", "D", "E"),
                              label_size = 16)
sfig2_panels_A <- plot_grid(pn_bis,
                            labels = "A",
                            label_size = 16)
sfig2_panels_A_E <- plot_grid(sfig2_panels_A, sfig2_panels_BCDE,
                              ncol = 1,
                              rel_heights = c(1/2, 1/2))
ggsave("sfigure2.png", height = 10, width = 6)

### for slide
sfig2_panels_ABCD_nolab <- plot_grid(RKO_fucci_SG2, c10T_fucci, NULL, NULL, WB_RKO_IDH2, WB_c10T_IDH2,
                               ncol = 2, rel_heights = c(0.6, 0.02, 0.4))
ggsave(plot = sfig2_panels_ABCD_nolab, "sfig2_panels_ABCD_nolab.png", height = 6, width = 8)







