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

load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201015 RKO + U2OS EdU gH2AX RPA32pS4S8/U2OS_WTv172_FdA/20201015 U2OS_IDH2 2hg.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201125 U2OS 2hg + pTurn IDH1 and 2 + RKO repeat 2hg/20201119_U2OS_IDH1.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201125 U2OS 2hg + pTurn IDH1 and 2 + RKO repeat 2hg/20201015 RKO 2hg aKg Aph.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201015 RKO + U2OS EdU gH2AX RPA32pS4S8/RKO_FdA/20201015 RKObleo.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201119 RKO 2hg akg + pTurnIDH/20201119_RKO_IDH2.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20201119 RKO 2hg akg + pTurnIDH/20201119_RKO_IDH1.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20220120 10T RKO WTvM +- dox/20201119_c10T_IDH2.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20211026 U2OS time course + guides/20201119_U2OS_IDH2.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20211029 U2OS sgCBX,2f,EHMT/20211029 U2OS guides rep1.rdata")
load("/Users/juanma/GDrive/LabData/Flow Cytometry/20220204 10T+RKO WTvM/20220204_RKO.rdata")

viridisE4 <- viridis_pal(option = "D", begin = 0, end = 0.95, direction = -1)(4)
viridisE4f <- viridisE4[c(2,1,4,3)]

###Figure 3
p1 <- U2OS_pT_IDH2_EdU_plot2 #U2OS_pT_IDH2_plot for EdU+/-
p2 <- U2OS_pT_IDH2_2dplot_bis
p3 <- U2OS_pT_IDH2_stats_plot
p4 <- U2OS_IDH2_plot_EdU2 + scale_fill_manual(values = viridisE4f, guide = guide_legend(reverse = TRUE)) +
                            scale_color_manual(values = viridisE4f, guide = guide_legend(reverse = TRUE))  +
                            theme(legend.position = "none") #U2OS_IDH2_plot for EdU+/-

fig3_panels_abc <- plot_grid(p1, p3, p4,
                             align = "h",
                             axis = "b",
                             nrow = 1,
                             rel_widths = c(0.4, 0.25, 0.35),
                             labels = c("A", "B", "C"),
                             label_size = 16)
ggsave("fig3_panels_A_C.png", height = 4, width = 10)

# p5 <- RKO_pT_IDH2_plot
# p6 <- RKO_pT_IDH2_2dplot
# p7 <- RKO_pT_IDH2_stats_plot
# p6b <- RKO_pT_IDH2_2dplot_bis
# p8 <- RKO_plot_select + scale_fill_manual(values = viridisE4f) + theme(legend.position = "none")
# 
# fig3_panels_efgh <- plot_grid(p5, p7, p8,
#                              align = "hv",
#                              axis = "b",
#                              nrow = 1,
#                              rel_widths = c(0.4, 0.25, 0.35),
#                              labels = c("D", "E", "F"),
#                              label_size = 16)
# ggsave("fig3_panels_E_H.png", height = 4, width = 10)

WB_CBX1 <- cowplot::ggdraw() + cowplot::draw_image("/Users/juanma/GDrive/LabData/Flow Cytometry/20211029 U2OS sgCBX,2f,EHMT/CBX1_WB.png")
WB_CBX5 <- cowplot::ggdraw() + cowplot::draw_image("/Users/juanma/GDrive/LabData/Flow Cytometry/20211029 U2OS sgCBX,2f,EHMT/CBX5_WB.png")
WB_EHMT1 <- cowplot::ggdraw() + cowplot::draw_image("/Users/juanma/GDrive/LabData/Flow Cytometry/20211029 U2OS sgCBX,2f,EHMT/EHMT1_WB.png")

fig3_panels_D <- plot_grid(U2OS_IDH2_all_guides_rep1_plot_final,
                           rel_widths = 1,
                           nrow = 1,
                           labels = "D",
                           label_size = 16)

fig3_panels_E <- plot_grid(NULL, WB_CBX1, NULL, WB_CBX5, NULL, WB_EHMT1, NULL,
                            rel_widths = c(0.05, 0.27, 0.01, 0.27, 0.01, 0.27, 0.15),
                            nrow = 1,
                            labels = "E",
                            label_size = 16)

fig3_panels_D_E <- plot_grid(fig3_panels_D, fig3_panels_E,
                             ncol = 1,
                             align = "hv",
                             axis = "tblr",
                             rel_heights = c(0.65, 0.35))
ggsave(plot = fig3_panels_D_E, "fig3_panels_D_E.png", height = 6, width = 10)

###for slide
fig3_panels_D_nolabs <- plot_grid(U2OS_IDH2_all_guides_rep1_plot_final,
                           rel_widths = 1,
                           nrow = 1)
fig3_panels_E_nolabs <- plot_grid(NULL, WB_CBX1, NULL, WB_CBX5, NULL, WB_EHMT1, NULL,
                           rel_widths = c(0.05, 0.27, 0.01, 0.27, 0.01, 0.27, 0.15),
                           nrow = 1)

fig3_panels_D_E_nolabs <- plot_grid(fig3_panels_D_nolabs, fig3_panels_E_nolabs,
                             ncol = 1,
                             align = "hv",
                             axis = "tblr",
                             rel_heights = c(0.65, 0.35))
ggsave(plot = fig3_panels_D_E_nolabs, "fig3_panels_D_E_nolabs.png", height = 6, width = 10)



fig3_panels_A_E <- plot_grid(fig3_panels_abc, NULL, fig3_panels_D, NULL, fig3_panels_E,
                             ncol = 1,
                             rel_heights = c(0.35, 0.01, 0.35, 0.01, 0.15))

ggsave("fig3_panels_final.png", height = 8, width = 10)


###Supplementary Fig3
sp1 <- U2OS_pT_IDH2_plot_allreps2
sp2 <- U2OS_pT_IDH2_stats_EdU_plot_2 + scale_x_discrete(labels = c("WT", "R172K"))
sp3 <- RKO_pT_IDH2_plot_allreps2
sp4 <- RKO_pT_IDH2_stats_EdU_plot_2 + scale_x_discrete(labels = c("WT", "R172K"))

sfig3_panels_abc <- plot_grid(plotlist = list(sp1, sp2, p2),
                             align = "hv",
                             axis = "b",
                             nrow = 1,
                             rel_widths = c(0.35,0.15,0.35, 0.15),
                             labels = c("A", "B", "C"))
ggsave(plot = sfig3_panels_abc, "sfig3_panels_abc.png", height = 4, width = 10)

###for slide
sfig1_panels_A_C_nolabs <- plot_grid(plotlist = list(sp1, sp2, p2),
                                     align = "hv",
                                     axis = "b",
                                     nrow = 1,
                                     rel_widths = c(0.4, 0.15, 0.3))
ggsave(plot = sfig1_panels_A_C_nolabs, "sfig1_panels_A_C_nolabs.png", height = 4, width = 10)
ggsave(plot = p4, "fig3_panel_C.png", height = 4, width = 5)



p5 <- RKO_pT_IDH2_EdU_plot2
p6 <- RKO_pT_IDH2_2dplot
p7 <- RKO_pT_IDH2_stats_plot
p6b <- RKO_pT_IDH2_2dplot_bis
p8 <- RKO_plot_select_EdU2 + scale_fill_manual(values = viridisE4f) +
                             scale_color_manual(values = viridisE4f) +
                             theme(legend.position = "none")

sfig3_panels_def <- plot_grid(p5, p7, p8,
                             align = "hv",
                             axis = "b",
                             nrow = 1,
                             rel_widths = c(0.4, 0.25, 0.35),
                             labels = c("D", "E", "F"),
                             label_size = 16)
ggsave("fig3_panels_def.png", height = 4, width = 10)

sfig3_panels_ghi <- plot_grid(plotlist = list(sp3, sp4, p6b),
                               align = "hv",
                               axis = "b",
                               nrow = 1,
                               rel_widths = c(0.35,0.15,0.35, 0.15),
                               labels = c("G", "H", "I"))
ggsave(plot = sfig3_panels_ghi, "sfig3_panels_ghi.png", height = 4, width = 10)

### Omit IDH1
# sp5 <- RKO_pT_IDH1_EdU_plot2
# sp6 <- RKO_pT_IDH1_2dplot_bis
# sp7 <- RKO_pT_IDH1_stats_plot
# 
# 
# 
# sfig3_panels_jkl <- plot_grid(sp5, sp6, sp7,
#                              align = "hv",
#                              axis = "b",
#                              nrow = 1,
#                              rel_widths = c(0.4, 0.4, 0.2),
#                              labels = c("J", "K", "L"))
# ggsave(plot = sfig3_panels_jkl, "sfig3_panels_jkl.png", height = 4, width = 10)


sp8 <- c10T_pT_IDH2_EdU_plot2
sp9 <- c10T_pT_IDH2_2dplot_bis
sp10 <- c10T_pT_IDH2_stats_plot

sfig3_panels_jkl <- plot_grid(sp8, sp9, sp10,
                              align = "hv",
                              axis = "b",
                              nrow = 1,
                              rel_widths = c(0.4, 0.4, 0.2),
                              labels = c("J", "K", "L"))
ggsave(plot = sfig3_panels_jkl, "sfig3_panels_jkl.png", height = 4, width = 10)


sp11 <- U2OS_IDH2_2dplot
sp12 <- RKO_2dplot

sfig3_panels_mn <- plot_grid(sp11 + theme(legend.position = "none"),
                             NULL,
                             sp12 + theme(legend.position = "none"),
                             align = "vh",
                             rel_widths = c(0.49, 0.02, 0.49),
                             labels = c("M", "", "N"), nrow = 1)
ggsave("sfig3_panels_mn.png", height = 9, width = 7.5)

sp13 <- RKObleo_plot_EdU2 + theme(legend.position = "none")
sp14 <- RKObleo_2dplot + theme(legend.position = "none")
sp15 <- RKObleo_stats_plot + xlab("Condition") + scale_x_discrete(labels = c("DMSO", "Aph", "Bleo", "oct-2HG"))

sfig3_panels_op <- plot_grid(sp13, NULL, sp15,
                             align = "v",
                             axis = "l",
                             rel_heights = c(0.49, 0.02, 0.49),
                             labels = c("O", "", "P"),
                             ncol = 1)
sfig3_panels_opq <- plot_grid(sfig3_panels_op, NULL, sp14,
                              align = "",
                              axis = "",
                              rel_widths = c(0.49, 0.02, 0.49),
                              labels = c("", "", "Q"),
                              nrow = 1)

sfig3_panels_mnopq <- plot_grid(sfig3_panels_mn, NULL, sfig3_panels_opq,
                             rel_widths = c(0.49, 0.02, 0.49),
                             nrow = 1)
ggsave(plot = sfig3_panels_mnopq, "sfig3_panels_mnopq.png", height = 6, width = 12)


sfig3_panels_A_Q <- plot_grid(sfig3_panels_abc, NULL,
                              sfig3_panels_def, NULL,
                              sfig3_panels_ghi, NULL,
                              sfig3_panels_jkl, NULL,
                              sfig3_panels_mnopq,
                              rel_heights = c(0.1, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1, 0.01, 0.3),
                              ncol = 1)

ggsave("sfig3_panels_final.png", height = 22, width = 14)











