rm(list = ls())
library(plotly)
library(dplyr)
library(ggpubr)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
trait1

# manhattan

p <- LD.plot(data_5)
p + xlim(0,30)

manhattan.plot(data = data_5, traits=c ("FA1_all", "MET_all"))
# qq.plot(data = data_5, trait = "MET_jul")


P1 <- manhattan.plot(data = data_5, traits=c ("may_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
P2 <- manhattan.plot(data = data_5, traits=c ("jun_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
P3 <- manhattan.plot(data = data_5, traits=c ("jul_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
P4 <- manhattan.plot(data = data_5, traits=c ("aug_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P5 <- manhattan.plot(data = data_5, traits=c ("sep_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P6 <- manhattan.plot(data = data_5, traits=c ("may_21_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P7 <- manhattan.plot(data = data_5, traits=c ("jun_21_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P8 <- manhattan.plot(data = data_5, traits=c ("jul_21_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P9 <- manhattan.plot(data = data_5, traits=c ("aug_21_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P10 <- manhattan.plot(data = data_5, traits=c ("sep_21_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P11 <- manhattan.plot(data = data_5, traits=c ("MET_may")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P12 <- manhattan.plot(data = data_5, traits=c ("MET_jun")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P13 <- manhattan.plot(data = data_5, traits=c ("MET_jul")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P14 <- manhattan.plot(data = data_5, traits=c ("MET_aug")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
P15 <- manhattan.plot(data = data_5, traits=c ("MET_sep")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())


Q1 <- manhattan.plot(data = data_5, traits=c ("MET_20")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q2 <- manhattan.plot(data = data_5, traits=c ("MET_21")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q3 <- manhattan.plot(data = data_5, traits=c ("FA1_all")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q4 <- manhattan.plot(data = data_5, traits=c ("PH_19_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q5 <- manhattan.plot(data = data_5, traits=c ("PH_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q6 <- manhattan.plot(data = data_5, traits=c ("FA1_PH")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
Q7 <- manhattan.plot(data = data_5, traits=c ("CT_20_1stage")) + theme_classic(base_family = "Arial", base_size = 12) + scale_color_manual(values=c("aquamarine4","azure4")) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())


myplot1 <- ggarrange(P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15,
                     labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"), ncol = 5, nrow = 3)
myplot2 <- ggarrange(Q1, Q2, Q3, Q4, Q5, Q6, Q7, 
                     labels = c("a", "b", "c", "d", "e", "f", "g"), ncol = 3, nrow = 3)

ggsave(filename = "myplot1.jpg", plot = myplot1, width = 15, height = 9)
ggsave(filename = "myplot2.jpg", plot = myplot2, width = 9, height = 9)

