setwd("D:/Sunflower_stability_GWAS/Stability/")

library(ggplot2)
library(ggpubr)


oil_stabilities <- read.csv("oil_stabilities.csv")

FAD_status <- read.csv("../Phenotype_Analysis/FAD2-1_mutant_status_renamed.csv")
colnames(FAD_status) <- c("Line", "FAD_Flag")
FAD_status$FAD_Flag <- as.factor(ifelse(FAD_status$FAD_Flag==0, "FAD WT", "FAD mutant"))
oil_stabilities <- na.omit(merge(oil_stabilities, FAD_status, by.x = "X", by.y = "Line"))

p1 <- ggplot(oil_stabilities, aes(y = oleic_stability, x= as.factor(FAD_Flag)#, fill = as.factor(FAD_Flag)
                                  ))+
  geom_boxplot() + 
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(title="Oleic Acid",y="Beta Stability", x = NULL)

p2 <- ggplot(oil_stabilities, aes(y = linoleic_stability, x= as.factor(FAD_Flag)#, fill = as.factor(FAD_Flag)
))+
  geom_boxplot() + 
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(title="Linoleic Acid",y="Beta Stability", x = NULL)
p2

p3 <- ggplot(oil_stabilities, aes(y = stearic_stability, x= as.factor(FAD_Flag)#, fill = as.factor(FAD_Flag)
))+
  geom_boxplot() + 
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(title="Stearic Acid",y="Beta Stability", x = NULL)
p3

p4 <- ggplot(oil_stabilities, aes(y = palmitic_stability, x= as.factor(FAD_Flag)#, fill = as.factor(FAD_Flag)
))+
  geom_boxplot() + 
  stat_compare_means(method = "t.test")+
  theme_bw()+
  labs(title="Palmitic Acid",y="Beta Stability", x = NULL)
p4

tiff("Figure_4.tiff", width = 500, height = 500)
ggarrange(p1,p2,p3,p4)
dev.off()