#!/apps/bin/Rscript
#Author: Yiran Hou

# Load required packages
library(ggplot2)
library(ggpubr)

# Load data and reformat
te <- read.table("all_ts_te_de_10.simp.txt", header=TRUE)
te.long <- melt(te, id.vars=c("TU","Group"))
colnames(te.long) <- c("TU","Group","Samples","value")

# Relevel and plot violin and box plot for the expression level of differentially expressed TEs across tissues. 
te.long %>% 
  mutate(Group = factor(Group, levels = c("Brain","Testis","Skin","Muscle","Heart","Kidney","Spleen","Blood","Liver","Intestine","Colon","e_Brain","e_Trunk")), Samples = factor(Samples, levels = c("Brain_rep1","Brain_rep2","Testis_rep1","Testis_rep2","Skin_rep1","Skin_rep2","Muscle_rep1","Muscle_rep2","Heart_rep1","Heart_rep2","Kidney_rep1","Kidney_rep2","Spleen_rep1","Spleen_rep2","Blood_rep1","Blood_rep2","Liver_rep1","Liver_rep2","Intestine_rep1","Intestine_rep2","Colon_rep1","Colon_rep2","e_Brain_rep1","e_Brain_rep2","e_Trunk_rep1","e_Trunk_rep2"))) %>%
  ggviolin(x = "Samples", y = "value", fill = "Samples") + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x=element_blank()) +
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA) + 
  geom_hline(yintercept = mean(te.long$value), linetype = 2) +
  facet_grid(rows=vars(Group)) + 
  stat_compare_means(method = "anova", label.y = 17, size=3)+        
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE) + ylab("vst-transformed counts")
ggsave("allTSte_vst_vln_groupby_tissue.pdf", width = 10, height = 12, dpi = 300)