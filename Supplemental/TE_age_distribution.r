#!/apps/bin/Rscript
#Author: Yiran Hou

# Load required packages
library(ggplot2)
library(ggpubr)
library(ggh4x)

## All active states
TEactive <- read.table("nm.anns", header = F, col.names = c("TE_Class","JC_GeneticDistance", "State"))
my_comparisons <- list(c("LMR","Others"),c("UMR","Others"),c("DistalATAC","Others"),c("ProximalATAC","Others"),c("WeakPromoter","Others"),c("ActivePromoter","Others"),c("ActiveEnhancer","Others"))

TEactive %>% 
  mutate(State = factor(State, levels = c("ActiveEnhancer","ActivePromoter","WeakPromoter","ProximalATAC","DistalATAC","UMR","LMR","Others"))) %>%
  ggviolin(x = "State", y = "JC_GeneticDistance", fill = "State") + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y = 1) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  scale_y_continuous(
    breaks = function(x) {
      x <- scales::extended_breaks()(x)
      x[x < 1.2]
    },
    guide = guide_axis_truncated(trunc_lower = -Inf, trunc_upper = 1.2))
ggsave("JC_GeneticDistance_vln_allTE_groupby_activity_allActive.stats.pdf", width = 6, height = 6, dpi = 300)

## Focus on testis TETSS
TEdiv <- read.table("danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann", header = F, comment.char= "",col.names = c("chr","start","stop","TEFamily_TEClass","TEFamily","TEClass","TE_Class","TE","TElength","#substitutions","ProportionOfSubstitutions","JC_GeneticDistance","act_enh","cre","lmr","umr_lmr","act_pro","distal_atac","proximal_atac","weak_pro","atac","heterochromatin","umr","testis_tss"), sep = "\t")

TEdiv %>% 
  mutate(testis_tss = factor(testis_tss, levels = c("Testis_TETSS", "Others"))) %>%
  ggviolin(x = "testis_tss", y = "JC_GeneticDistance", fill = "testis_tss") +
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA) + 
  # stat_compare_means(label = "p.signif") +
  stat_compare_means(label.y = 1) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(axis.title.x = element_blank())
ggsave("JC_GeneticDistance_vln_innerBox_allTE_groupby_activity.testis_tss.stats.pdf", width = 3.5, height = 5, dpi = 300)