######################################################
### Krys Kibler 2025-05-27
### Making the figures
######################################################

library(tidyverse)
library(ggpubr)
library(ggtext)
library(patchwork)
library(zoo)


######################################################
### Figure 2
### Differential ani on mapping and nucl diversity

# sum mapped reads
dataframe <- instrain_genome %>% select(SampleDate,count,sum_reads_92,sum_reads_95) %>% 
  pivot_longer(-c(1:2), names_to = "ANI_level", values_to = "sum_reads") %>% distinct()

dataframe %>% filter(ANI_level == "sum_reads_92") %>% summarise(sum=sum(count)) # 371 samples
dataframe %>% na.omit() %>% filter(ANI_level == "sum_reads_95") %>% summarise(sum=sum(count)) # 242 samples
dataframe$ANI_level <- gsub("sum_reads_", "", dataframe$ANI_level)
dataframe$ANI_level <- gsub("2", "2%", dataframe$ANI_level)
dataframe$ANI_level <- gsub("5", "5%", dataframe$ANI_level)

p1 <- dataframe %>% ggplot(aes(x = ANI_level, y = sum_reads, fill = ANI_level)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox", label.y = 1.5E7, size = 6) +
  geom_text(aes(x = "92%", y = 1.25E7, label = "n = 371", size = 40)) +
  geom_text(aes(x = "95%", y = 1.25E7, label = "n = 242", size = 40)) +
  coord_cartesian(ylim = c(0,1.75E7)) +
  scale_fill_manual(values = c("92%" = "#F8766D",
                               "95%" = "#00BFC4")) +
  labs(x = "Minimum read ANI percentage", y = "Total mapped reads", title = "Total mapped reads between 92% and 95% read ANI") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "none")

# sum normalized coverage between the ani level
dataframe <- instrain_genome %>% select(SampleDate,count,sum_cov.92,sum_cov.95) %>% 
  pivot_longer(-c(1:2), names_to = "ANI_level", values_to = "cov") %>% distinct()

dataframe$ANI_level <- gsub("sum_cov.", "", dataframe$ANI_level)
dataframe$ANI_level <- gsub("2", "2%", dataframe$ANI_level)
dataframe$ANI_level <- gsub("5", "5%", dataframe$ANI_level)

dataframe %>% filter(ANI_level == "92%") %>% summarise(sum=sum(count)) # 371 samples
dataframe %>% na.omit() %>% filter(ANI_level == "95%") %>% summarise(sum=sum(count)) # 242 samples


p2 <- dataframe %>% ggplot(aes(x = ANI_level, y = cov, fill = ANI_level)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox", label.y = 825, size = 6) +
  geom_text(aes(x = "92%", y = 725, label = "n = 371", size = 40)) +
  geom_text(aes(x = "95%", y = 725, label = "n = 242", size = 40)) +
  coord_cartesian(ylim = c(0,950)) +
  labs(x = "Minimum read ANI percentage", y = "Total normalized coverage", title = "Cumulative coverage between 92% and 95% read ANI") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "none")

# coverage variability across the genome between ani levels
dataframe <- instrain_genome %>% select(SampleDate,tax_abbr_name,tax_abbr_name.col,coverage_SEM_92ANI,coverage_SEM_95ANI) %>% na.omit() %>% 
  pivot_longer(-c(1:3), names_to = "ANI_level", values_to = "coverage_SEM")

dataframe$ANI_level <- gsub("coverage_SEM_", "", dataframe$ANI_level)
dataframe$ANI_level <- gsub("2ANI", "2%", dataframe$ANI_level)
dataframe$ANI_level <- gsub("5ANI", "5%", dataframe$ANI_level)

p3 <- dataframe %>% ggplot(aes(x = factor(tax_abbr_name.col, levels = levels), y = coverage_SEM, fill = ANI_level)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(label = after_stat(p.signif)),paired = TRUE, label.y = .06, size = 10, method = "wilcox.test") +
  labs(x = "Cyano-mOTUs", y = "Coverage SEM", title = "Coverage SEM between 92% and 95% read ANI") +
  coord_cartesian(ylim = c(0,0.075), expand = FALSE) +
  #facet_wrap(~group, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle=45, hjust=1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "none")

# coverage variability but with breadth
dataframe <- instrain_genome %>% select(SampleDate,tax_abbr_name,tax_abbr_name.col,breadth_92ANI,breadth_95ANI,breadth_expected_92ANI,breadth_expected_95ANI)

dataframe$difference_92 <- dataframe %>% with(breadth_expected_92ANI - breadth_92ANI)
dataframe$difference_95 <- dataframe %>% with(breadth_expected_95ANI - breadth_95ANI)

dataframe <- dataframe %>% select(SampleDate,tax_abbr_name,tax_abbr_name.col,difference_92,difference_95) %>% 
  pivot_longer(-c(1:3), names_to = "ANI_level", values_to = "difference")

dataframe$ANI_level <- gsub("difference_", "", dataframe$ANI_level)
dataframe$ANI_level <- gsub("2", "2%", dataframe$ANI_level)
dataframe$ANI_level <- gsub("5", "5%", dataframe$ANI_level)

pvalue_labels <- c(.55, .60, .55, .60, .55,
                   .60, 55, .60, .55, .60,
                   .55, .60, .55, .60, .55,
                   .60, 55, .60, .55, .60,
                   .55, .60, .55, .60, .55,
                   .60, 55, .60, .55, .60,
                   .55, .60, .55, .60, .55,
                   .60, 55, .60, .55, .60,
                   .55, .60, .55, .60, .55,
                   .60, 55, .60, .55, .60)

p3.1 <- dataframe %>% na.omit() %>% ggplot(aes(x = factor(tax_abbr_name.col, levels = levels), y = difference, fill = ANI_level)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(label = after_stat(p.signif)), paired = TRUE, label.y = pvalue_labels, size = 7) +
  labs(x = "Cyano-mOTUs", y = "Coverage Uneveness", title = "Difference in expected and calculated breadth between 92% and 95% read ANI") +
  coord_cartesian(ylim = c(0,0.7), expand = FALSE) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle=45, hjust=1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "none")

# nucleotide variability between ani levels
dataframe <- instrain_genome %>% select(SampleDate,tax_abbr_name,tax_abbr_name.col,nucl_diversity_92ANI,nucl_diversity_95ANI) %>% na.omit() %>% 
  pivot_longer(-c(1:3), names_to = "ANI_level", values_to = "nucl_diversity")

dataframe$ANI_level <- gsub("nucl_diversity_", "", dataframe$ANI_level)
dataframe$ANI_level <- gsub("2ANI", "2%", dataframe$ANI_level)
dataframe$ANI_level <- gsub("5ANI", "5%", dataframe$ANI_level)

pvalue_labels <- c(.035, 0.039, .035, 0.039, .035,
                   0.039, .035, 0.039, .035, 0.039,
                   .035, 0.039, .035, 0.039, .035,
                   0.039, .035, 0.039, .035, 0.039,
                   .035, 0.039, .035, 0.039, .035,
                   0.039, .035, 0.039, .035, 0.039,
                   .035, 0.039, .035, 0.039, .035,
                   0.039, .035, 0.039, .035, 0.039,
                   .035, 0.039, .035, 0.039, .035,
                   0.039, .035, 0.039, .035, 0.039)

p4 <- dataframe %>% ggplot(aes(x = factor(tax_abbr_name.col, levels = levels), y = nucl_diversity, fill = ANI_level)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(label = after_stat(p.signif)),paired = TRUE, label.y = pvalue_labels, size = 7, method = "wilcox.test") +
  labs(x = "Cyano-mOTUs", y = "Nucleotide diversity", title = "Nucleotide diversity between 92% and 95% read ANI") +
  #coord_cartesian(ylim = c(0,0.075), expand = FALSE) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle=45, hjust=1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "none")

  
layout <- "
AB
CC
DD
"

p1+p2+p3.1+p4+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

rm(dataframe,p1,p2,p3,p3.1,p4,pvalue_labels,layout)
######################################################
### Figure 3
### Normalized coverage by motu

dataframe <- instrain_genome %>% select(tax_abbr_name,tax_abbr_name.col,group,count,SampleDate,norm_coverage_92,RA_by.cov_92)

count <- dataframe %>% group_by(tax_abbr_name) %>% summarize(sum=sum(count,na.rm =TRUE))
dataframe <- left_join(dataframe, count, by = "tax_abbr_name")

dataframe %>% group_by(tax_abbr_name) %>% summarize(median = median(norm_coverage_92))

# boxplot of coverage
p1 <- dataframe %>% ggplot(aes(fct_rev(factor(tax_abbr_name.col, levels = levels)), y = norm_coverage_92)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip(ylim = c(0,265)) +
  labs(x = "Cyano-mOTUs", y = "Normalized Coverage", title = "Normalized Coverage at 92% Read ANI") +
  theme_bw() +
  theme(axis.text.y = ggtext::element_markdown(size = 14),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=16, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin = margin(t = 10, r = 0, b = 0, l = 10, "pt"),
        title = element_text(size = 16),
        legend.position = "none") +
  annotate("text",label = paste("n =",dataframe$sum), y = 250, x = fct_rev(factor(dataframe$tax_abbr_name.col, levels = levels)), size = 5)

# sum RA by coverage
p2 <- dataframe %>% group_by(tax_abbr_name.col) %>% summarise(Sum_RA=sum(RA_by.cov_92,na.rm=TRUE)) %>% 
  ggplot(aes(fct_rev(factor(tax_abbr_name.col, levels = levels)), y = Sum_RA)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Sum RA") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=20, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        axis.ticks.length.x = unit(0, "pt"),
        title = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "pt"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border=element_blank(), axis.line.x=element_blank())

# averaged coverage by group, time series plot shaded lines
dataframe$jday <- yday(dataframe$SampleDate) 
dataframe_backup <- dataframe


dataframe <- dataframe_backup
dataframe <- dataframe %>% group_by(SampleDate,jday,group) %>% summarise(sum_cov=sum(norm_coverage_92,na.rm = TRUE))

dataframe <- dataframe %>% group_by(jday, group) %>% summarise(avg_cov=mean(sum_cov, na.rm = TRUE))

dataframe <- dataframe %>% group_by(group) %>% mutate(value = rollmean(avg_cov, k = 6, fill = FALSE))

dataframe$date <- as.Date(dataframe$jday, origin=as.Date("2000-01-01"))


p3 <- dataframe %>% ggplot(aes(x = date, y = value, group = group, fill = group)) +
  geom_area() +
  coord_cartesian(xlim = as.Date(c("2000-05-13","2000-11-14")), expand = FALSE) +
  scale_fill_manual(values = c("Cyanobium" = "#00dede",
                               "Filamentous-1" = "#236a00",
                               "Microcystis" = "#e39f0c",
                               "Nodosilinea" = "#3cb306",
                               "Pseudanabaena" = "#4cd505",
                               "Snowella" = "#f4659d",
                               "Solitary-1" = "#008bff",
                               "Vulcanococcus" = "#7a05c6")) +
  labs(x = "Month", y = "Cumulative normalized coverage", title = "Rolling weekly average of combined normalized\ncoverage by cyano-mOTU group", fill = "Cyano-mOTU groups") +
  theme_bw() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=20, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=20, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "bottom",
        legend.box="horizontal",
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size = 10)))

layout <- "
AAB
"

patch <- wrap_elements((p1+p2) + plot_layout(design = layout))

layout <- "
A#
AB
AB
AB
A#
" 

patch+free(p3)+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

rm(layout,pvalue_labels,p1,p2,p3,p3.1,patch,count)

######################################################
### Nucleotide diversity

dataframe <- instrain_genome %>% select(tax_abbr_name.col, tax_abbr_name, group, morphology, SampleDate, breadth_92ANI, coverage_92ANI, nucl_diversity_92ANI, RA_by.cov_92, norm_coverage_92)
dataframe <- dataframe %>% filter(breadth_92ANI >= 0.8) %>% filter(coverage_92ANI >= 10)
# pi by morphology 

# Statistics
# https://statsandr.com/blog/anova-in-r/ #link for normal anova
# https://statsandr.com/blog/kruskal-wallis-test-nonparametric-version-anova/ for non-normal anova

# anova with normal distribution assumption
res_aov <- aov(nucl_diversity_92ANI ~ morphology,
               data = dataframe)

par(mfrow = c(1, 2)) # combine plots
# check distribution - #NOT normal
# the look
# histogram
hist(res_aov$residuals)
# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE) # id = FALSE to remove point identification
# the numbers
shapiro.test(res_aov$residuals) # NOT normal

# test for variance with Levene # NOT even variance
leveneTest(nucl_diversity_92ANI ~ morphology,
           data = dataframe) 

# Anova with non-normal distribution assumption
kruskal.test(nucl_diversity_92ANI ~ morphology,
             data = dataframe)
#	Kruskal-Wallis rank sum test
#data:  nucl_diversity_92ANI by morphology
#Kruskal-Wallis chi-squared = 902.06, df = 3, p-value < 2.2e-16

# post-hoc Dunn Test
library(FSA)
stats <- dunnTest(nucl_diversity_92ANI ~ morphology,
                  data = dataframe,
                  method = "holm")
stats <- stats$res
stats <- stats %>% mutate(across(Comparison, ~factor(., levels = c("Colonial - Filamentous", "Filamentous - Solitary", 
                                                                   "Colonial - Stalked", "Solitary - Stalked",
                                                                   "Colonial - Solitary", "Filamentous - Stalked")))) %>% arrange(Comparison)
stats <- stats %>% separate(Comparison, c("group1", "group2"), sep = " - ")
stars <- c("****", "****", "**", "****", "****", "****")
stats <- cbind(stats, stars)

#Dunn (1964) Kruskal-Wallis multiple comparison
#p-values adjusted with the Holm method.

#group1      group2          Z       P.unadj         P.adj
#1    Colonial Filamentous  29.640842 4.450849e-193 2.670510e-192 0.00000
#2 Filamentous    Solitary -16.360228  3.677493e-60  1.470997e-59 0.00000
#3    Colonial     Stalked   3.270963  1.071818e-03  1.071818e-03 0.00107
#4    Solitary     Stalked  -4.715564  2.410419e-06  4.820838e-06 0.00000
#5    Colonial    Solitary  19.264995  1.056750e-82  5.283751e-82 0.00000
#6 Filamentous     Stalked -10.582358  3.597867e-26  1.079360e-25 0.00000

library(ggsignif)

dataframe <- dataframe %>% mutate(across(morphology, ~factor(., levels=c("Filamentous","Colonial","Stalked","Solitary")))) 

dataframe %>% group_by(morphology) %>% 
  summarize(mean = mean(nucl_diversity_92ANI), sd = sd(nucl_diversity_92ANI))

#morphology     mean      sd
#<fct>         <dbl>   <dbl>
#1 Filamentous 0.00501 0.00274
#2 Colonial    0.0175  0.00335
#3 Solitary    0.0102  0.00685
#4 Branched    0.0138  0.00300

dataframe %>% group_by(group) %>% 
  summarize(mean = mean(nucl_diversity_92ANI), sd = sd(nucl_diversity_92ANI))

#group            mean      sd
#<chr>           <dbl>   <dbl>
#1 Cyanobium     0.00995 0.00598
#2 Filamentous-1 0.00472 0.00220
#3 Microcystis   0.0175  0.00335
#4 Nodosilinea   0.00657 0.00235
#5 Pseudanabaena 0.00448 0.00398
#6 Snowella      0.0138  0.00300
#7 Solitary-1    0.00680 0.00383
#8 Vulcanococcus 0.0143  0.00822

# boxplot between morphologies
p1 <- dataframe %>% 
  ggplot(aes(x = morphology, y = nucl_diversity_92ANI, color = morphology)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.5) +
  geom_violin(width=1, alpha = 0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.8, outlier.shape = NA) +
  stat_pvalue_manual(stats, label = "stars", y.position = .045, 
                     step.increase = 0.045, tip.length = 0.01, size = 6) +
  scale_color_manual(values = c("Filamentous" = "#0ab053",
                                "Colonial" = "#e39f0c",
                                "Solitary" = "#7660f6",
                                "Stalked" = "#f4659d")) +
  labs(y = "Nucleotide Diversity (π)", x = "Cyano-mOTUs Morphology", title = "     Nucleotide diversity between morphologies") +
  theme_bw() +
  theme(title = element_text(size = 22),
        axis.text.x = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=26, face="bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=26, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        plot.margin = margin(t = 10, r = 0, b = 10, l = 5),
        aspect.ratio = 3/2,
        legend.position = "none")

# boxplot of filamentous groups
p2 <- dataframe %>% filter(morphology == "Filamentous") %>%  
  mutate(across(group, ~factor(., levels=c("Pseudanabaena","Filamentous-1","Nodosilinea")))) %>% 
  ggplot(aes(x = group, y = nucl_diversity_92ANI, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(comparisons = list(c("Pseudanabaena","Filamentous-1"),
                                 c("Filamentous-1","Nodosilinea"),
                                 c("Pseudanabaena","Nodosilinea")),
              map_signif_level = TRUE, textsize = 6, test = "wilcox.test",
              margin_top = 0.01,step_increase = 0.05, tip_length = 0.01, vjust = 0.5, y_position = 0.015) +  
  coord_cartesian(ylim = c(0,0.018), expand = TRUE) +
  scale_fill_manual(values = c("Pseudanabaena" = "#4cd505",
                               "Filamentous-1" = "#236a00",
                               "Nodosilinea" = "#3cb306")) +
  labs(x = "Filamentous clades", y = "π") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        aspect.ratio=1,
        legend.position = "none")

p3 <- dataframe %>% filter(morphology == "Filamentous") %>% filter(group == "Filamentous-1") %>% 
  ggplot(aes(x = fct_reorder(tax_abbr_name,nucl_diversity_92ANI,median), y = nucl_diversity_92ANI)) +
  geom_boxplot(outlier.shape = NA, color = "#236a00") +
  coord_cartesian(ylim = c(0,0.018), expand = TRUE) +
  labs(x = "Filamentous-1 mOTUs", y = "π") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 20, r = 10, b = 0, l = 10)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),,
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        aspect.ratio = 1,
        legend.position = "none")

layout <- "
AB
"
patch <- wrap_elements((p2+p3) + plot_layout(design = layout)) 

# Colonial mOTUs
p4 <- dataframe %>% filter(morphology == "Colonial") %>% 
  ggplot(aes(x = tax_abbr_name, y = nucl_diversity_92ANI)) +
  geom_boxplot(outlier.shape = NA, color = "#e39f0c") +
  coord_cartesian(ylim = c(0.005,0.028), expand = TRUE) +
  geom_signif(comparisons = list(c("MCYST_2","MCYST_31"),
                                 c("MCYST_31","MCYST_62"),
                                 c("MCYST_2","MCYST_62")),
              map_signif_level = TRUE, textsize = 6, test = "wilcox.test",
              margin_top = 0.01,step_increase = 0.05, tip_length = 0.01, vjust = 0.5, y_position = 0.025) +
  labs(x = "Microcystis mOTUs", y = "π") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin = margin(t = 10, r = 0, b = 10, l = 10, "pt"),
        aspect.ratio=1,
        legend.position = "none")

# Solitary mOTUs

p5 <- dataframe %>% filter(morphology == "Solitary") %>%  
  mutate(across(group, ~factor(., levels=c("Cyanobium","Solitary-1","Vulcanococcus")))) %>% 
  ggplot(aes(x = group, y = nucl_diversity_92ANI, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(comparisons = list(c("Cyanobium","Solitary-1"),
                                 c("Solitary-1","Vulcanococcus"),
                                 c("Cyanobium","Vulcanococcus")),
              map_signif_level = TRUE, textsize = 6, test = "wilcox.test",
              margin_top = 0.01,step_increase = 0.05, tip_length = 0.01, vjust = 0.5, y_position = 0.03) +  
  coord_cartesian(ylim = c(0,0.036), expand = TRUE) +
  scale_fill_manual(values = c("Cyanobium" = "#00dede",
                               "Solitary-1" = "#008bff",
                               "Vulcanococcus" = "#7a05c6")) +
  labs(x = "Solitary clades", y = "π") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        aspect.ratio = 1,
        legend.position = "none")

# How does average abundance rank pair with nucleotide diversity
dataframe_backup <- dataframe

dataframe <- dataframe_backup
dataframe <- dataframe %>% group_by(tax_abbr_name, tax_abbr_name.col, group) %>% 
  summarise(avg_nucl=mean(nucl_diversity_92ANI), avg_cov=mean(norm_coverage_92))

p6 <- dataframe %>% mutate(across(group, ~factor(., levels=c("Pseudanabaena","Nodosilinea","Filamentous-1","Snowella","Microcystis","Cyanobium","Solitary-1","Vulcanococcus")))) %>% 
  ggplot(aes(x = avg_cov, y = avg_nucl, color = group)) +
  geom_point(aes(size = 2)) +
  scale_colour_manual(values = c("Pseudanabaena" = "#4cd505",
                                 "Nodosilinea" = "#3cb306",
                                 "Filamentous-1" = "#236a00",
                                 "Snowella" = "#f4659d",
                                 "Microcystis" = "#e39f0c",
                                 "Cyanobium" = "#00dede",
                                 "Solitary-1" = "#008bff",
                                 "Vulcanococcus" = "#7a05c6")) +
  labs(x = "Average norm coverage", y = "Average π", color = "Cyano-mOTU Clades") +
  theme_bw() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=20, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=20, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        title = element_text(size=16),
        legend.position = "right",
        legend.box="vertical",
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size = 7))) +
  scale_size(guide = 'none')
 

# combination station

layout <- "
ABC
ADE
AFF
"
p1+p2+p3+p4+p5+free(p6)+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

rm(p1,p2,p3,p4,p5,p6,layout,patch,res_aov,stats,stars)

######################################################
### Interannual variability

dataframe <- instrain_genome %>% select(tax_abbr_name.col, tax_abbr_name, group, morphology, SampleDate, breadth_92ANI, coverage_92ANI, nucl_diversity_92ANI)
dataframe <- dataframe %>% filter(breadth_92ANI >= 0.8) %>% filter(coverage_92ANI >= 10)
dataframe_backup <- dataframe

dataframe <- dataframe_backup

# clear out mOTUs with less than 20 samples
dataframe$count <- 1
count <- dataframe %>% group_by(tax_abbr_name) %>% summarise(count_samples=sum(count))
dataframe <- left_join(dataframe, count, by = "tax_abbr_name")
dataframe <- dataframe %>% filter(count_samples >= 20) # 30 remain out of 50

dataframe %>% ggplot(aes(x = SampleDate, y = nucl_diversity_92ANI, color = group)) +
  geom_line() +
  facet_wrap(~factor(tax_abbr_name.col,levels=levels), scales = "free_y") +
  scale_color_manual(values = c("Pseudanabaena" = "#4cd505",
                                "Nodosilinea" = "#3cb306",
                                "Filamentous-1" = "#236a00",
                                "Snowella" = "#f4659d",
                                "Microcystis" = "#e39f0c",
                                "Cyanobium" = "#00dede",
                                "Solitary-1" = "#008bff",
                                "Vulcanococcus" = "#7a05c6")) +
  scale_x_date(date_labels = "'%y", breaks = "2 year") +
  labs(x = "Year", y = "π") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white"))



dataframe$year4 <- year(dataframe$SampleDate)
dataframe <- dataframe %>% group_by(tax_abbr_name, tax_abbr_name.col, group, morphology, year4) %>% summarise(mean_nucl_92=mean(nucl_diversity_92ANI))
dataframe$date <- as_date(paste(dataframe$year4,"-01-01",sep = ""))

# time series of nuc diversity
p1 <- dataframe %>% 
  ggplot(aes(x = date, y = mean_nucl_92, group = tax_abbr_name, color = group)) +
  geom_line() +
  facet_wrap(~factor(group,levels=c("Pseudanabaena","Nodosilinea","Filamentous-1","Microcystis","Snowella","Vulcanococcus","Cyanobium","Solitary-1")), ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Pseudanabaena" = "#4cd505",
                                "Nodosilinea" = "#3cb306",
                                "Filamentous-1" = "#236a00",
                                "Snowella" = "#f4659d",
                                "Microcystis" = "#e39f0c",
                                "Cyanobium" = "#00dede",
                                "Solitary-1" = "#008bff",
                                "Vulcanococcus" = "#7a05c6")) +
  scale_x_date(date_labels = "'%y", breaks = "year") +
  labs(x = "Year", y = "Annual mean π") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white"))

# trend line
library(broom)

linear_fit <- dataframe %>%
  group_by(tax_abbr_name) %>% 
  do(tidy(lm(mean_nucl_92~year4, data = .))) %>% 
  select(variable = tax_abbr_name, term, t_stat = statistic, p.value, slope = estimate)


fitted_models <- linear_fit
fitted_models <- mutate(fitted_models, pass = ifelse(p.value < 0.05, 'pass', 'fail'))
fitted_models <- mutate(fitted_models, trend = ifelse(slope > 0, 'pos', 'neg'))

fitted_models <- left_join(fitted_models, cyano.mOTUs, by = c("variable" = "tax_abbr_name"))
fitted_models$category <- paste(fitted_models$pass, fitted_models$trend, sep = "-")

fitted_models$category <- gsub("fail-neg", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("fail-pos", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("NA.+", "No-sig change", fitted_models$category)

fitted_models$category <- gsub("pass-neg", "Negative, Sig", fitted_models$category)
fitted_models$category <- gsub("pass-pos", "Positive, Sig", fitted_models$category)

fitted_models$count <- 1

fitted_models <- fitted_models %>% filter(term == "year4")  

fitted_models %>% 
  group_by(category) %>% 
  summarise(sum = sum(count))
#pass-pos = 0
#pass-neg = 2
#no-trend = 28


library(gt)
fitted_models$p.value <- signif(fitted_models$p.value,2)
fitted_models$slope <- signif(fitted_models$slope,2)

save_table <- fitted_models %>% ungroup() %>% select(variable, t_stat, p.value, slope, pass, trend, category, group)
write_csv(save_table, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/interannual_trend.csv")

fitted_models <- fitted_models %>% ungroup() %>% select(variable,p.value,slope,category,group)

colnames(fitted_models) <- c("Cyano-mOTU","p.value","Slope","Trend","Group")

table <- fitted_models %>% select(c("Cyano-mOTU","Group","p.value","Slope","Trend")) %>% 
  arrange(factor(Group, levels=c("Pseudanabaena", "Nodosilinea", "Filamentous-1", "Microcystis", "Snowella", "Vulcanococcus", "Cyanobium", "Solitary-1")))

table %>% gt()
table <- wrap_table(table, panel = "full", space = "fixed")

free(p1) + table

rm(p1, table, linear_fit, fitted_models, count)

######################################################
### Nucleotide diversity and iRep
dataframe <- instrain_genome %>% select(tax_abbr_name,tax_abbr_name.col,coverage_92ANI,breadth_92ANI,nucl_diversity_92ANI,iRep_92ANI)
dataframe <- dataframe %>% filter(breadth_92ANI >= 0.8) %>% filter(coverage_92ANI >= 10)

dataframe$count <- 1
count <- dataframe %>% group_by(tax_abbr_name) %>% summarise(count_samples=sum(count))
dataframe <- left_join(dataframe, count, by = "tax_abbr_name")
dataframe <- dataframe %>% filter(count_samples >= 20) # 30 remain out of 50

dataframe_backup <- dataframe


# make sure theres enough data points in figures
dataframe <- dataframe %>% na.omit(iRep_92ANI)
count_irep <- dataframe %>% group_by(tax_abbr_name) %>% 
  summarise(count_irep = sum(count))
dataframe <- left_join(dataframe, count_irep, by = "tax_abbr_name")
dataframe <- dataframe %>% filter(count_irep >= 10) # 26 out of 30 out of 50

filt_cyanos <- unique(dataframe$tax_abbr_name)

# Set the threshold value
threshold <- 100

# Create a new variable that truncates the size values
dataframe$truncated_size <- ifelse(dataframe$coverage_92ANI >= threshold, threshold, dataframe$coverage_92ANI)

# Cyano-mOTUs - APHAN and Mcyst_2
p5 <- dataframe %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(tax_abbr_name != "MCYST_2" &
           tax_abbr_name != "APHAN_134") %>% 
  ggplot(aes(x = log10(iRep_92ANI), y = nucl_diversity_92ANI, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 6) +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        #legend.box.margin = margin(l = 0),
        plot.margin = margin(t = 10, r = 0, b = 10, l = -100, "pt"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"),
        aspect.ratio = 1/1.5) +
  guides(alpha = guide_legend(title.position="top", title.hjust = 0.5,
                              override.aes = list(size = 7)))

# Mcyst_31
p3 <- dataframe %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(tax_abbr_name == "MCYST_2") %>% 
  ggplot(aes(x = log10(iRep_92ANI), y = nucl_diversity_92ANI, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 5) +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        strip.text = ggtext::element_markdown(size = 18),
        strip.background = element_rect(fill="white"),
        aspect.ratio = 1)

# APHAN_134
p4 <- dataframe %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(tax_abbr_name == "APHAN_134") %>% 
  ggplot(aes(x = log10(iRep_92ANI), y = nucl_diversity_92ANI, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 5) +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        strip.text = ggtext::element_markdown(size = 18),
        strip.background = element_rect(fill="white"),
        aspect.ratio = 1)


# Hypothetical figures
seq <- seq(-1, 1, by=0.01)
seq <- as.data.frame(seq)



calculate_exponential <- function(x) {
  # Calculate the value of the equation y = -1 *((0.1 - x) / x )
  y <- -1 *((0.1 - x) / x )
  # Return the calculated value y
  return(y)
}
seq$calc_exp <- lapply(seq$seq, calculate_exponential)
seq <- as.data.frame(seq)

seq$calc_exp <- as.numeric(as.character(unlist(seq[[2]])))


p1 <- seq %>% filter(seq > 0 ) %>% 
  filter(calc_exp >= -5 &
           seq <= 0.5) %>% 
  ggplot(aes(x = seq, y = calc_exp)) +
  geom_line(linewidth = 2) +
  labs(x = "Replication Index (iRep)", y = "Nucleotide Diversity (π)", title = "Diversity Maintenance") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 45)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        plot.title = element_text(size = 15, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1)

p2 <- seq %>% filter(seq < 0 ) %>% 
  filter(calc_exp <= 5 &
           seq >= -0.5) %>% 
  ggplot(aes(x = seq, y = calc_exp)) +
  geom_line(linewidth = 2) +
  labs(x = "iRep", y = "π", title ="Clonal Expansion") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1)

# combination station

patch <- wrap_elements(p1+p2)
patch1 <- wrap_elements(p3+p4) 

layout <- "
AC
BC
"
patch + patch1 + p5 + plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

rm(p1,p2,p3,p4,p5,patch,patch1,layout,seq,threshold,count,count_irep)

######################################################
### Linkage Decay

dataframe <- instrain_92_linkage 

# filter for only cyanos that made it through iRep analysis
dataframe <- dataframe %>% filter(tax_abbr_name %in% filt_cyanos)
dataframe$count <- 1
count <- dataframe %>% select(tax_abbr_name,SampleDate,count) %>% distinct() %>% group_by(tax_abbr_name) %>% summarize(count_samples = sum(count))
dataframe <- dataframe %>% na.omit() %>% group_by(tax_abbr_name,tax_abbr_name.col,distance) %>% 
                                                                  summarize(r2_avg = mean(r2),
                                                                            r2_norm_avg = mean(r2_normalized),
                                                                            d_prime_avg = mean(d_prime),
                                                                            d_prime_norm_avg = mean(d_prime_normalized),
                                                                            count = sum(count))

dataframe$distance <- dataframe %>% with(distance + 1)
dataframe <- dataframe %>% ungroup() %>% group_by(tax_abbr_name, tax_abbr_name.col) %>% mutate(r2 = rollmean(r2_avg, k = 5, fill = FALSE),
                                                                         r2_norm = rollmean(r2_norm_avg, k = 5, fill = FALSE),
                                                                         d_prime = rollmean(d_prime_avg, k = 5, fill = FALSE),
                                                                         d_prime_norm = rollmean(d_prime_norm_avg, k = 5, fill = FALSE))

dataframe <- left_join(dataframe, count, by = "tax_abbr_name")
dataframe <- dataframe %>% select(-c("r2_avg", "r2_norm_avg", "d_prime_avg", "d_prime_norm_avg"))

dist = seq(5,250, by = 5)
dataframe_backup <- dataframe

dataframe <- dataframe_backup
dataframe <- dataframe %>% filter(distance %in% dist) %>% 
  pivot_longer(-c(1:3,"count","count_samples"), names_to = "stat", values_to = "linkage") 

dataframe$stat <- gsub("d_prime", "d'", dataframe$stat)

dataframe_filt <- dataframe %>% filter(tax_abbr_name == "APHAN_134" |
                                         tax_abbr_name == "MCYST_2")

p1 <- dataframe_filt %>%  
  ggplot(aes(x = distance, y = linkage, color = stat)) +
  geom_line(linewidth = 1.25) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 250, by = 50)) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.2)) +
  coord_cartesian(ylim = c(0,1.0), xlim = c(0,260), expand = FALSE) +
  facet_wrap(~factor(tax_abbr_name.col, levels=c('<span style="color: #e39f0c">MCYST_2</span>', '<span style="color: #236a00">APHAN_134</span>')), nrow = 2) +
  labs(x = "Distance between SNPs (bp)", y = "SNP linkage", color = "") +
  geom_text(label = paste("n = ", dataframe_filt$count_samples), y = 0.1, x = 50, size = 5, color = "black") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        plot.title = element_text(size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        #aspect.ratio = 1,
        legend.text = element_text(size = 16),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dataframe <- dataframe %>% filter(linkage != 0)
dataframe %>% ungroup() %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  ggplot(aes(x = distance, y = linkage, color = stat)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 200, by = 50)) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.2)) +
  coord_cartesian(ylim = c(0,1.0), xlim = c(0,260), expand = FALSE) +
  facet_wrap(~tax_abbr_name.col, ncol = 6) +
  labs(x = "Distance between SNPs (bp)", y = "SNP linkage", color = "") +
  geom_text(label = paste("n = ", dataframe$count_samples), y = 0.1, x = 75, size = 5, color = "black") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        plot.title = element_text(size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        #aspect.ratio = 1,
        legend.text = element_text(size = 16),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white")) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# trend line with d' normalized
library(broom)

linear_fit <- dataframe %>% filter(stat == "d'_norm") %>% 
  group_by(tax_abbr_name) %>% 
  do(tidy(lm(linkage~distance, data = .))) %>% 
  select(variable = tax_abbr_name, term, t_stat = statistic, p.value, slope = estimate)


fitted_models <- linear_fit
fitted_models <- mutate(fitted_models, pass = ifelse(p.value < 0.05, 'pass', 'fail'))
fitted_models <- mutate(fitted_models, trend = ifelse(slope > 0, 'pos', 'neg'))

fitted_models <- left_join(fitted_models, cyano.mOTUs, by = c("variable" = "tax_abbr_name"))
fitted_models$category <- paste(fitted_models$pass, fitted_models$trend, sep = "-")

fitted_models$category <- gsub("fail-neg", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("fail-pos", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("NA.+", "No-sig change", fitted_models$category)

fitted_models$category <- gsub("pass-neg", "Negative, Sig", fitted_models$category)
fitted_models$category <- gsub("pass-pos", "Positive, Sig", fitted_models$category)

fitted_models <- fitted_models %>% filter(term == "distance")

fitted_models$p.value <- signif(fitted_models$p.value,2)
fitted_models$slope <- signif(fitted_models$slope,2)

stat_table <- fitted_models %>% select(variable, t_stat, p.value, slope, pass, trend, group)
write_csv(stat_table, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/dprime-recombination_trend.csv")

# linear regresssion of r2
linear_fit2 <- dataframe %>% filter(stat == "r2_norm") %>% 
  group_by(tax_abbr_name) %>% 
  do(tidy(lm(linkage~distance, data = .))) %>% 
  select(variable = tax_abbr_name, term, t_stat = statistic, p.value, slope = estimate)


fitted_models2 <- linear_fit2
fitted_models2 <- mutate(fitted_models2, pass = ifelse(p.value < 0.05, 'pass', 'fail'))
fitted_models2 <- mutate(fitted_models2, trend = ifelse(slope > 0, 'pos', 'neg'))

fitted_models2 <- left_join(fitted_models2, cyano.mOTUs, by = c("variable" = "tax_abbr_name"))
fitted_models2$category <- paste(fitted_models2$pass, fitted_models2$trend, sep = "-")

fitted_models2$category <- gsub("fail-neg", "No-sig change", fitted_models2$category)
fitted_models2$category <- gsub("fail-pos", "No-sig change", fitted_models2$category)
fitted_models2$category <- gsub("NA.+", "No-sig change", fitted_models2$category)

fitted_models2$category <- gsub("pass-neg", "Negative, Sig", fitted_models2$category)
fitted_models2$category <- gsub("pass-pos", "Positive, Sig", fitted_models2$category)

fitted_models2 <- fitted_models2 %>% filter(term == "distance")

fitted_models2$p.value <- signif(fitted_models2$p.value,2)
fitted_models2$slope <- signif(fitted_models2$slope,2)

stat_table2 <- fitted_models %>% select(variable, t_stat, p.value, slope, pass, trend, group)
write_csv(stat_table2, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/r2-recombination_trend.csv")


# average and sd of r2 using rolling average dataframe from above

r2 <- dataframe %>% filter(stat == "r2_norm") %>% 
  group_by(tax_abbr_name, tax_abbr_name.col) %>% 
  summarise(mean = median(linkage),
            sd = sd(linkage))
r2 <- left_join(r2, cyano.mOTUs, by = c("tax_abbr_name", "tax_abbr_name.col"))

# average and sd of r2 using all data
dataframe <- instrain_92_linkage 

dataframe <- dataframe %>% filter(tax_abbr_name %in% filt_cyanos)
dataframe$count <- 1
count <- dataframe %>% select(tax_abbr_name,SampleDate,count) %>% distinct() %>% group_by(tax_abbr_name) %>% summarize(count_samples = sum(count))

dataframe$distance <- dataframe %>% with(distance + 1)
r2_1 <- dataframe %>% na.omit() %>% filter(distance <= 250) %>% 
  group_by(tax_abbr_name,tax_abbr_name.col,group) %>% 
  summarize(mean = mean(r2_normalized),
            sd = sd(r2_normalized),
            count = sum(count))

r2_1 <- left_join(r2_1, count, by = "tax_abbr_name")



# who has the steepest slopes
p2 <- fitted_models %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  mutate(across(group, ~factor(., levels = c("Pseudanabaena","Nodosilinea","Filamentous-1","Microcystis","Snowella","Vulcanococcus","Cyanobium","Solitary-1")))) %>% 
  ggplot(aes(x = tax_abbr_name.col, y = slope, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Pseudanabaena" = "#4cd505",
                                "Nodosilinea" = "#3cb306",
                                "Filamentous-1" = "#236a00",
                                "Snowella" = "#f4659d",
                                "Microcystis" = "#e39f0c",
                                "Cyanobium" = "#00dede",
                                "Solitary-1" = "#008bff",
                                "Vulcanococcus" = "#7a05c6")) +
  labs(title = "Slope of SNP linkage (d') across 5 - 250 bp", x = "Cyano-mOTUs", y = "d' normalized slope", fill = "") +
  theme_bw()+
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        plot.title = element_text(size = 15, face = "bold", vjust = -7), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        #aspect.ratio = 1,
        legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dataframe <- dataframe_backup
r2_range <- dataframe %>% filter(distance == 5 | distance == 250) %>% select(c("tax_abbr_name", "tax_abbr_name.col", "distance", "r2_norm"))
r2_range <- r2_range %>% pivot_wider(values_from = r2_norm, names_from = distance)
r2 <- left_join(r2, r2_range, by = c("tax_abbr_name", "tax_abbr_name.col"))

stat_table2 <- r2 %>% select(tax_abbr_name, mean, sd, `5.x`, `250.x`, group)
write_csv(stat_table2, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/r2-recombination_trend.csv")

p3 <- r2 %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  mutate(across(group, ~factor(., levels = c("Pseudanabaena","Nodosilinea","Filamentous-1","Microcystis","Snowella","Vulcanococcus","Cyanobium","Solitary-1")))) %>% 
  ggplot(aes(x = tax_abbr_name.col, y = mean, fill = group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x=tax_abbr_name.col, ymin=`250`, ymax=`5`), linewidth=0.4, colour="black", size=1.3) +
  scale_fill_manual(values = c("Pseudanabaena" = "#4cd505",
                               "Nodosilinea" = "#3cb306",
                               "Filamentous-1" = "#236a00",
                               "Snowella" = "#f4659d",
                               "Microcystis" = "#e39f0c",
                               "Cyanobium" = "#00dede",
                               "Solitary-1" = "#008bff",
                               "Vulcanococcus" = "#7a05c6")) +
  labs(title = "Average SNP linkage (r2) across 5 - 250 bp", x = "Cyano-mOTUs", y = "AVG r2 normalized", fill = "") +
  theme_bw()+
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        plot.title = element_text(size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        #aspect.ratio = 1,
        legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)))




layout <- "
ABB
ACC
"
p1+p2+p3+ plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))




p3 <- fitted_models2 %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  ggplot(aes(x = tax_abbr_name.col, y = slope, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Pseudanabaena" = "#4cd505",
                               "Nodosilinea" = "#3cb306",
                               "Filamentous-1" = "#236a00",
                               "Snowella" = "#f4659d",
                               "Microcystis" = "#e39f0c",
                               "Cyanobium" = "#00dede",
                               "Solitary-1" = "#008bff",
                               "Vulcanococcus" = "#7a05c6")) +
  labs(x = "Cyano-mOTUs", y = "r2 normalized slope", fill = "") +
  theme_bw()+
  theme(axis.text.x = ggtext::element_markdown(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=16, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 5)),
        plot.title = element_text(size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, "pt"),
        #aspect.ratio = 1,
        legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# keep fitted models
rm(count, fitted_models2, linear_fit, linear_fit2, p1,p2,p3, r2, r2_1, r2_range, 
   stat_table, stat_table2, dist, dataframe, dataframe_filt, dataframe_backup)

### Fun triangle plot 

# dprime slope
tri_slopedprime <- fitted_models %>% select(variable, tax_abbr_name.col, slope)
tri_slopedprime <- left_join(tri_slopedprime, count, by = c("variable" = "tax_abbr_name"))

# avg nucleotide diversity and irep
tri_pi.irep <- instrain_genome %>% filter(tax_abbr_name %in% filt_cyanos) %>% 
  select(tax_abbr_name, tax_abbr_name.col, group, morphology, SampleDate, nucl_diversity_92ANI, iRep_92ANI, coverage_92ANI, breadth_92ANI)
tri_pi.irep$count <- 1

tri_pi.irep <- tri_pi.irep %>% filter(coverage_92ANI >= 10 & breadth_92ANI >= 0.8)

tri_pi <- tri_pi.irep %>% group_by(tax_abbr_name, tax_abbr_name.col, group, morphology) %>% 
  summarise(mean_pi = mean(nucl_diversity_92ANI),
            count.pi = sum(count))

tri_irep <- tri_pi.irep %>% na.omit(iRep_92ANI) %>% 
  group_by(tax_abbr_name, tax_abbr_name.col, group, morphology) %>% 
  summarise(mean_irep = mean(iRep_92ANI, na.rm = TRUE),
            count = sum(count))

stat_table <- left_join(tri_slopedprime, tri_pi, by = c("variable" = "tax_abbr_name", "tax_abbr_name.col"))
stat_table <- left_join(stat_table, tri_irep, by = c("variable" = "tax_abbr_name", "tax_abbr_name.col", "group", "morphology"))
write_csv(stat_table, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/trenary_data.csv")


# normalize (min max)

# custom function to implement min max scaling
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#normalise data using custom function
tri_slopedprime["slope"] <- lapply(tri_slopedprime["slope"], minMax)
tri_pi["mean_pi"] <- lapply(tri_pi["mean_pi"], minMax)
tri_irep["mean_irep"] <- lapply(tri_irep["mean_irep"], minMax)

tri <- left_join(tri_slopedprime, tri_pi, by =c("variable" = "tax_abbr_name", "tax_abbr_name.col"))
tri <- left_join(tri, tri_irep, by = c("variable" = "tax_abbr_name", "tax_abbr_name.col", "group", "morphology"))


library(ggtern)
tri_plot <- tri %>% #filter(morphology == "Solitary") %>% 
  mutate(across(group, ~factor(., levels = c("Pseudanabaena","Nodosilinea","Filamentous-1","Microcystis","Snowella","Vulcanococcus","Cyanobium","Solitary-1")))) %>% 
  ggtern(aes(z=mean_pi,y=mean_irep, x=slope)) +
  scale_x_continuous(breaks=seq(-100,0,20)) +
  geom_point(aes(color = group), size = 2) +
  geom_text(aes(label = variable), vjust=1) +
  labs(title="", z = "AVG π", y = "AVG iRep", x = "d' slope", color = "") +
  scale_color_manual(values = c("Pseudanabaena" = "#4cd505",
                               "Nodosilinea" = "#3cb306",
                               "Filamentous-1" = "#236a00",
                               "Snowella" = "#f4659d",
                               "Microcystis" = "#e39f0c",
                               "Cyanobium" = "#00dede",
                               "Solitary-1" = "#008bff",
                               "Vulcanococcus" = "#7a05c6")) +
  theme_rgbg(base_size = 18, base_family = "") +
  theme(legend.text = element_text(size = 15),
        legend.position = "bottom",
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, "pt"),
        title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

tri <- left_join(tri, count, by = c("variable" = "tax_abbr_name"))


tri_plot <- tri %>% #filter(morphology == "Solitary") %>% 
  mutate(across(group, ~factor(., levels = c("Pseudanabaena","Nodosilinea","Filamentous-1","Microcystis","Snowella","Vulcanococcus","Cyanobium","Solitary-1")))) %>% 
  ggtern(aes(y=mean_pi,z=mean_irep, x=slope)) +
  geom_point(aes(color = group), size = 2) +
  #geom_text(aes(label = variable), vjust=1) +
  labs(title="", y = "AVG π", z = "AVG iRep", x = "d' slope", color = "") +
  scale_color_manual(values = c("Pseudanabaena" = "#4cd505",
                                "Nodosilinea" = "#3cb306",
                                "Filamentous-1" = "#236a00",
                                "Snowella" = "#f4659d",
                                "Microcystis" = "#e39f0c",
                                "Cyanobium" = "#00dede",
                                "Solitary-1" = "#008bff",
                                "Vulcanococcus" = "#7a05c6")) +
  theme_rgbg(base_size = 18, base_family = "") +
  theme(legend.text = element_text(size = 15),
        legend.position = "bottom",
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, "pt"),
        title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))



