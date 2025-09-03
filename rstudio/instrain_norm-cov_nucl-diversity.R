######################################################
### Krys Kibler 
### Select for final (analyzed) cyano-mOTUs
######################################################


######################################################
library(tidyverse)


######################################################
### Taxonomy
cyano.mOTUs <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/data/files/cyanos_genome_summary.csv"
cyano.mOTUs <- read_delim(cyano.mOTUs)


# create genera column
cyano.mOTUs$genera <- gsub("d__Bacteria;.+g__", "", cyano.mOTUs$taxonomy)
cyano.mOTUs$genera <- paste("g__", cyano.mOTUs$genera, sep = "")
cyano.mOTUs$genera <- gsub(";s.+", "", cyano.mOTUs$genera)
cyano.mOTUs$genera <- gsub("g__", "", cyano.mOTUs$genera)

######################################################
### instrain dataframe

d_instrain_genome <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/output/all_genome_info.tsv"
d_metadata <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/data/files/all_JGI_metadata_471_samples.txt"


instrain_genome <- read_delim(d_instrain_genome)
metadata <- read_delim(d_metadata)

metadata <- metadata %>% select(SampleDate, SampleID)
metadata$SampleDate <- substr(metadata$SampleDate, 3, 11)
metadata$SampleDate <- dmy(metadata$SampleDate)
metadata <- metadata %>% na.omit()

instrain_genome$sampleID <- gsub("cyano_76-vs-", "", instrain_genome$sampleID)
instrain_genome$sampleID <- gsub(".IS_genome_info", "", instrain_genome$sampleID)

instrain_genome <- left_join(instrain_genome, metadata, by = c("sampleID" = "SampleID"))

instrain_genome$genome <- gsub(".fna", "", instrain_genome$genome)

# Full dataframe for instrain genome level diversity
instrain_genome <- left_join(cyano.mOTUs, instrain_genome, by = c("tax_abbr_name" = "genome"))

instrain_genome <- instrain_genome %>% filter(morphology != "parasite")
######################################################
### normalized coverage from instrain_genome

# total read pairs from cyano in samples
metag_mapped <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_reads = sum(filtered_read_pair_count))

mean_mapped = mean(metag_mapped$sum_reads) # avg mapped number of read pairs = 3046701

instrain_genome <- left_join(instrain_genome, metag_mapped, by = "SampleDate")

# normalized coverage calc
instrain_genome$norm_coverage <- instrain_genome %>% with(coverage * (3046701 / sum_reads))
instrain_genome$norm_ratio <- instrain_genome %>% with(3046701 / sum_reads)

# relative abundance calc
instrain_genome$RA_by.reads <- instrain_genome %>% with((filtered_read_pair_count / sum_reads) * 100)

RA <- instrain_genome %>% 
  filter(norm_coverage >= 5) %>% 
  group_by(SampleDate) %>% 
  summarise(sum_norm.cov = sum(norm_coverage))
instrain_genome <- left_join(instrain_genome, RA, by = "SampleDate")

instrain_genome$RA_by.cov <- instrain_genome %>% with((norm_coverage / sum_norm.cov) * 100)


### Figure of normalized coverage for filter cyano-mOTUs

genome <- c("PSEUDA_13","PSEUDA_70","PSEUDA_281","PSEUDA_31","PSEUDA_157",
  "NODOS_103","NODOS_105","MCYST_56","MCYST_62","MCYST_31",
  "MCYST_2","CYANO_45","CYANO_51_1","CYANO_106","CYANO_84",
  "APHAN_134","DOLIS_105","DOLIS_187","VULCA_20","VULCA_120",
  "VULCA_28","VULCA_96","CYBIM_76","CYBIM_15","CYBIM_77",
  "CYBIM_22","CYBIM_90","CYBIM_39","CYBIM_157","CYBIM_104",
  "CYBIM_89","CYBIM_60_1","CYBIM_101","CYBIM_93","CYBIM_63",
  "CYBIM_73_1","CYBIM_190","CYBIM_200","CYBIM_51","CYBIM_31",
  "CYBIM_119_1","CYBIM_80","CYBIM_130_1","CYBIM_28","CYBIM_160",
  "CYBIM_52","CYBIM_136_1","CYBIM_94","CYBIM_119","CYBIM_73")


col <- c("#52f000","#52f000","#52f000","#52f000","#52f000",
            "#52f000","#52f000","#f4659d","#e39f0c","#e39f0c",
            "#e39f0c","#308e00","#308e00","#308e00","#308e00",
            "#308e00","#308e00","#308e00","#7a05c6","#7a05c6",
            "#7a05c6","#7a05c6","#008bff","#008bff","#008bff",
            "#008bff","#008bff","#008bff","#008bff","#008bff",
            "#008bff","#008bff","#008bff","#008bff","#008bff",
            "#008bff","#008bff","#00dede","#00dede","#00dede",
            "#00dede","#00dede","#00dede","#00dede","#00dede",
            "#00dede","#00dede","#00dede","#00dede","#00dede")


cols <- as.data.frame(cbind(genome,col))

instrain_genome <- left_join(instrain_genome, cols, by = c("tax_abbr_name"="genome"))

instrain_genome$tax_abbr_name.col<- paste0("<span style=\"color: ", instrain_genome$col, "\">", instrain_genome$tax_abbr_name, "</span>")

levels <- c("<span style=\"color: #52f000\">PSEUDA_13</span>","<span style=\"color: #52f000\">PSEUDA_70</span>","<span style=\"color: #52f000\">PSEUDA_281</span>","<span style=\"color: #52f000\">PSEUDA_31</span>","<span style=\"color: #52f000\">PSEUDA_157</span>",
            "<span style=\"color: #52f000\">NODOS_103</span>","<span style=\"color: #52f000\">NODOS_105</span>","<span style=\"color: #f4659d\">MCYST_56</span>","<span style=\"color: #e39f0c\">MCYST_62</span>","<span style=\"color: #e39f0c\">MCYST_31</span>",
            "<span style=\"color: #e39f0c\">MCYST_2</span>","<span style=\"color: #308e00\">CYANO_45</span>","<span style=\"color: #308e00\">CYANO_51_1</span>","<span style=\"color: #308e00\">CYANO_106</span>","<span style=\"color: #308e00\">CYANO_84</span>",
            "<span style=\"color: #308e00\">APHAN_134</span>","<span style=\"color: #308e00\">DOLIS_105</span>","<span style=\"color: #308e00\">DOLIS_187</span>","<span style=\"color: #7a05c6\">VULCA_20</span>","<span style=\"color: #7a05c6\">VULCA_120</span>",
            "<span style=\"color: #7a05c6\">VULCA_28</span>","<span style=\"color: #7a05c6\">VULCA_96</span>","<span style=\"color: #008bff\">CYBIM_76</span>","<span style=\"color: #008bff\">CYBIM_15</span>","<span style=\"color: #008bff\">CYBIM_77</span>",
            "<span style=\"color: #008bff\">CYBIM_22</span>","<span style=\"color: #008bff\">CYBIM_90</span>","<span style=\"color: #008bff\">CYBIM_39</span>","<span style=\"color: #008bff\">CYBIM_157</span>","<span style=\"color: #008bff\">CYBIM_104</span>",
            "<span style=\"color: #008bff\">CYBIM_89</span>","<span style=\"color: #008bff\">CYBIM_60_1</span>","<span style=\"color: #008bff\">CYBIM_101</span>","<span style=\"color: #008bff\">CYBIM_93</span>","<span style=\"color: #008bff\">CYBIM_63</span>",
            "<span style=\"color: #008bff\">CYBIM_73_1</span>","<span style=\"color: #008bff\">CYBIM_190</span>","<span style=\"color: #00dede\">CYBIM_200</span>","<span style=\"color: #00dede\">CYBIM_51</span>","<span style=\"color: #00dede\">CYBIM_31</span>",
            "<span style=\"color: #00dede\">CYBIM_119_1</span>","<span style=\"color: #00dede\">CYBIM_80</span>","<span style=\"color: #00dede\">CYBIM_130_1</span>","<span style=\"color: #00dede\">CYBIM_28</span>","<span style=\"color: #00dede\">CYBIM_160</span>",
            "<span style=\"color: #00dede\">CYBIM_52</span>","<span style=\"color: #00dede\">CYBIM_136_1</span>","<span style=\"color: #00dede\">CYBIM_94</span>","<span style=\"color: #00dede\">CYBIM_119</span>","<span style=\"color: #00dede\">CYBIM_73</span>")


# Sum RA by mOTU
RA_mOTU <- instrain_genome %>% group_by(tax_abbr_name, tax_abbr_name.col) %>% 
  summarise(RA_sum = sum(RA_by.cov))

### Add number of samples to boxplot 
instrain_genome$count <- 1

count <- instrain_genome %>% 
  filter(norm_coverage >= 5) %>% 
  group_by(tax_abbr_name.col) %>% 
  summarise(count_above5cov = sum(count))

instrain_genome <- left_join(instrain_genome, count, by = "tax_abbr_name.col")


library(ggtext)
# ggplot boxplot for normalized coverage
p1 <- instrain_genome %>% filter(norm_coverage >= 5) %>% 
  ggplot(aes(x = fct_rev(factor(tax_abbr_name.col, levels = levels)), y = norm_coverage)) +
  geom_boxplot() +
  geom_text(aes(label = count_above5cov, y = 225, x = fct_rev(factor(tax_abbr_name.col, levels = levels)))) +
  coord_flip() +
  labs(x = "Cyano-mOTUs", y = "Normalized Coverage") +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = ggtext::element_markdown(size = 12),
        axis.title.x = element_text(size=20, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=20, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.position = "right",
        title = element_text(size = 30, face="bold"),
        plot.margin = margin(t = 10, r = 0, b = 10, l = 5))

library(scales)
# ggplot bar chart for summed RA
p2 <- RA_mOTU %>% 
  ggplot(aes(x = fct_rev(factor(tax_abbr_name.col, levels = levels)), y = RA_sum)) +
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
        plot.margin = margin(t = 10, r = 10, b = 10, l = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border=element_blank(), axis.line.x=element_blank())



library(patchwork)

layout <- "
AAAB
AAAB
"
png(filename="/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/final-plots/norm_coverage.png", width=700, height=1200)

p1+p2+plot_layout(design = layout)

dev.off()

# save data, remove dataframes captured already in instrain_genome
rm(metag_mapped)
rm(metadata)
rm(cols)
rm(RA_mOTU)
rm(p1,p2)

# save data, remove values that were used to load erased dataframes or create dfs
rm(col)
rm(d_instrain_genome,d_metadata)
rm(mean_mapped)
rm(layout)


##########################ND############################ 
### Nucleotide diversity 
instrain_genome_plotdf <- instrain_genome %>% 
  filter(coverage >= 10) %>% 
  filter(breadth_minCov >= 0.8)

instrain_genome_plotdf$count <- 1
count <- instrain_genome_plotdf %>% 
  group_by(tax_abbr_name) %>% 
  summarise(count_samps10cov = sum(count))

instrain_genome_plotdf <- left_join(instrain_genome_plotdf, count, by = "tax_abbr_name")
rm(count)

# Statistics
# https://statsandr.com/blog/anova-in-r/ #link for normal anova
# https://statsandr.com/blog/kruskal-wallis-test-nonparametric-version-anova/ for non-normal anova

# anova with normal distribution assumption
res_aov <- aov(nucl_diversity ~ morphology,
               data = instrain_genome_plotdf)

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
leveneTest(nucl_diversity ~ morphology,
           data = instrain_genome_plotdf) 

# Anova with non-normal distribution assumption
kruskal.test(nucl_diversity ~ morphology,
             data = instrain_genome_plotdf)
#	Kruskal-Wallis rank sum test
#data:  nucl_diversity by morphology
#Kruskal-Wallis chi-squared = 370.09, df = 3, p-value < 2.2e-16

# post-hoc Dunn Test
library(FSA)
stats <- dunnTest(nucl_diversity ~ morphology,
         data = instrain_genome_plotdf,
         method = "holm")
stats <- stats$res
stats <- stats %>% mutate(across(Comparison, ~factor(., levels = c("Chain - Colonial", "Chain - Single", 
                                                                   "Branched - Colonial", "Branched - Single",
                                                                   "Colonial - Single", "Branched - Chain")))) %>% arrange(Comparison)
stats <- stats %>% separate(Comparison, c("group1", "group2"), sep = " - ")
stars <- c("****", "****", "**", "ns", "****", "**")
stats <- cbind(stats, stars)


#Dunn (1964) Kruskal-Wallis multiple comparison
#p-values adjusted with the Holm method.
#Comparison          Z      P.unadj        P.adj
#1    Branched - Chain   3.4192565 6.279250e-04 1.883775e-03 .0018 **
#2 Branched - Colonial  -2.9550796 3.125884e-03 6.251768e-03 .0062 **
#3    Chain - Colonial -17.9278171 7.153352e-72 4.292011e-71 .000000 ****
#4   Branched - Single   0.4997803 6.172298e-01 6.172298e-01 .62 ns
#5      Chain - Single -13.3893666 6.977324e-41 3.488662e-40 .000000 ****
#6   Colonial - Single  10.4233220 1.940544e-25 7.762174e-25 .000000 ****

library(ggsignif)

instrain_genome_plotdf <- instrain_genome_plotdf %>% mutate(across(morphology, ~factor(., levels=c("Chain","Colonial","Single","Branched")))) 

instrain_genome_plotdf %>% group_by(morphology) %>% 
  summarize(mean = mean(nucl_diversity), sd = sd(nucl_diversity))

#morphology    mean      sd
#<fct>        <dbl>   <dbl>
#1 Chain      0.00353 0.00177
#2 Colonial   0.00968 0.00270
#3 Single     0.00611 0.00388
#4 Branched   0.00661 0.00347

instrain_genome_plotdf %>% group_by(clade) %>% 
  summarize(mean = mean(nucl_diversity), sd = sd(nucl_diversity))

#clade            mean      sd
#<chr>           <dbl>   <dbl>
#1 Chain-1       0.00376 0.00228
#2 Chain-2       0.00338 0.00135
#3 Cyanobium-1   0.00417 0.00266
#4 Cyanobium-2   0.00450 0.00234
#5 Microcystis   0.00968 0.00270
#6 Snowella      0.00661 0.00347
#7 Vulcanococcus 0.00978 0.00368

# boxplot between morphologies
p1 <- instrain_genome_plotdf %>% 
  ggplot(aes(x = morphology, y = nucl_diversity, color = morphology)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.5) +
  geom_violin(width=1, alpha = 0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.8, outlier.shape = NA) +
  stat_pvalue_manual(stats, label = "stars", y.position = .019, 
                     step.increase = 0.04, tip.length = 0.01) +
  scale_color_manual(values = c("Chain" = "#0ab053",
                               "Colonial" = "#e39f0c",
                               "Single" = "#7660f6",
                               "Branched" = "#f4659d")) +
  labs(y = "Nucleotide Diversity (π)", x = "Cyano-mOTUs Morphology") +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold"),
                 axis.text.y = element_text(size=15),
                 axis.title.x = element_text(size=26, face="bold", margin = margin(t = -20, r = 0, b = 0, l = 0)),
                 axis.title.y = element_text(size=26, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = "none")
  
  
# boxplot chain clades
p2 <- instrain_genome_plotdf %>% filter(morphology == "Chain") %>% 
  ggplot(aes(x = clade, y = nucl_diversity, fill = clade)) +
  geom_boxplot() +
  labs(x = "Chain Clades", y = "π") +
  scale_y_continuous(limits = c(0,0.02)) +
  scale_fill_manual(values = c("Chain-1" = "#52f000",
                               "Chain-2" = "#308e00")) +
  #geom_signif(map_signif_level = TRUE, comparisons = list(c("Chain-1, Chain-2"))) +
  stat_compare_means(method = "t.test", label.y = 0.013, label = "p.signif", size = 12) +      # Add global p-value
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none")

# boxplot of chain mOTUs
p3 <- instrain_genome_plotdf %>% filter(morphology == "Chain") %>% 
  ggplot(aes(x = fct_reorder(tax_abbr_name,nucl_diversity,median), y = nucl_diversity, fill = clade)) +
  geom_boxplot() +
  labs(x = "'Chain' cyano-mOTUs", y = "π") +
  scale_y_continuous(limits = c(0,0.02)) +
  scale_fill_manual(values = c("Chain-1" = "#52f000",
                                "Chain-2" = "#308e00")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none")

# boxplot of colonial mOTUs
p4 <- instrain_genome_plotdf %>% filter(morphology == "Colonial") %>% 
  ggplot(aes(x = fct_reorder(tax_abbr_name,nucl_diversity,median), y = nucl_diversity, fill = clade)) +
  geom_boxplot() +
  labs(x = "'Colonial' cyano-mOTUs", y = "π") +
  scale_y_continuous(limits = c(0,0.025)) +
  scale_fill_manual(values = c("Microcystis" = "#e39f0c")) +
  stat_compare_means(label.y = .02, label = "p.signif", size = 7) +      # Add global p-value
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none")

library(ggsignif)
# boxplot of single clades
p5 <- instrain_genome_plotdf %>% filter(morphology == "Single") %>% 
  ggplot(aes(x = fct_reorder(clade,nucl_diversity,median), y = nucl_diversity, fill = clade)) +
  geom_boxplot() +
  labs(x = "Single Clades", y = "π") +
  scale_y_continuous(limits = c(0,0.025)) +
  scale_fill_manual(values = c("Cyanobium-1" = "#008bff",
                               "Cyanobium-2" = "#00dede",
                               "Vulcanococcus" = "#7a05c6")) +
  geom_signif(comparisons = list(c("Cyanobium-1","Cyanobium-2"),
                                  c("Cyanobium-1","Vulcanococcus"),
                                  c("Cyanobium-2","Vulcanococcus")),
              map_signif_level = TRUE, textsize = 6, test = "wilcox.test",
              margin_top = 0.01,step_increase = 0.12, tip_length = 0.01, vjust = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none")

# boxplot of single mOTUs
p6 <- instrain_genome_plotdf %>% filter(morphology == "Single") %>% 
  ggplot(aes(x = fct_reorder(tax_abbr_name,nucl_diversity,median), y = nucl_diversity, fill = clade)) +
  geom_boxplot() +
  labs(x = "'Single' cyano-mOTUs", y = "π") +
  scale_y_continuous(limits = c(0,0.02)) +
  scale_fill_manual(values = c("Cyanobium-1" = "#008bff",
                               "Cyanobium-2" = "#00dede",
                               "Vulcanococcus" = "#7a05c6")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none")


# https://patchwork.data-imaginist.com/articles/guides/layout.html
layout <- "
AABCD
AAEEE
AAFFF
"
png(filename="/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/final-plots/norm_coverage.png", width=700, height=1200)

p1+p2+p5+p4+p3+p6+plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

dev.off()


rm(p1,p2,p3,p4,p5,p6)
rm(res_aov,stats)
rm(layout)
rm(stars)


#######################iRep###############################
### Nucleotide diversity and iRep

instrain_genome_plotdf %>% filter(tax_abbr_name == "APHAN_134") %>% 
  filter(coverage >= 10) %>% 
  ggplot(aes(x = iRep, y = nucl_diversity)) +
  geom_point()

instrain_genome_plotdf %>% filter(tax_abbr_name == "MCYST_2") %>% 
  filter(coverage >= 10) %>% 
  ggplot(aes(x = iRep, y = nucl_diversity)) +
  geom_point()


# setup iRep focused dataframe
instrain_genome_plotdf$month <- month(instrain_genome_plotdf$SampleDate)
instrain_genome_plotdf_irep <- instrain_genome_plotdf %>% na.omit(iRep)
count_irep <- instrain_genome_plotdf_irep %>% group_by(tax_abbr_name) %>% 
  summarise(count_irep = sum(count))
instrain_genome_plotdf_irep <- left_join(instrain_genome_plotdf_irep, count_irep, by = "tax_abbr_name")


# color points by month
instrain_genome_plotdf_irep %>% #filter(tax_abbr_name == "APHAN_134") %>% 
  #filter(count_samps10cov >= 20) %>% 
  filter(count_irep >= 10) %>% 
  #filter(month >= 5 &
  #         month <= 9) %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, color = month)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free") +
  scale_color_viridis_c() +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


# Set the threshold value
threshold <- 100

# Create a new variable that truncates the size values
instrain_genome_plotdf_irep$truncated_size <- ifelse(instrain_genome_plotdf_irep$coverage >= threshold, threshold, instrain_genome_plotdf_irep$coverage)


### mOTUs plotted individually
p9 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(count_irep >= 10) %>% 
  filter(tax_abbr_name != "MCYST_2" &
           tax_abbr_name != "APHAN_134") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 4) +
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
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


p1 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(count_irep >= 10) %>% 
  filter(tax_abbr_name == "MCYST_2") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 4) +
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
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white"))

p2 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(count_irep >= 10) %>% 
  filter(tax_abbr_name == "APHAN_134") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, scales = "free", nrow = 4) +
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
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 16),
        strip.background = element_rect(fill="white"))




### Grouped by category
instrain_genome_plotdf_irep %>% 
  # filter(cat_figureH != "Chain") %>% 
  filter(coverage >= 10) %>% 
  #filter(count_samps10cov >= 20) %>% 
  filter(count_irep >= 10) %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, size = truncated_size)) +
  geom_point() +
  facet_wrap(~clade) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 200), 
                        breaks = c(10, 25, 50, 75, 100, 200),
                        labels = c("10x", "25x", "50x", "75x", "100x", ">200x")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", size = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


### Chains
# Chain-1
p1 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Chain-1") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "Nucleotide Diversity (π)", alpha = "Coverage", title = "Chain-1") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(vjust = -7))


p2 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Chain-1") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, size = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, nrow = 1) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 700), 
                        breaks = c(10, 25, 50, 75, 100, 200),
                        labels = c("10x", "25x", "50x", "75x", "100x", ">200x")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", size = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))

# Chain-2
p3 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Chain-2") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage", title = "Chain-2") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"),
        plot.title = element_text(vjust = -7))


p4 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Chain-2") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, size = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, nrow = 1) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 700), 
                        breaks = c(10, 25, 50, 75, 100, 200),
                        labels = c("10x", "25x", "50x", "75x", "100x", ">200x")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", size = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


# Microcystis
p5 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Microcystis") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage", title = "Microcystis") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


p6 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(clade == "Microcystis" |
           clade == "Snowella") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, size = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, nrow = 1) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 700), 
                        breaks = c(10, 25, 50, 75, 100, 200),
                        labels = c("10x", "25x", "50x", "75x", "100x", ">200x")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", size = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


### Single
p7 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(morphology == "Single") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, alpha = truncated_size)) +
  geom_point() +
  scale_alpha_continuous(range  = c(0.1, 6), 
                         limits = c(0, 700), 
                         breaks = c(10, 25, 50, 100),
                         labels = c("10x", "25x", "50x", ">100x")) +
  labs(x = "Log(iRep)", y = "π", alpha = "Coverage", title = "Single Clades") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))


p8 <- instrain_genome_plotdf_irep %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  filter(count_irep >= 10) %>% 
  filter(morphology == "Single") %>% 
  ggplot(aes(x = log10(iRep), y = log10(nucl_diversity), size = truncated_size)) +
  geom_point() +
  facet_wrap(~tax_abbr_name.col, nrow = 2) +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 700), 
                        breaks = c(10, 25, 50, 75, 100, 200),
                        labels = c("10x", "25x", "50x", "75x", "100x", ">200x")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", size = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"),
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))



layout <- "
ABEEEEE
CDEEEEE
"

png(filename="/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/final-plots/norm_coverage.png", width=700, height=1200)

p1+p3+p5+p7+p9+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

dev.off()

### More impressive figure
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


p10 <- seq %>% filter(seq > 0 ) %>% 
  filter(calc_exp >= -5 &
           seq <= 0.5) %>% 
  ggplot(aes(x = seq, y = calc_exp)) +
  geom_line(linewidth = 2) +
  labs(x = "Replication Index (iRep)", y = "Nucleotide Diversity (π)", title = "Diversity Maintenance") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        plot.title = element_text(vjust = -7, size = 15, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
p11 <- seq %>% filter(seq < 0 ) %>% 
  filter(calc_exp <= 5 &
           seq >= -0.5) %>% 
  ggplot(aes(x = seq, y = calc_exp)) +
  geom_line(linewidth = 2) +
  labs(x = "iRep", y = "π", title ="Clonal Expansion") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(vjust = -7, size = 15, face = "bold"), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



library(patchwork)
layout <- "
AC
BD
EE
EE
"

png(filename="/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/final-plots/iRep_pi.png", width=1000, height=2000)

p10+p11+p1+p2+p9+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

dev.off()

rm(p1,p2,p3,p4,p5,p6,p7,p8,p9)
rm(count_irep)
rm(instrain_genome_plotdf_irep)
rm(threshold, layout)

# coverage and nucl diversity
instrain_genome_plotdf %>% 
  filter(coverage >= 10) %>% 
  ggplot(aes(x = log10(coverage), y = nucl_diversity, color = factor(cat_figureH, levels = c("Chain","Colonial","Branched","Cyanobium","Vulcanococcus", "Other")))) +
  geom_point() +
  scale_color_manual(values = c("Chain" = "#0ab053",
                                "Colonial" = "#e39f0c",
                                "Branched" = "#f4659d",
                                "Cyanobium" = "#4c2ef1",
                                "Vulcanococcus" = "#7660f6",
                                "Other" = "gray")) +
  labs(x = "log_10(iRep)", y = "Nucleotide Diversity (π)", color = "Cyano-mOTU Categories") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"))


###########################inter###########################
### Interannual nucleotide diversity variances 

# minimum-presence applied

p1 <- instrain_genome_plotdf %>% 
  mutate(across(tax_abbr_name, ~factor(., levels=genome))) %>% 
  filter(count_samps10cov >= 20) %>% 
  ggplot(aes(x = SampleDate, y = nucl_diversity)) +
  geom_line() +
  geom_point() +
  scale_x_date(date_labels = "'%y", breaks = "2 year") +
  facet_wrap(~tax_abbr_name.col, nrow = 5) +
  labs(x = "Year", y = "Nucleotide Diversity (π)", color = "Morphology") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        #legend.title = element_text(size = 15),
        #legend.text = element_text(size = 14), 
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

p2 <- instrain_genome_plotdf %>% 
  filter(genera != "Snowella") %>% 
  mutate(across(clade, ~factor(., levels=c("Chain-1", "Chain-2", "Microcystis","Snowella","Vulcanococcus","Cyanobium-1", "Cyanobium-2")))) %>% 
  filter(count_samps10cov >= 20) %>%
  ggplot(aes(x = SampleDate, y = nucl_diversity, color = clade)) +
  geom_line(aes(group = tax_abbr_name, alpha = 0.5, color = "darkgrey")) +
  geom_point(aes(group = tax_abbr_name, alpha = 0.5)) +
  geom_smooth(se=TRUE, color = "black", size=1) +
  facet_wrap(~clade, nrow = 3, scales = "free") +
  scale_color_manual(values = c("Chain-1" = "#52f000",
                                "Chain-2" = "#308e00",
                                "Microcystis" = "#e39f0c",
                                "Cyanobium-1" = "#008bff",
                                "Cyanobium-2" = "#00dede",
                                "Vulcanococcus" = "#7a05c6")) +
  scale_x_date(date_labels = "'%y", breaks = "year") +
  labs(x = "Year", y = "Nucleotide Diversity (π)") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=22, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=22, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        strip.text = element_text(size=15, face = "bold"),
        strip.background = element_rect(fill="white"),
        plot.margin = unit(c(0.25, 0.33, 0.25, 0.25), "inches"))

library("plotly")
ggplotly(p2)

# linear regression models
library(broom)

linear_fit <- instrain_genome_plotdf %>% 
  filter(count_samps10cov >= 20) %>% 
  group_by(tax_abbr_name) %>% 
  do(tidy(lm(nucl_diversity~SampleDate, data = .))) %>% 
  select(variable = tax_abbr_name, term, t_stat = statistic, p.value, slope = estimate)


fitted_models <- linear_fit
fitted_models <- mutate(fitted_models, pass = ifelse(p.value < 0.05, 'pass', 'fail'))
fitted_models <- mutate(fitted_models, trend = ifelse(slope > 0, 'pos', 'neg'))

fitted_models <- left_join(fitted_models, cyano.mOTUs, by = c("variable" = "tax_abbr_name"))
fitted_models$category <- paste(fitted_models$pass, fitted_models$trend, sep = "-")

fitted_models$category <- gsub("fail-neg", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("fail-pos", "No-sig change", fitted_models$category)
fitted_models$category <- gsub("NA-NA", "No-sig change", fitted_models$category)

fitted_models$category <- gsub("pass-neg", "Negative, Sig", fitted_models$category)
fitted_models$category <- gsub("pass-pos", "Positive, Sig", fitted_models$category)

fitted_models$count <- 1

numeric <- fitted_models %>% filter(term == "SampleDate") %>% 
  #filter(trend == "neg") %>% 
  filter(pass == "fail")
max(numeric$slope)
fitted_models %>% filter(term == "SampleDate") %>% 
  group_by(category) %>% 
  summarise(sum = sum(count))
#pass-pos = 6
#pass-neg = 7
#no-trend = 8

p3 <- fitted_models %>% filter(term == "SampleDate") %>% 
  mutate(across(clade, ~factor(., levels=c("Chain-1", "Chain-2", "Microcystis","Snowella","Vulcanococcus","Cyanobium-1", "Cyanobium-2")))) %>% 
  ggplot(aes(x = factor(category, levels = c("Positive, Sig", "Negative, Sig", "No-sig change")), y = count, fill = clade)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Chain-1" = "#52f000",
                               "Chain-2" = "#308e00",
                               "Microcystis" = "#e39f0c",
                               "Snowella" = "#f4659d",
                               "Cyanobium-1" = "#008bff",
                               "Cyanobium-2" = "#00dede",
                               "Vulcanococcus" = "#7a05c6")) +
  labs(y = "No of cyanomOTUs", x = "Interannual Trend", fill = "Cyano-mOTU Clades") +
  scale_x_discrete(labels= c("Pos", "Neg", "No change")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 3),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

# prepare sup table
fitted_models <- fitted_models %>% filter(term == "SampleDate") %>% 
  select(c(1, 3:7, 13))
write_csv(fitted_models, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/interannual_trend.csv")

# interannual variation
variation <- instrain_genome_plotdf %>%
  filter(count_samps10cov >= 20) %>%
  group_by(tax_abbr_name, tax_abbr_name.col, clade) %>% 
  summarise(mean = mean(nucl_diversity), sd = sd(nucl_diversity))

variation$CV <- variation %>% with((sd / mean) * 100)

library(scales)
p4 <- variation %>% 
  ggplot(aes(x = factor(tax_abbr_name.col, level = levels), y = CV)) +
  geom_rect(aes(ymin=13.5, ymax=45, xmin=0.5, xmax=9.5), alpha = .01) +
  geom_rect(aes(ymin=3.5, ymax=33, xmin=9.5, xmax=12.5), alpha = .02) +
  geom_rect(aes(ymin=29.5, ymax=67.5, xmin=12.5, xmax=21.5), alpha = .03) +
  geom_point(aes(size = 2)) +
  labs(x = "cyano-mOTUs", y = "π Coeff of Variation") +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        strip.text = element_text(size=15, face = "bold"),
        plot.margin = unit(c(0.25, 0.33, 0.25, 0.25), "inches"))


layout <- "
AAA
AAA
BCC
"

p2+p3+p4+plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))


















######################################################
### Intra-annual nucleotide diversity 

# variability nucl diversity within a year using CV
instrain_genome_plotdf$year4 <- year(instrain_genome_plotdf$SampleDate)

count_yearsum <- instrain_genome_plotdf %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(coverage >= 10) %>% 
  group_by(tax_abbr_name, year4) %>% 
  summarise(count_samps10cov_byyear = sum(count))
instrain_genome_plotdf <- left_join(instrain_genome_plotdf, count_yearsum, by = c("tax_abbr_name", "year4"))

instrain_genome_plotdf_coeff <- instrain_genome_plotdf %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(coverage >= 10) %>% 
  filter(count_samps10cov_byyear >= 5) %>% 
  group_by(tax_abbr_name, tax_abbr_name.col, morphology, cat_figureH, year4, count_samps10cov_byyear) %>% 
  summarise(mean_nucl.div = mean(nucl_diversity), sd_nucl.div = sd(nucl_diversity))

instrain_genome_plotdf_coeff$CV <- instrain_genome_plotdf_coeff %>% with((sd_nucl.div / mean_nucl.div) * 100)
instrain_genome_plotdf_coeff$count <- 1

count_yearsumtot <- instrain_genome_plotdf_coeff %>% 
  group_by(tax_abbr_name) %>% 
  summarise(count_year = sum(count))

instrain_genome_plotdf_coeff <- left_join(instrain_genome_plotdf_coeff, count_yearsumtot, by = "tax_abbr_name")


my_comparisons <- list( c("Chain", "Colonial"), c("Chain", "Cyanobium"), c("Chain", "Vulcanococcus"),
                        c("Cyanobium", "Vulcanococcus"))

axis_labels <- c("Chain (7)", "Colonial (3)", "Cyanobium (9)", "Vulcanococcus (3)", "Other (5)")
p1 <- instrain_genome_plotdf_coeff %>%
  filter(count_year >= 3) %>% 
  ggplot(aes(x = factor(cat_figureH, levels=c("Chain","Colonial","Cyanobium", "Vulcanococcus", "Branched", "Other")), y = CV, fill = cat_figureH)) + 
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", vjust = 0.5,
                     label.y = c(75, 80, 85, 90)) + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 80) +    # Add global p-value
  #stat_compare_means(method = "anova", label.y = 55) +
  labs(x = "Cyano mOTUs", y = "Intrannual Pi Variation (CV)") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = c("Chain" = "#0ab053",
                               "Colonial" = "#e39f0c",
                               "Cyanobium" = "#4c2ef1",
                               "Vulcanococcus" = "#7660f6",
                               "Other" = "gray")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"))

p2 <- instrain_genome_plotdf_coeff %>%
  filter(count_year >= 3) %>% 
  ggplot(aes(x = tax_abbr_name.col, y = CV)) + 
  geom_boxplot() +
  geom_point() +
  labs(x = "Cyano mOTUs", y = "Intrannual Pi Variation (CV)") +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"))

p3 <- instrain_genome_plotdf_coeff %>%
  filter(count_year >= 5) %>% 
  ggplot(aes(x = fct_reorder(tax_abbr_name,CV,median), y = CV, fill = cat_figureH)) + 
  geom_boxplot() +
  geom_point() +
  labs(x = "Cyano mOTUs", y = "Intrannual Pi Variation (CV)") +
  scale_fill_manual(values = c("Chain" = "#0ab053",
                               "Colonial" = "#e39f0c",
                               "Cyanobium" = "#4c2ef1",
                               "Vulcanococcus" = "#7660f6",
                               "Other" = "gray")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"))

ggarrange(p1,p3, ncol = 2, labels = c("A", "B"), align = "h")


### Intra-annual trend

instrain_genome_plotdf$jday <- yday(instrain_genome_plotdf$SampleDate)


instrain_genome_plotdf %>% filter(count_samps10cov_byyear >= 5) %>% 
  filter(tax_abbr_name == "APHAN_134" |
           tax_abbr_name == "MCYST_2") %>% 
  ggplot(aes(x = jday, y = nucl_diversity, group = year4)) +
  geom_line() +
  facet_wrap(~tax_abbr_name) +
  theme_bw()
  












instrain_genome_plotdf %>% 
  filter(coverage >= 10) %>% 
  filter(count_samps10cov >= 20) %>% 
  filter(morphology == "Chain") %>% 
  ggplot(aes(x = log10(iRep), y = nucl_diversity, color = genera)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "right",
        plot.margin = unit(c(0.25, 0.25, 0.1, 0.1), "inches"))
  


instrain_genome %>% 
  filter(tax_abbr_name == "DOLIS_105" |
           tax_abbr_name == "DOLIS_187" |
           tax_abbr_name == "APHAN_134") %>% 
  ggplot(aes(x = SampleDate, y = log10(coverage), groups = year4)) +
  geom_point() +
  geom_line() +
  scale_x_date(breaks = "year") +
  facet_wrap(~tax_abbr_name, nrow=3)
