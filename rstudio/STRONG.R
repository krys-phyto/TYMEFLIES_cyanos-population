######################################################
### Krys Kibler
### Strong


######################################################
library(tidyverse)
library(readr)
library(dplyr)

######################################################
### load the strong data

df_MCYST.2_strong <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/STRONG/run1/MCYST_2/mag_and_haplo_cov.tsv"
df_APHAN.134_strong <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/STRONG/run1/APHAN_134/mag_and_haplo_cov.tsv"

MCYST.2_strong <- read_delim(df_MCYST.2_strong)
APHAN.134_strong <- read_delim(df_APHAN.134_strong)

MCYST.2_strong <- MCYST.2_strong %>% pivot_longer(-c(1), values_to = "STRONG_cov", names_to = "SampleID")
APHAN.134_strong <- APHAN.134_strong %>% pivot_longer(-c(1), values_to = "STRONG_cov", names_to = "SampleID")

# 81 missing data points for MCYST_2
# 96 missing data points for APHAN_134

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
# pull from diamond-50_instrain_tidy.R

dataframe <- instrain_genome %>% select(tax_abbr_name, tax_abbr_name.col, SampleDate, norm_coverage_92, nucl_diversity_92ANI)

######################################################
# Prep for figures
MCYST.2_strong <- MCYST.2_strong %>% na.omit()
MCYST.2_strong <- left_join(MCYST.2_strong, metadata, by = "SampleID")
MCYST.2_strong$year4 <- year(MCYST.2_strong$SampleDate)
colnames(MCYST.2_strong) <- c("MCYST.2_haplotypes", "SampleID", "STRONG_cov", "SampleDate", "year4")

APHAN.134_strong <- APHAN.134_strong %>% na.omit()
APHAN.134_strong <- left_join(APHAN.134_strong, metadata, by = "SampleID")
APHAN.134_strong$year4 <- year(APHAN.134_strong$SampleDate)
colnames(APHAN.134_strong) <- c("APHAN.134_haplotypes", "SampleID", "STRONG_cov", "SampleDate", "year4")

# Fix leading zeroes in haplotype names (option to shorten as well)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Bin_11_", "", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("haplo_", "Haplo_", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_1", "Haplo_01", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_2", "Haplo_02", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_3", "Haplo_03", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_4", "Haplo_04", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_5", "Haplo_05", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_6", "Haplo_06", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_7", "Haplo_07", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_8", "Haplo_08", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_9", "Haplo_09", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_010", "Haplo_10", APHAN.134_strong$APHAN.134_haplotypes)
APHAN.134_strong$APHAN.134_haplotypes <- gsub("Haplo_011", "Haplo_11", APHAN.134_strong$APHAN.134_haplotypes)


MCYST.2_strong$MCYST.2_haplotypes <- gsub("Bin_2_", "", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("haplo_", "Haplo_", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_1", "Haplo_01", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_2", "Haplo_02", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_3", "Haplo_03", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_4", "Haplo_04", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_5", "Haplo_05", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_6", "Haplo_06", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_7", "Haplo_07", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_8", "Haplo_08", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_9", "Haplo_09", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_010", "Haplo_10", MCYST.2_strong$MCYST.2_haplotypes)
MCYST.2_strong$MCYST.2_haplotypes <- gsub("Haplo_011", "Haplo_11", MCYST.2_strong$MCYST.2_haplotypes)



dataframe$year4 <- year(dataframe$SampleDate)

######################################################
library(scales)

p1 <- dataframe %>% filter(tax_abbr_name == "MCYST_2") %>% 
  filter(year4 == 2015) %>% 
   ggplot(aes(x = SampleDate, y = norm_coverage_92)) +
  geom_line(aes(color = "Normalized Coverage")) +
  geom_point(aes(color = "Normalized Coverage")) +
  geom_line(aes(y = nucl_diversity_92ANI *1000, color = "Nucleotide Diversity"), linetype = "dashed") +
  geom_point(aes(y = nucl_diversity_92ANI*1000, color = "Nucleotide Diversity")) +
  scale_y_continuous(sec.axis = sec_axis(~./1250, 
                                         name = "Nucleotide Diversity")) +
  scale_x_date(limits = as.Date(c("2015-06-19","2015-08-03"))) +
             #  breaks = "1 month",
             # labels = date_format("%b")) +
  labs(x = "", y = "Normalized Coverage", color = "", size = "iRep", title = "MCYST_2") +
  scale_color_manual(values = c("gray30","orange2")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0, 0.25), "inches"))

p2 <- instrain_genome %>% filter(tax_abbr_name == "APHAN_134") %>% 
  filter(year4 == 2015) %>% 
  ggplot(aes(x = SampleDate, y = norm_coverage_92)) +
  geom_point(aes(color = "Normalized Coverage")) +
  geom_line(aes(color = "Normalized Coverage")) +
  geom_line(aes(y = nucl_diversity_92ANI*10000, color = "Nucleotide Diversity"), linetype = "dashed") +
  geom_point(aes(y = nucl_diversity_92ANI*10000, color = "Nucleotide Diversity")) +
  scale_y_continuous(sec.axis = sec_axis(~./10000, name = "Nucleotide Diversity")) +
  scale_x_date(limits = as.Date(c("2015-06-04","2015-07-21"))) +
                  # breaks = "1 month",
                  # labels = date_format("%b")) +
  labs(x = "", y = "Normalized Coverage", color = "", size = "iRep", title = "APHAN_134") +
  scale_color_manual(values = c("gray30","#0ab053")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0, 0.25), "inches"))



p3 <- MCYST.2_strong %>% filter(year4 == 2015) %>% 
  filter(MCYST.2_haplotypes != "Bin_2") %>% 
  ggplot(aes(x = as.factor(SampleDate), y = STRONG_cov, fill = MCYST.2_haplotypes)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "2015", y = "STRONG weighted coverage", fill = "MCYST_2 Haplotypes", title = "MCYST_2") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))


p4 <- APHAN.134_strong %>% filter(year4 == 2015) %>% 
  filter(APHAN.134_haplotypes != "Bin_11") %>% 
  ggplot(aes(x = as.factor(SampleDate), y = STRONG_cov, fill = APHAN.134_haplotypes)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "2015", y = "STRONG weighted coverage", fill = "APHAN_134 Haplotypes", title = "APHAN_134") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

layout <- "
A
B
B
"

# combination stations
patch1 <- wrap_elements(p1 / free(p3) + plot_layout(design = layout))
patch2 <- wrap_elements(p2 / free(p4) + plot_layout(design = layout))

layout <- "
AB
"
patch1+patch2+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))




instrain_genome %>% filter(tax_abbr_name == "APHAN_134") %>% 
  ggplot(aes(x = SampleDate, y = norm_coverage)) +
  geom_line() +
  geom_point(aes(size = iRep))


instrain_genome %>% filter(tax_abbr_name == "MCYST_2") %>% 
  ggplot(aes(x = SampleDate, y = norm_coverage)) +
  geom_line() +
  geom_point(aes(size = iRep))


######################################################
### Blooms

instrain_genome_APHAN <- instrain_genome %>% filter(tax_abbr_name == "APHAN_134") 
instrain_genome_MCYST <- instrain_genome %>% filter(tax_abbr_name == "MCYST_2")

instrain_genome_APHAN <- instrain_genome_APHAN %>% group_by(year4) %>% 
  slice(which.max(coverage))

instrain_genome_MCYST <- instrain_genome_MCYST %>% group_by(year4) %>% 
  slice(which.max(coverage))

APHAN_blooms <- unique(instrain_genome_APHAN$SampleDate)
MCYST_blooms <- unique(instrain_genome_MCYST$SampleDate)

p1 <- instrain_genome %>% filter(tax_abbr_name == "MCYST_2") %>% 
  filter(SampleDate %in% MCYST_blooms) %>% 
  ggplot(aes(x = as.factor(SampleDate), y = norm_coverage)) +
  #geom_line(aes(color = "Normalized Coverage")) +
  geom_point(aes(color = "Normalized Coverage", size = 2)) +
  #geom_line(aes(y = nucl_diversity*1000, color = "Nucleotide Diversity"), linetype = "dashed") +
  #geom_point(aes(y = nucl_diversity*1000, color = "Nucleotide Diversity")) +
  #scale_y_continuous(sec.axis = sec_axis(~./3000, 
  #                                       name = "Nucleotide Diversity")) +
  labs(x = "", y = "Normalized Coverage", color = "", size = "iRep", title = "MCYST_2") +
  scale_color_manual(values = c("gray30","orange2")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, angle = 45, hjust = 1),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0, 0.25), "inches"))

p2 <- instrain_genome %>% filter(tax_abbr_name == "APHAN_134") %>% 
  filter(SampleDate %in% APHAN_blooms) %>% 
  ggplot(aes(x = as.factor(SampleDate), y = norm_coverage)) +
  #geom_line(aes(color = "Normalized Coverage")) +
  geom_point(aes(color = "Normalized Coverage", size = 2)) +
  #geom_line(aes(y = nucl_diversity*1000, color = "Nucleotide Diversity"), linetype = "dashed") +
  #geom_point(aes(y = nucl_diversity*1000, color = "Nucleotide Diversity")) +
  #scale_y_continuous(sec.axis = sec_axis(~./3000, 
  #                                       name = "Nucleotide Diversity")) +
  labs(x = "", y = "Normalized Coverage", color = "", size = "iRep", title = "APHAN_134") +
  scale_color_manual(values = c("gray30","#0ab053")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, angle = 45, hjust = 1),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0, 0.25), "inches"))

p3 <- MCYST.2_strong %>% filter(SampleDate %in% MCYST_blooms) %>% 
  filter(MCYST.2_haplotypes != "Bin_2") %>% 
  ggplot(aes(x = as.factor(SampleDate), y = STRONG_cov, fill = MCYST.2_haplotypes)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "2015", y = "STRONG weighted coverage", fill = "MCYST_2 Haplotypes", title = "MCYST_2") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))


p4 <- APHAN.134_strong %>% filter(SampleDate %in% APHAN_blooms) %>%
  filter(APHAN.134_haplotypes != "Bin_11") %>% 
  ggplot(aes(x = as.factor(SampleDate), y = STRONG_cov, fill = APHAN.134_haplotypes)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "2015", y = "STRONG weighted coverage", fill = "APHAN_134 Haplotypes", title = "APHAN_134") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=15, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size = 18),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))

layout <- "
AC
BD
BD
"
p1+p3+p2+p4+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))


######################################################
### Species Evenness

library(vegan)

MCYST.2_strong_wide <- MCYST.2_strong %>% filter(MCYST.2_haplotypes != "Bin_2") %>% 
