######################################################
### Krys Kibler
### How variable are the SCGs compared to the rest of the genome

library(tidyverse)
library(readr)
library(dplyr)

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

instrain_genome$count <- 1
count <- instrain_genome %>% filter(coverage >= 10) %>% group_by(tax_abbr_name) %>% 
  summarise_at(c("count"), sum, na.rm = TRUE)
colnames(count) <- c("tax_abbr_name", "count_samps10cov")
instrain_genome <- left_join(instrain_genome, count, by = "tax_abbr_name")

count <- instrain_genome %>% filter(coverage >= 20) %>% group_by(tax_abbr_name) %>% 
  summarise_at(c("count"), sum, na.rm = TRUE)
colnames(count) <- c("tax_abbr_name", "count_samps20cov")
instrain_genome <- left_join(instrain_genome, count, by = "tax_abbr_name")

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



######################################################
d_anvio.genefxns <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene_functions.txt"
d_anvio.genecalls <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene-calls.txt"

anvio.genefxns <- read_delim(d_anvio.genefxns)
anvio.genecalls <- read_delim(d_anvio.genecalls)

anvio.genefxns <- left_join(anvio.genefxns, select(anvio.genecalls, c(1:7)), by = "gene_callers_id")
anvio.genefxns$genome <- gsub(".Contig.+", "", anvio.genefxns$contig)

anvio.genefxns_function <- anvio.genefxns %>% filter(source == "COG20_FUNCTION")
anvio.genefxns_category <- anvio.genefxns %>% filter(source == "COG20_CATEGORY")

anvio.genefxns <- left_join(anvio.genefxns_function, select(anvio.genefxns_category, c(1,3,4)), by = "gene_callers_id")
anvio.genefxns$end <- anvio.genefxns %>% with(stop -1)

anvio.genefxns$direction <- gsub("r", "-1", anvio.genefxns$direction)
anvio.genefxns$direction <- gsub("f", "1", anvio.genefxns$direction)
anvio.genefxns$direction <- as.numeric(anvio.genefxns$direction)

# cogs id'ed in APHAN_134
APHAN_strong_cogs <- c("COG0016","COG0048","COG0049","COG0051","COG0052",
                       "COG0080","COG0087","COG0089","COG0090","COG0091",
                       "COG0092","COG0102","COG0103","COG0130","COG0185",
                       "COG0198","COG0256","COG0504")

# cogs id'ed in MCYST_2
MCYST_strong_cogs <- c("COG0048","COG0049","COG0052","COG0080","COG0081",
                       "COG0087","COG0089","COG0090","COG0091","COG0092",
                       "COG0093","COG0094","COG0100","COG0103","COG0184",
                       "COG0185","COG0186","COG0198","COG0201","COG0244",
                       "COG0541")


APHAN_strong_cogs.loc <- anvio.genefxns %>% filter(genome == "APHAN_134") %>% 
  filter(accession.x %in% APHAN_strong_cogs) %>% 
  filter(source == "COG20_FUNCTION")

MCYST_strong_cogs.loc <- anvio.genefxns %>% filter(genome == "MCYST_2") %>% 
  filter(accession.x %in% MCYST_strong_cogs) %>% 
  filter(source == "COG20_FUNCTION")


######################################################
d_instrain_gene <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/output/all_gene_info.tsv"
instrain_gene <- read_delim(d_instrain_gene)

instrain_gene$sampleID <- gsub("cyano_76-vs-", "", instrain_gene$sampleID)
instrain_gene$sampleID <- gsub(".IS_gene_info", "", instrain_gene$sampleID)

instrain_gene <- left_join(instrain_gene, metadata, by = c("sampleID" = "SampleID"))

instrain_gene$genome <- gsub(".Contig.+", "", instrain_gene$scaffold)

# Full dataframe for instrain genome level diversity
instrain_gene <- left_join(cyano.mOTUs, instrain_gene, by = c("tax_abbr_name" = "genome"))

instrain_gene <- instrain_gene %>% filter(morphology != "parasite")

instrain_gene <- left_join(instrain_gene, anvio.genefxns, by = c("tax_abbr_name"="genome", "scaffold" = "contig", "start", "end", "direction"))


######################################################
### Mean and sd nucl diversities

# strong cogs
# Aphan
instrain_gene_aph.strong.cogs <- instrain_gene %>% 
  filter(tax_abbr_name == "APHAN_134") %>% 
  filter(accession.x %in% APHAN_strong_cogs)
instrain_gene_aph.strong.cogs_mean.sd <- instrain_gene_aph.strong.cogs %>% filter(coverage >= 10) %>% 
  group_by(SampleDate) %>% 
  summarise(nucl_diversity=mean(nucl_diversity))
instrain_gene_aph.strong.cogs_mean.sd$category <- "STRONG"
instrain_gene_aph.strong.cogs_mean.sd$tax_abbr_name <- "APHAN_134"

# Mcyst
instrain_gene_mcyst.strong.cogs <- instrain_gene %>% 
  filter(tax_abbr_name == "MCYST_2") %>% 
  filter(accession.x %in% MCYST_strong_cogs)
instrain_gene_mcyst.strong.cogs_mean.sd <- instrain_gene_mcyst.strong.cogs %>% filter(coverage >= 10) %>% 
  group_by(SampleDate) %>% 
  summarise(nucl_diversity=mean(nucl_diversity))
instrain_gene_mcyst.strong.cogs_mean.sd$category <- "STRONG"
instrain_gene_mcyst.strong.cogs_mean.sd$tax_abbr_name <- "MCYST_2"

# global
instrain_genome_mean.sd <- instrain_genome %>% 
  filter(coverage >= 10) %>% 
  filter(breadth_minCov >= 0.8) %>% 
  filter(count_samps10cov >= 20) %>% 
  select(tax_abbr_name, SampleDate, nucl_diversity)
  #group_by(tax_abbr_name) %>% 
  #summarize(mean=mean(nucl_diversity), sd=sd(nucl_diversity))
instrain_genome_mean.sd$category <- "Global"

# just SCGs (NCBI COG J)
instrain_gene_J.mean.sd <- instrain_gene %>% 
  filter(accession.y == "J") %>% 
  filter(coverage >= 10) %>% 
  group_by(tax_abbr_name, SampleDate) %>% 
  summarize(nucl_diversity=mean(nucl_diversity))
instrain_gene_J.mean.sd$category <- "NCBI_COG - J"

# all genes not J
instrain_gene_mean.sd <- instrain_gene %>% 
  filter(accession.y != "J") %>% 
  filter(coverage >= 10) %>% 
  group_by(tax_abbr_name, SampleDate) %>% 
  summarize(nucl_diversity=mean(nucl_diversity))
instrain_gene_mean.sd$category <- "Gene"



dataframe <- rbind(instrain_genome_mean.sd, instrain_gene_mean.sd, instrain_gene_J.mean.sd, instrain_gene_aph.strong.cogs_mean.sd, instrain_gene_mcyst.strong.cogs_mean.sd)


######################################################
dataframe %>% ggplot(aes(x = category, y = mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_point(size=2) +
  facet_wrap(~tax_abbr_name)

p1 <- dataframe %>% mutate(across(category, ~factor(., levels=c("Global", "Gene", "NCBI_COG - J", "STRONG")))) %>% 
  filter(tax_abbr_name == "APHAN_134") %>% 
  ggplot(aes(x = category, y = nucl_diversity)) +
  geom_boxplot() +
  stat_compare_means(label.y = .015, label = "p.signif", size = 7) + 
  labs(x = "", y = "Nucleotide Diversity (π)", title = "APHAN_134") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(0,1,0,0), "cm"),
        plot.title = element_text(vjust = -6))

p2 <- dataframe %>% mutate(across(category, ~factor(., levels=c("Global", "Gene", "NCBI_COG - J", "STRONG")))) %>% 
  filter(tax_abbr_name == "MCYST_2") %>% 
  ggplot(aes(x = category, y = nucl_diversity)) +
  geom_boxplot() +
  stat_compare_means(label.y = .015, label = "p.signif", size = 7) + 
  labs(x = "", y = "π", title = "MCYST_2") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(0,1,0,0), "cm"))

genomes_interest <- c("CYANO_106","CYANO_45","CYANO_51_1","CYBIM_101","CYBIM_104","CYBIM_119","CYBIM_119_1",
"CYBIM_157","CYBIM_200","CYBIM_31","CYBIM_63","CYBIM_73","CYBIM_73_1","CYBIM_89","CYBIM_90",
"MCYST_31","MCYST_56","MCYST_62","NODOS_105","PSEUDA_13","PSEUDA_70","VULCA_20",
"VULCA_28","VULCA_96")

genomes_interest <- unique(instrain_genome_mean.sd$tax_abbr_name)

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
dataframe <- left_join(dataframe, cols, by = c("tax_abbr_name" = "genome"))

dataframe$tax_abbr_name.col<- paste0("<span style=\"color: ", dataframe$col, "\">", dataframe$tax_abbr_name, "</span>")

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


library(ggtext)
library(ggpubr)
p3 <- dataframe %>% mutate(across(tax_abbr_name.col, ~factor(., levels=levels))) %>% 
  mutate(across(category, ~factor(., levels=c("Global", "Gene", "NCBI_COG - J")))) %>%
  filter(tax_abbr_name %in% genomes_interest) %>% 
  filter(tax_abbr_name != "APHAN_134" &
           tax_abbr_name != "MCYST_2") %>% 
  ggplot(aes(x = category, y = nucl_diversity)) +
  geom_boxplot() +
  stat_compare_means(label.y = .015, label = "p.signif", size = 7) + 
  facet_wrap(~tax_abbr_name.col, ncol = 6) +
  labs(x = "", y = "π", title = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15, face = "bold", margin = margin(t = 5, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size=15, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 10)),
        legend.position = "none",
        strip.text = ggtext::element_markdown(size = 12),
        strip.background = element_rect(fill="white"))

library(patchwork)

layout <- "
ACCCC
BCCCC
"

p1+p2+p3+plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))
