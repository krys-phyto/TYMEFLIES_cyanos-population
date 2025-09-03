#################################
# DNA mechanisms present absence to go next to tree
# Krys Kibler 2025-03-04
#################################


#################################
# Libraries the basics
library(tidyverse)
library(readr)
library(dplyr)

library(ggdendro)
library(grid)
#################################


#################################

d_anvio.genefxns <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene_functions.txt"
d_anvio.genecalls <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene-calls.txt"

anvio.genefxns <- read_delim(d_anvio.genefxns)
anvio.genecalls <- read_delim(d_anvio.genecalls)

anvio.genefxns <- left_join(anvio.genefxns, select(anvio.genecalls, c(1,2)), by = "gene_callers_id")
anvio.genefxns$genome <- gsub(".Contig.+", "", anvio.genefxns$contig)


Uvr <- anvio.genefxns[grepl("Uvr", anvio.genefxns$`function`),]
RecA <- anvio.genefxns[grepl("RecA", anvio.genefxns$`function`),]
RecB <- anvio.genefxns[grepl("RecB", anvio.genefxns$`function`),]
RecC <- anvio.genefxns[grepl("RecC", anvio.genefxns$`function`),]
RecD <- anvio.genefxns[grepl("RecD", anvio.genefxns$`function`),]
RecF <- anvio.genefxns[grepl("RecF", anvio.genefxns$`function`),]
RecG <- anvio.genefxns[grepl("RecG", anvio.genefxns$`function`),]
RecJ <- anvio.genefxns[grepl("RecJ", anvio.genefxns$`function`),]
RecN <- anvio.genefxns[grepl("RecN", anvio.genefxns$`function`),]
RecQ <- anvio.genefxns[grepl("RecQ", anvio.genefxns$`function`),]
RecR <- anvio.genefxns[grepl("RecR", anvio.genefxns$`function`),]
RecO <- anvio.genefxns[grepl("RecO", anvio.genefxns$`function`),]

Ruv <- anvio.genefxns[grepl("RuvABC", anvio.genefxns$`function`),]


#MutH <- anvio.genefxns[grepl("", anvio.genefxns$`function`),]     # no mutH
MutL <- anvio.genefxns[grepl("MutL", anvio.genefxns$`function`),]
#MutM <- anvio.genefxns[grepl("", anvio.genefxns$`function`),]     # no mutM
MutS <- anvio.genefxns[grepl("MutS", anvio.genefxns$`function`),] # both S1 and S2
MutT <- anvio.genefxns[grepl("MutT", anvio.genefxns$`function`),]
MutY <- anvio.genefxns[grepl("MutY", anvio.genefxns$`function`),]

# no umuCD or dinB

Ssb <- anvio.genefxns[grepl("Ssb", anvio.genefxns$`function`),]

Com <- anvio.genefxns[grepl("DNA uptake", anvio.genefxns$`function`),]

PhrB <- anvio.genefxns[grepl("PhrB", anvio.genefxns$`function`),]
SplB <- anvio.genefxns[grepl("SplB", anvio.genefxns$`function`),]

AlkB <- anvio.genefxns[grepl("AlkB", anvio.genefxns$`function`),]
AlkC <- anvio.genefxns[grepl("AlkC", anvio.genefxns$`function`),]
AlkD <- anvio.genefxns[grepl("AlkD", anvio.genefxns$`function`),]

Xer <- anvio.genefxns[grepl("Xer", anvio.genefxns$`function`),]

AdaA <- anvio.genefxns[grepl("AdaA", anvio.genefxns$`function`),]
AdaB <- anvio.genefxns[grepl("AdaB", anvio.genefxns$`function`),]

RadA <- anvio.genefxns[grepl("RadA/Sms", anvio.genefxns$`function`),]
Rad52 <- anvio.genefxns[grepl("RAD52", anvio.genefxns$`function`),]

transposase <- anvio.genefxns[grepl("Transposase", anvio.genefxns$`function`, ignore.case = TRUE),]

dna_mech <- rbind(Uvr, RecA, RecB, RecC, RecD, RecF, RecG, RecJ, RecN, RecO, RecQ, RecR, 
      Ruv, MutL, MutS, MutT, MutY, 
      Ssb, Com, PhrB, SplB, AlkB, AlkC, AlkD, Xer, 
      AdaA, AdaB, RadA, Rad52, transposase)

dna_mech <- dna_mech %>% distinct()
write_csv(dna_mech, file = "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/anvio_dna-mech.csv")

dna_mech_plot <- read_csv("~/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/anvio_dna-mech_EDIT.csv")


dna_mech_plot$accession_short <- gsub("!!!.+", "*", dna_mech_plot$accession)
dna_mech_plot$count <- 1


dna_mech_plot_sum <- dna_mech_plot %>% group_by(genome, gene_sh, gene_facet) %>% 
  summarise(sum_count = sum(count))


# from instrain_norm-cov_nucl-diversity
cyano.mOTUs # already prefiltered

cyano.mOTUs <- cyano.mOTUs %>% filter(clade != "Vampire")

###################################################
# Present/Absence

dna_mech <- left_join(dna_mech, select(dna_mech_plot, c("function", "accession", "genome", "gene_sh", "gene_facet", "gene_callers_id")), by = c("genome", "accession", "function", "gene_callers_id"))
dna_mech <- left_join(cyano.mOTUs, dna_mech, by = c("tax_abbr_name" = "genome"))


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

dna_mech <- left_join(dna_mech, cols, by = c("tax_abbr_name"="genome"))

dna_mech$tax_abbr_name.col<- paste0("<span style=\"color: ", dna_mech$col, "\">", dna_mech$tax_abbr_name, "</span>")

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

dna_mech$count <- 1

p1 <- dna_mech %>% filter(gene_facet == "Alk") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Alk", y = "Cyano mOTUs") + theme.full

p2 <- dna_mech %>% filter(gene_facet == "DNA Uptake") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "DNA Uptake", y = "Cyano mOTUs") + theme.full

p4 <- dna_mech %>% filter(gene_facet == "Mut") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Mut", y = "Cyano mOTUs") + theme.skim

p6 <- dna_mech %>% filter(gene_facet == "RAD") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "RAD", y = "Cyano mOTUs") + theme.skim

p7 <- dna_mech %>% filter(gene_facet == "Rec") %>% 
  filter(gene_sh != "RepA") %>%
  filter(gene_sh != "RepA*") %>% 
  filter(gene_sh != "Cas4") %>% 
  filter(gene_sh != "Cas4*") %>% 
  #subset(!grepl("COG", gene_sh)) %>% 
  ggplot(aes(x = factor(gene_sh, levels = c("LexA", "RecA", "RecB", "RecB*", "NucS", "COG2251", "COG3372",
                                            "RecC", "RecD", "RecF", "RecF*", "RecG", "RecJ", "RecN", "RecO",
                                            "RecQ", "RecQ*", "RecR")),
             y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Rec", y = "Cyano mOTUs") + theme.skim

p8 <- dna_mech %>% filter(gene_facet == "Ruv") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Ruv", y = "Cyano mOTUs") + theme.skim

interest <- c("PSEUDA_13", "PSEUDA_70", "DOLIS_187", "APHAN_134", "MCYST_162_1", "MCYST_62", "MCYST_79", "MCYST_100", 
              "MCYST_2", "MCYST_162", "MCYST_24", "MCYST_31", "MCYST_162_2", "MCYST_157", "MCYST_56", "CYANO_51_1")
p9 <- dna_mech %>%
  #filter(genome %in% interest) %>% 
  filter(gene_facet == "Transposases") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Transposases", y = "Cyano mOTUs") + theme.skim

p10 <- dna_mech %>% filter(gene_facet == "Uvr") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Uvr", y = "Cyano mOTUs") + theme.skim

p11 <- dna_mech %>% filter(gene_facet == "Xer") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Xer", y = "Cyano mOTUs") + theme.skim


p3 <- dna_mech %>% filter(gene_facet == "Endonuclease") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Endonuclease", y = "Cyano mOTUs") + theme.skim
  

p5 <- dna_mech %>% filter(gene_facet == "Photolyase") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels), fill = count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Photolyase", y = "Cyano mOTUs") + theme.skim


layout <- "
ABCCDEF
"

p2+p4+p7+p8+p10+p11+plot_layout(design = layout)


theme.full <-   theme_bw() +
  theme(axis.text.x = element_text(size=15, angle =90, vjust = 0.5),
        axis.text.y = ggtext::element_markdown(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "inches"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        title = element_text(size = 16, face = "bold"))

theme.skim <-   theme_bw() +
  theme(axis.text.x = element_text(size=15, angle =90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "inches"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        title = element_text(size = 16, face = "bold"))


dna_mech_table <- dna_mech %>% select("accession", "function", "gene_sh", "gene_facet")
dna_mech_table <- distinct(dna_mech_table)
dna_mech_table <- dna_mech_table %>% filter(gene_facet == "DNA Uptake" |
                                              gene_facet == "Mut" |
                                              gene_facet == "Rec" |
                                              gene_facet == "Ruv" |
                                              gene_facet == "Uvr" |
                                              gene_facet == "Xer")
write_csv(dna_mech_table, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/Paper/dna_mech_genes.csv")








dna_mech_plot_sum <- left_join(dna_mech_plot_sum, cyano.mOTUs, by = c("genome" = "tax_abbr_name"))


d_instrain_genome.ts_avg <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/generated_dfs/cyano_instrain-genome_ts-avgs.csv"
instrain_genome.ts_avg <- read.csv(d_instrain_genome.ts_avg)

dna_mech_plot_sum <- left_join(dna_mech_plot_sum, instrain_genome.ts_avg, by = "genome") 


levels = c("PSEUDA_13","PSEUDA_70","PSEUDA_281","PSEUDA_31","PSEUDA_157",
           "NODOS_103","NODOS_105","MCYST_56","MCYST_62","MCYST_31","MCYST_2",
           "CYANO_45","CYANO_51_1","CYANO_106","CYANO_84","APHAN_134",
           "DOLIS_105","DOLIS_187","VULCA_20","VULCA_120","VULCA_28","VULCA_96",
           "CYBIM_200","CYBIM_51","CYBIM_31","CYBIM_119_1","CYBIM_80","CYBIM_130_1",
           "CYBIM_28","CYBIM_160","CYBIM_52","CYBIM_136_1","CYBIM_94","CYBIM_119",
           "CYBIM_73","CYBIM_76","CYBIM_15","CYBIM_77","CYBIM_22","CYBIM_90",
           "CYBIM_39","CYBIM_157","CYBIM_104","CYBIM_89","CYBIM_60_1","CYBIM_101",
           "CYBIM_93","CYBIM_63","CYBIM_73_1","CYBIM_190")

col <- c("#0ab053","#0ab053","#0ab053","#0ab053","#0ab053",
         "#0ab053","#0ab053","#f4659d","#e39f0c","#e39f0c","#e39f0c",
         "#0ab053","#0ab053","#0ab053","#0ab053","#0ab053",
         "#0ab053","#0ab053","#5166f0","#5166f0","#5166f0","#5166f0",
         "#5166f0","#5166f0","#5166f0","#5166f0","#5166f0","#5166f0",
         "#5166f0","#5166f0","#5166f0","#5166f0","#5166f0","#5166f0",
         "#5166f0","#5166f0","#5166f0","#5166f0","#5166f0","#5166f0",
         "#5166f0","#5166f0","#5166f0","#5166f0","#5166f0","#5166f0",
         "#5166f0","#5166f0","#5166f0","#5166f0")

cols <- as.data.frame(cbind(levels,col))

cols$tax_abbr_name.col<- paste0("<span style=\"color: ", cols$col, "\">", cols$levels, "</span>")

levels_2 = c("<span style=\"color: #0ab053\">PSEUDA_13</span>","<span style=\"color: #0ab053\">PSEUDA_70</span>","<span style=\"color: #0ab053\">PSEUDA_281</span>","<span style=\"color: #0ab053\">PSEUDA_31</span>","<span style=\"color: #0ab053\">PSEUDA_157</span>",
             "<span style=\"color: #0ab053\">NODOS_103</span>","<span style=\"color: #0ab053\">NODOS_105</span>","<span style=\"color: #f4659d\">MCYST_56</span>","<span style=\"color: #e39f0c\">MCYST_62</span>","<span style=\"color: #e39f0c\">MCYST_31</span>","<span style=\"color: #e39f0c\">MCYST_2</span>",
             "<span style=\"color: #0ab053\">CYANO_45</span>","<span style=\"color: #0ab053\">CYANO_51_1</span>","<span style=\"color: #0ab053\">CYANO_106</span>","<span style=\"color: #0ab053\">CYANO_84</span>","<span style=\"color: #0ab053\">APHAN_134</span>",
             "<span style=\"color: #0ab053\">DOLIS_105</span>","<span style=\"color: #0ab053\">DOLIS_187</span>","<span style=\"color: #5166f0\">VULCA_20</span>","<span style=\"color: #5166f0\">VULCA_120</span>","<span style=\"color: #5166f0\">VULCA_28</span>","<span style=\"color: #5166f0\">VULCA_96</span>",
             "<span style=\"color: #5166f0\">CYBIM_200</span>","<span style=\"color: #5166f0\">CYBIM_51</span>","<span style=\"color: #5166f0\">CYBIM_31</span>","<span style=\"color: #5166f0\">CYBIM_119_1</span>","<span style=\"color: #5166f0\">CYBIM_80</span>","<span style=\"color: #5166f0\">CYBIM_130_1</span>",
             "<span style=\"color: #5166f0\">CYBIM_28</span>","<span style=\"color: #5166f0\">CYBIM_160</span>","<span style=\"color: #5166f0\">CYBIM_52</span>","<span style=\"color: #5166f0\">CYBIM_136_1</span>","<span style=\"color: #5166f0\">CYBIM_94</span>","<span style=\"color: #5166f0\">CYBIM_119</span>",
             "<span style=\"color: #5166f0\">CYBIM_73</span>","<span style=\"color: #5166f0\">CYBIM_76</span>","<span style=\"color: #5166f0\">CYBIM_15</span>","<span style=\"color: #5166f0\">CYBIM_77</span>","<span style=\"color: #5166f0\">CYBIM_22</span>","<span style=\"color: #5166f0\">CYBIM_90</span>",
             "<span style=\"color: #5166f0\">CYBIM_39</span>","<span style=\"color: #5166f0\">CYBIM_157</span>","<span style=\"color: #5166f0\">CYBIM_104</span>","<span style=\"color: #5166f0\">CYBIM_89</span>","<span style=\"color: #5166f0\">CYBIM_60_1</span>","<span style=\"color: #5166f0\">CYBIM_101</span>",
             "<span style=\"color: #5166f0\">CYBIM_93</span>","<span style=\"color: #5166f0\">CYBIM_63</span>","<span style=\"color: #5166f0\">CYBIM_73_1</span>","<span style=\"color: #5166f0\">CYBIM_190</span>")

dna_mech_plot_sum <- left_join(dna_mech_plot_sum, cols, by = c("genome" = "levels"))

dna_mech_plot_sum %>% ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  facet_wrap(~gene_facet, ncol = 15, drop = T)+
  theme(axis.text.x = element_text(size=7, angle =90, vjust = 0.5),
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))


dna_mech_plot_sum <- dna_mech_plot_sum %>% filter(count_samps10cov >= 20)
dna_mech_plot_sum$log_gene.count <- log10(dna_mech_plot_sum$sum_count)

#ordered_by_phylo y = factor(genome, levels = levels)
#ordered_by_nucl.div y = fct_reorder(genome,mean_instrain.above10cov_nucldiv,median)

dna_mech_plot_sum <- dna_mech_plot_sum %>% filter(genome != "VAMPV_156" |
                                                    genome != "VAMPV_261")
dna_mech_plot_sum <- dna_mech_plot_sum %>% filter(tax_abbr_name.col %in% levels_2)

p1 <- dna_mech_plot_sum %>% filter(gene_facet == "Alk") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Alk", y = "Cyano mOTUs", fill = "Gene Count") 
  
p2 <- dna_mech_plot_sum %>% filter(gene_facet == "DNA Uptake") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "DNA Uptake", y = "Cyano mOTUs", fill = "Gene Count")

p4 <- dna_mech_plot_sum %>% filter(gene_facet == "Mut") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Mut", y = "Cyano mOTUs", fill = "Gene Count") 

p6 <- dna_mech_plot_sum %>% filter(gene_facet == "RAD") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "RAD", y = "Cyano mOTUs", fill = "Gene Count")

p7 <- dna_mech_plot_sum %>% filter(gene_facet == "Rec") %>% 
  filter(gene_sh != "RepA") %>%
  filter(gene_sh != "RepA*") %>% 
  filter(gene_sh != "Cas4") %>% 
  filter(gene_sh != "Cas4*") %>% 
  subset(!grepl("COG", gene_sh)) %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Rec", y = "Cyano mOTUs", fill = "Gene Count") 

p8 <- dna_mech_plot_sum %>% filter(gene_facet == "Ruv") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Ruv", y = "Cyano mOTUs", fill = "Gene Count") 

interest <- c("PSEUDA_13", "PSEUDA_70", "DOLIS_187", "APHAN_134", "MCYST_162_1", "MCYST_62", "MCYST_79", "MCYST_100", 
              "MCYST_2", "MCYST_162", "MCYST_24", "MCYST_31", "MCYST_162_2", "MCYST_157", "MCYST_56", "CYANO_51_1")
p9 <- dna_mech_plot_sum %>%
  #filter(genome %in% interest) %>% 
  filter(gene_facet == "Transposases") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = log10(sum_count))) +
  #ggplot(aes(x = gene_sh, y = fct_reorder(genome,mean_instrain.above10cov_nucldiv,median), fill = log10(sum_count))) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  labs(title = "Transposases", y = "Cyano mOTUs", fill = "Gene Count")

p10 <- dna_mech_plot_sum %>% filter(gene_facet == "Uvr") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Uvr", y = "Cyano mOTUs", fill = "Gene Count")

p11 <- dna_mech_plot_sum %>% filter(gene_facet == "Xer") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Xer", y = "Cyano mOTUs", fill = "Gene Count")

p3 <- dna_mech_plot_sum %>% filter(gene_facet == "Endonuclease") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Endonucleases", y = "Cyano mOTUs", fill = "Gene Count") 


p5 <- dna_mech_plot_sum %>% filter(gene_facet == "Photolyase") %>% 
  ggplot(aes(x = gene_sh, y = factor(tax_abbr_name.col, levels = levels_2), fill = sum_count)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  labs(title = "Photolyase", y = "Cyano mOTUs", fill = "Gene Count") 



library(ggpubr)
ggarrange(p1, p2, p4, p6, p7, p8, p9, p10, p11, p3, p5,
          nrow = 1, align = "h")

plot <- ggarrange(p2 + theme.full,
          p4 + theme.skim,
          p7 + theme.skim,
          p8 + theme.skim,
          p10 + theme.skim,
          p11 + theme.skim,
          align = "h",
          nrow = 1)

annotate_figure(plot, top = text_grob("Cyano mOTUs orderd by Tree", 
                                      color = "black", face = "bold", size = 14))

# transposases 
p9 + theme.full
### Presentable figure

theme.full <-   theme_bw() +
  theme(axis.text.x = element_text(size=15, angle =90, vjust = 0.5),
        axis.text.y = ggtext::element_markdown(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "inches"))

theme.skim <-   theme_bw() +
  theme(axis.text.x = element_text(size=15, angle =90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "inches"))





#################################
### Using Anvio cog stuff (L)

d_anvio.genefxns <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene_functions.txt"
d_anvio.genecalls <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/gene-calls.txt"

anvio.genefxns <- read_delim(d_anvio.genefxns)
anvio.genecalls <- read_delim(d_anvio.genecalls)

anvio.genefxns <- left_join(anvio.genefxns, select(anvio.genecalls, c(1,2)), by = "gene_callers_id")
anvio.genefxns$genome <- gsub(".Contig.+", "", anvio.genefxns$contig)

anvio.genefxns_L_gene.calls <- anvio.genefxns[grepl("L", anvio.genefxns$accession),]

anvio.genefxns_L_gene.calls <- left_join(anvio.genefxns_L_gene.calls, filter(anvio.genefxns, source == "COG20_FUNCTION"), by = "gene_callers_id")

write_csv(anvio.genefxns_L_gene.calls, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/anvio_cogL.csv")

df_L_cogs_sh <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/anvio_cogL.csv"
L_cogs_sh <- read.delim(df_L_cogs_sh, sep = ",")

anvio.genefxns_L_gene.calls <- left_join(anvio.genefxns_L_gene.calls, L_cogs_sh, by = "function.y")

anvio.genefxns_L_gene.calls$present <- 1

anvio.genefxns_L_gene.calls <- anvio.genefxns_L_gene.calls %>% select(c(genome.x, gene_short, present))


anvio.genefxns_L_gene.calls_sum <- anvio.genefxns_L_gene.calls %>% 
  group_by(genome.x, gene_short) %>% 
  summarise(sum = sum(present))

levels = c("CYBIM_77","CYBIM_231","CYBIM_76","CYBIM_15","CYBIM_130","CYBIM_114","CYBIM_101","CYBIM_93",
           "CYBIM_63","CYBIM_73_1","CYBIM_190","CYBIM_96","CYBIM_101_1","CYBIM_60_1","CYBIM_89",
           "CYBIM_157","CYBIM_134_1","CYBIM_39","CYBIM_22","CYBIM_90","CYBIM_104",
           "CYBIM_31","CYBIM_51","CYBIM_200","CYBIM_80","CYBIM_130_1",
           "CYBIM_119","CYBIM_73","CYBIM_94","CYBIM_52","CYBIM_136_1","CYBIM_28","CYBIM_160","CYBIM_119_1","CYBIM_173",
           "VULCA_43","VULCA_242","VULCA_51","VULCA_120","VULCA_20","VULCA_96","VULCA_28","VULCA_209","VULCA_54", # VULCA_54 not in tree
           "CYBIM_67","CYANO_79",
           "NODOS_105","NODOS_42","NODOS_103","CYANO_45","CYANO_51_1",
           "MCYST_56","MCYST_157","MCYST_162_2","MCYST_31","MCYST_24","MCYST_162", #MCYST_162 not in tree
           "MCYST_2","MCYST_100","MCYST_79","MCYST_62","MCYST_162_1",
           "APHAN_134","DOLIS_187","DOLIS_36","DOLIS_53","DOLIS_105","CYANO_84","CYANO_106",
           "PSEUDA_31","PSEUDA_157","PSEUDA_281","PSEUDA_70","PSEUDA_13",
           "VAMPV_261","VAMPV_156")

anvio.genefxns_L_gene.calls_sum %>% ggplot(aes(x = gene_short, y = factor(genome.x, levels = levels), fill = sum)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(size=7, angle =90, vjust = 0.5),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))




# Run clustering
anvio.genefxns_L_gene.calls <- anvio.genefxns_L_gene.calls %>% distinct()
anvio.genefxns_L_gene.calls_wide <- pivot_wider(anvio.genefxns_L_gene.calls, names_from = "genome.x", values_from = "present")

anvio.genefxns_L_gene.calls_wide[is.na(anvio.genefxns_L_gene.calls_wide)] <- 0

matrix <- as.matrix(select(anvio.genefxns_L_gene.calls_wide, -c("gene_short")))
rownames(matrix) <- anvio.genefxns_L_gene.calls_wide$gene_short
matrix_dendro <- as.dendrogram(hclust(d = dist(x = matrix)))

# Create dendro
ggdendrogram(data = matrix_dendro, rotate = TRUE)+ theme(axis.text.y = element_text(size = 6))


anvio.genefxns_L_gene.calls <- pivot_longer(anvio.genefxns_L_gene.calls_wide, -c("gene_short"), 
                                            values_to = "present", names_to = "genome.x")


matrix_order <- order.dendrogram(matrix_dendro)

anvio.genefxns_L_gene.calls$gene_short <- factor(x = anvio.genefxns_L_gene.calls$gene_short,
                                                 levels = anvio.genefxns_L_gene.calls_wide$gene_short[matrix_order], 
                                                 ordered = TRUE)

anvio.genefxns_L_gene.calls <- left_join(anvio.genefxns_L_gene.calls, anvio.genefxns_L_gene.calls_sum, by = c("genome.x", "gene_short"))

anvio.genefxns_L_gene.calls %>% ggplot(aes(x = factor(genome.x, levels = levels), y = gene_short, fill = sum)) +
  geom_tile() +
  scale_fill_gradient2() +
  labs() +
  theme_bw() +
  theme(axis.text.x = element_text(size=7, angle =90, vjust = 0.5),
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=15, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=15, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"))













#############################
anvio.genefxns_L_gene.calls_sum <- anvio.genefxns_L_gene.calls_sum[!grepl("/", anvio.genefxns_L_gene.calls_sum$gene_short),]

library(pheatmap)

anvio.genefxns_L_gene.calls_sum <- anvio.genefxns_L_gene.calls_sum %>% pivot_wider(names_from = gene_short, values_from = sum)

anvio.genefxns_L_gene.calls_sum[is.na(anvio.genefxns_L_gene.calls_sum)] <- 0
anvio.genefxns_L_gene.calls_sum %>% select(where(~n_distinct(.) > 1))

#vapply(anvio.genefxns_L_gene.calls_sum, function(x) length(unique(x)) > 1, logical(1L))
#anvio.genefxns_L_gene.calls_sum[vapply(anvio.genefxns_L_gene.calls_sum, function(x) length(unique(x)) > 1, logical(1L))]


m <- as.matrix(select(anvio.genefxns_L_gene.calls_sum, -c("RecQ", "Ssb", "SSL2", "MutS2", "GyrA", "RecG", "TopA")))
pheatmap(m)

anvio.genefxns_L_gene.calls_simp <- anvio.genefxns_L_gene.calls %>% distinct()

anvio.genefxns_L_gene.calls_simp <- anvio.genefxns_L_gene.calls_simp %>% 
  pivot_wider(names_from = "gene_short.x", values_from = "present")

anvio.genefxns_L_gene.calls_simp[is.na(anvio.genefxns_L_gene.calls_simp)] <- 0
anvio.genefxns_L_gene.calls_simp <- anvio.genefxns_L_gene.calls_simp %>% select(where(~n_distinct(.) > 1))



m <- as.matrix(anvio.genefxns_L_gene.calls_simp)
pheatmap(m)


### ggtree
library(ape)
library(ggtree)

tree <- read.tree("/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/anvio/tree/16ribosomal-fasta_tree.txt")
tree
tree_root <- ape::root(tree, node = 128, resolve.root=TRUE)
tree_toot <- root(tree, 128)

ggplot(tree) + geom_tree() + geom_tiplab() + theme_tree()
ggtree(tree) + geom_text(aes(label=node), hjust=-.3) + geom_tiplab()


ptree <- ggtree(tree_root) + geom_tiplab() + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches")) +
  geom_text(aes(label=node), hjust=-.3)

flip(ptree, 147, 124) %>% flip(124)

ptree <- ggtree(tree_root) + geom_tiplab() + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches")) %>% gheatmap(anvio.genefxns_L_gene.calls_simp)



get_taxa_name(p)


################################
# Data
df_blast <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/cyano-diamond_dna-mech_blast.txt"
blast <- read.delim(df_blast)

df_ncbi_seq <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/blastn_seqIDs.txt"
ncbi_seq <- read.delim(df_ncbi_seq)

# NZ_LAHB01000001.1:131028-132905 uvrC [organism=Nostoc linckia z6] [GeneID=57091776] [chromosome=]
ncbi_seq <- separate(ncbi_seq, col = gene_seq, into = c("Seq_ID", "Gene_norm", "organism", "GeneID", "chr"), sep = "\t")


ncbi_seq$Seq_ID <- gsub(" .+", "", ncbi_seq$gene_seq)

ncbi_seq$Gene_norm <- gsub(" .organism.+", "", ncbi_seq$gene_seq)
ncbi_seq$Gene_norm <- gsub(".+ ", "", ncbi_seq$Gene_norm)

ncbi_seq$organism <- gsub(". .GeneID.+", "", ncbi_seq$gene_seq)
ncbi_seq$organism <- gsub(".+organism=", "", ncbi_seq$organism)

ncbi_seq$GeneID <- gsub(". .chrom.+", "", ncbi_seq$gene_seq)
ncbi_seq$GeneID <- gsub(".+GeneID=", "", ncbi_seq$GeneID)


### combine
blast <- left_join(blast, ncbi_seq, by = c("sseqid" = "Seq_ID"))

write_csv(blast, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/dna-mech-ncbi/blast_results.csv")


blast %>% group_by(qseqid, Gene_norm) %>% unique()

blast_set <- blast %>% 
  group_by(qseqid, Gene_norm) %>% 
  filter(n() == 1)

blast_set <- blast %>% 
  group_by(qseqid, Gene_norm) %>% 
  slice_max(order_by = pident)

blast_set$present <- 1
blast_set$genome <- gsub(".Contig.+", "", blast_set$qseqid)


blast_set %>% ggplot(aes(x = Gene_norm, y = genome, fill = present)) +
  geom_tile()


