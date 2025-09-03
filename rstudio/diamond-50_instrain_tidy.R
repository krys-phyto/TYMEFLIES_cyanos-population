######################################################
### Krys Kibler 2025-05-27
### Tidying instrain output files
######################################################

library(tidyverse)
library(data.table)

######################################################
### Taxonomy
cyano.mOTUs <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/data/files/cyanos_genome_summary.csv"
cyano.mOTUs <- read_delim(cyano.mOTUs)

# create genera column
cyano.mOTUs$genera <- gsub("d__Bacteria;.+g__", "", cyano.mOTUs$taxonomy)
cyano.mOTUs$genera <- paste("g__", cyano.mOTUs$genera, sep = "")
cyano.mOTUs$genera <- gsub(";s.+", "", cyano.mOTUs$genera)
cyano.mOTUs$genera <- gsub("g__", "", cyano.mOTUs$genera)

### TYMEFLIES Metadata
metadata <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/data/files/all_JGI_metadata_471_samples.txt"
metadata <- read_delim(metadata)

# obtain a sampledata column
metadata <- metadata %>% select(SampleDate, SampleID)
metadata$SampleDate <- substr(metadata$SampleDate, 3, 11)
metadata$SampleDate <- dmy(metadata$SampleDate)
metadata <- metadata %>% na.omit()

### Instrain
#instrain_95_mapping <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-95/cyano_50_readANI95_mapping.tsv"
#instrain_92_mapping <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-92/cyano_50_readANI92_mapping.tsv"

instrain_95_genome <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-96/cyano_50_readANI95_genome.tsv"
instrain_92_genome <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-92/cyano_50_readANI92_genome.tsv"

#instrain_95_linkage <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-96/cyano_50_readANI95_linkage_trimmed.tsv"
instrain_92_linkage <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/analysis/instrain/50_cyano-motu_instrain-92/cyano_50_readANI92_linkage_trimmed.tsv"

#instrain_95_mapping <- read_delim(instrain_95_mapping)
#instrain_92_mapping <- read_delim(instrain_92_mapping)

instrain_95_genome <- fread(instrain_95_genome)
instrain_92_genome <- fread(instrain_92_genome)

#instrain_95_linkage <- fread(instrain_95_linkage)
instrain_92_linkage <- fread(instrain_92_linkage)


######################################################


######################################################
### Tidying

#False header line
#instrain_92_mapping <- instrain_92_mapping[-1, ]
#instrain_95_mapping <- instrain_95_mapping[-1, ]

instrain_92_genome <- instrain_92_genome[-1, ]
instrain_95_genome <- instrain_95_genome[-1, ]

colnames(instrain_92_linkage) <- c("scaffold", "distance", "r2", "d_prime", "r2_normalized", "d_prime_normalized", "total", "sampleID")

# Merge the sampledate
#instrain_92_mapping <- left_join(instrain_92_mapping, metadata, by = c("sampleID" = "SampleID"))
#instrain_95_mapping <- left_join(instrain_95_mapping, metadata, by = c("sampleID" = "SampleID"))

instrain_92_genome <- left_join(instrain_92_genome, metadata, by = c("sampleID" = "SampleID"))
instrain_95_genome <- left_join(instrain_95_genome, metadata, by = c("sampleID" = "SampleID"))

instrain_92_linkage <- left_join(instrain_92_linkage, metadata, by = c("sampleID" = "SampleID"))

# Fix the genome column
instrain_92_genome$genome <- gsub(".fna", "", instrain_92_genome$genome)
instrain_95_genome$genome <- gsub(".fna", "", instrain_95_genome$genome)

instrain_92_linkage$genome <- gsub(".Contig_.*", "", instrain_92_linkage$scaffold)

#instrain_92_mapping$genome <- gsub(".Contig_.+", "", instrain_92_mapping$scaffold)
#instrain_95_mapping$genome <- gsub(".Contig_.+", "", instrain_95_mapping$scaffold)

# Add readANI value to column names
colnames_92 <- colnames(instrain_92_genome)
colnames_92 <- paste(colnames_92,"_92ANI",sep = "")

colnames_95 <- colnames(instrain_95_genome)
colnames_95 <- paste(colnames_95,"_95ANI",sep = "")

colnames(instrain_92_genome) <- c("genome","coverage_92ANI","breadth_92ANI","nucl_diversity_92ANI","length","true_scaffolds_92ANI","detected_scaffolds_92ANI","coverage_median_92ANI","coverage_std_92ANI","coverage_SEM_92ANI","breadth_minCov_92ANI","breadth_expected_92ANI","nucl_diversity_rarefied_92ANI","conANI_reference_92ANI","popANI_reference_92ANI","iRep_92ANI","iRep_GC_corrected_92ANI","linked_SNV_count_92ANI","SNV_distance_mean_92ANI","r2_mean_92ANI","d_prime_mean_92ANI","consensus_divergent_sites_92ANI","population_divergent_sites_92ANI","SNS_count_92ANI","SNV_count_92ANI","filtered_read_pair_count_92ANI","reads_unfiltered_pairs_92ANI","reads_mean_PID_92ANI","divergent_site_count_92ANI","reads_unfiltered_reads_92ANI","sampleID","SampleDate")    
colnames(instrain_95_genome) <- c("genome","coverage_95ANI","breadth_95ANI","nucl_diversity_95ANI","length","true_scaffolds_95ANI","detected_scaffolds_95ANI","coverage_median_95ANI","coverage_std_95ANI","coverage_SEM_95ANI","breadth_minCov_95ANI","breadth_expected_95ANI","nucl_diversity_rarefied_95ANI","conANI_reference_95ANI","popANI_reference_95ANI","iRep_95ANI","iRep_GC_corrected_95ANI","linked_SNV_count_95ANI","SNV_distance_mean_95ANI","r2_mean_95ANI","d_prime_mean_95ANI","consensus_divergent_sites_95ANI","population_divergent_sites_95ANI","SNS_count_95ANI","SNV_count_95ANI","filtered_read_pair_count_95ANI","reads_unfiltered_pairs_95ANI","reads_mean_PID_95ANI","divergent_site_count_95ANI","reads_unfiltered_reads_95ANI","sampleID","SampleDate")    

# add the two instrain genome dataframes
instrain_genome <- left_join(instrain_92_genome, instrain_95_genome, by = c("genome","length","sampleID","SampleDate"))
instrain_genome$count <- 1


# Full dataframe for instrain genome level diversity
instrain_genome <- left_join(cyano.mOTUs, instrain_genome, by = c("tax_abbr_name" = "genome"))
instrain_genome <- instrain_genome %>% filter(morphology != "parasite")

# change character columns into numeric
cols.num <-  c("coverage_92ANI","breadth_92ANI","nucl_diversity_92ANI","length","true_scaffolds_92ANI","detected_scaffolds_92ANI","coverage_median_92ANI","coverage_std_92ANI","coverage_SEM_92ANI","breadth_minCov_92ANI","breadth_expected_92ANI","nucl_diversity_rarefied_92ANI","conANI_reference_92ANI","popANI_reference_92ANI","iRep_92ANI","linked_SNV_count_92ANI","SNV_distance_mean_92ANI","r2_mean_92ANI","d_prime_mean_92ANI","consensus_divergent_sites_92ANI","population_divergent_sites_92ANI","SNS_count_92ANI","SNV_count_92ANI","filtered_read_pair_count_92ANI","reads_unfiltered_pairs_92ANI","reads_mean_PID_92ANI","divergent_site_count_92ANI","reads_unfiltered_reads_92ANI",
               "coverage_95ANI","breadth_95ANI","nucl_diversity_95ANI","true_scaffolds_95ANI","detected_scaffolds_95ANI","coverage_median_95ANI","coverage_std_95ANI","coverage_SEM_95ANI","breadth_minCov_95ANI","breadth_expected_95ANI","nucl_diversity_rarefied_95ANI","conANI_reference_95ANI","popANI_reference_95ANI","iRep_95ANI","linked_SNV_count_95ANI","SNV_distance_mean_95ANI","r2_mean_95ANI","d_prime_mean_95ANI","consensus_divergent_sites_95ANI","population_divergent_sites_95ANI","SNS_count_95ANI","SNV_count_95ANI","filtered_read_pair_count_95ANI","reads_unfiltered_pairs_95ANI","reads_mean_PID_95ANI","divergent_site_count_95ANI","reads_unfiltered_reads_95ANI")
instrain_genome[cols.num] <- sapply(instrain_genome[cols.num],as.numeric)

# Full dataframe for instrain linkage
instrain_92_linkage <- left_join(cyano.mOTUs, instrain_92_linkage, by = c("tax_abbr_name" = "genome"))
instrain_92_linkage <- instrain_92_linkage %>% filter(morphology != "parasite")



######################################################
### Normalized Coverage

# total read pairs from cyano in samples
metag_mapped_92 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_reads_92 = sum(filtered_read_pair_count_92ANI))
metag_mapped_95 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_reads_95 = sum(filtered_read_pair_count_95ANI))

mean_mapped_92 = mean(metag_mapped_92$sum_reads_92) # avg mapped number of read pairs = 3464910
mean_mapped_95 = mean(metag_mapped_95$sum_reads_95, na.rm = TRUE) # avg mapped number of read pairs = 3097187

instrain_genome <- left_join(instrain_genome, metag_mapped_92, by = "SampleDate")
instrain_genome <- left_join(instrain_genome, metag_mapped_95, by = "SampleDate")

# normalized coverage calc
instrain_genome$norm_coverage_92 <- instrain_genome %>% with(coverage_92ANI * (mean_mapped_92 / sum_reads_92))
instrain_genome$norm_coverage_95 <- instrain_genome %>% with(coverage_95ANI * (mean_mapped_95 / sum_reads_95))

instrain_genome$norm_ratio_92 <- instrain_genome %>% with(mean_mapped_92 / sum_reads_92)
instrain_genome$norm_ratio_95 <- instrain_genome %>% with(mean_mapped_95 / sum_reads_95)

# relative abundance calc by coverage
sum_norm_coverage_92 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_norm.cov_92 = sum(norm_coverage_92))
sum_norm_coverage_95 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_norm.cov_95 = sum(norm_coverage_95))

instrain_genome <- left_join(instrain_genome, sum_norm_coverage_92, by = "SampleDate")
instrain_genome <- left_join(instrain_genome, sum_norm_coverage_95, by = "SampleDate")

instrain_genome$RA_by.cov_92 <- instrain_genome %>% with((norm_coverage_92 / sum_norm.cov_92) * 100)
instrain_genome$RA_by.cov_95 <- instrain_genome %>% with((norm_coverage_95 / sum_norm.cov_95) * 100)

sum_coverage_92 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_cov.92 = sum(coverage_92ANI))
sum_coverage_95 <- instrain_genome %>% group_by(SampleDate) %>% 
  summarise(sum_cov.95 = sum(coverage_95ANI))

instrain_genome <- left_join(instrain_genome, sum_coverage_92, by = "SampleDate")
instrain_genome <- left_join(instrain_genome, sum_coverage_95, by = "SampleDate")

rm(instrain_92_genome, instrain_92_mapping, instrain_95_genome, instrain_95_mapping, metadata, 
   metag_mapped_92, metag_mapped_95, sum_norm_coverage_92, sum_norm_coverage_95, sum_coverage_92, sum_coverage_95)
rm(cols.num, colnames_92, colnames_95, mean_mapped_92, mean_mapped_95)


genome <- c("PSEUDA_13","PSEUDA_70","PSEUDA_281","PSEUDA_157","PSEUDA_31",
            "NODOS_103","NODOS_105","CYANO_45","CYANO_51_1","CYANO_106",
            "CYANO_84","APHAN_134","DOLIS_187","DOLIS_105","MCYST_56",
            "MCYST_62","MCYST_31","MCYST_2","CYBIM_15","CYBIM_76",
            "CYBIM_77","CYBIM_22","CYBIM_90","CYBIM_39","CYBIM_157",
            "CYBIM_104","CYBIM_89","CYBIM_60_1","CYBIM_101","CYBIM_93",
            "CYBIM_63","CYBIM_190","CYBIM_73_1","CYBIM_200","CYBIM_31",
            "CYBIM_51","CYBIM_119_1","CYBIM_80","CYBIM_130_1","CYBIM_160",
            "CYBIM_28","CYBIM_52","CYBIM_136_1","CYBIM_94","CYBIM_119",
            "CYBIM_73","VULCA_20","VULCA_120","VULCA_28","VULCA_96")

col <- c("#4cd505","#4cd505","#4cd505","#4cd505","#4cd505",
            "#3cb306","#3cb306","#236a00","#236a00","#236a00",
            "#236a00","#236a00","#236a00","#236a00","#f4559d",
            "#e39f0c","#e39f0c","#e39f0c","#00dede","#00dede",
            "#00dede","#00dede","#00dede","#00dede","#00dede",
            "#00dede","#00dede","#00dede","#00dede","#00dede",
            "#00dede","#00dede","#00dede","#008bff","#008bff",
            "#008bff","#008bff","#008bff","#008bff","#008bff",
            "#008bff","#008bff","#008bff","#008bff","#008bff",
            "#008bff","#7a05c5","#7a05c5","#7a05c5","#7a05c5")


cols <- as.data.frame(cbind(genome,col))
cyano.mOTUs <- left_join(cyano.mOTUs, cols, by = c("tax_abbr_name" = "genome"))
instrain_genome <- left_join(instrain_genome, cols, by = c("tax_abbr_name" = "genome"))
instrain_92_linkage <- left_join(instrain_92_linkage, cols, by = c("tax_abbr_name" = "genome"))

cyano.mOTUs$tax_abbr_name.col <- paste0("<span style=\"color: ", cyano.mOTUs$col, "\">", cyano.mOTUs$tax_abbr_name, "</span>")
instrain_genome$tax_abbr_name.col<- paste0("<span style=\"color: ", instrain_genome$col, "\">", instrain_genome$tax_abbr_name, "</span>")
write_csv(instrain_genome, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/cyanoSTRONG/krys/rstudio/generated_dfs/cyano-mOTU_instrain-genome.csv")

instrain_92_linkage$tax_abbr_name.col <- paste0("<span style=\"color: ", instrain_92_linkage$col, "\">", instrain_92_linkage$tax_abbr_name, "</span>")

levels <- paste0("<span style=\"color: ", col, "\">", genome, "</span>")

rm(col,cols)
instrain_genome

# instrain_genome pop vs con ani 92% ANI table

dataframe <- instrain_genome %>% select(SampleDate, tax_abbr_name, tax_abbr_name.col, popANI_reference_92ANI, conANI_reference_92ANI)

dataframe$count <- 1
dataframe <- dataframe %>% group_by(tax_abbr_name, tax_abbr_name.col) %>% 
  summarise(mean_pop = mean(popANI_reference_92ANI),
            sd_pop = sd(popANI_reference_92ANI),
            mean_con = mean(conANI_reference_92ANI),
            sd_con = sd(conANI_reference_92ANI),
            count = sum(count))



### Find the sampleIDs for 2015 and bloom samples
dataframe <- instrain_genome %>% select(sampleID, SampleDate, tax_abbr_name, coverage_92ANI, norm_coverage_92) %>% filter(tax_abbr_name == "APHAN_134" | 
                                                                                                                            tax_abbr_name == "MCYST_31" |
                                                                                                                            tax_abbr_name == "MCYST_2")
dataframe$year4 <- year(dataframe$SampleDate)
dataframe$month <- month(dataframe$SampleDate)

samples2015 <- dataframe %>% filter(year4 == 2015)

samplesbloom <- dataframe %>%
  filter(month >= 5 &
           month <= 8) %>% 
  filter(year4 != 2015) %>% 
  group_by(year4, tax_abbr_name) %>%
  filter(coverage_92ANI == max(coverage_92ANI, na.rm = TRUE))

dataframe <- rbind(samples2015, samplesbloom)

dataframe1 <- dataframe %>% filter(tax_abbr_name == "APHAN_134")
unique(dataframe1$sampleID)

dataframe1 <- dataframe %>% filter(tax_abbr_name == "MCYST_31")
unique(dataframe1$sampleID)

dataframe1 <- dataframe %>% filter(tax_abbr_name == "MCYST_2")
unique(dataframe1$sampleID)

