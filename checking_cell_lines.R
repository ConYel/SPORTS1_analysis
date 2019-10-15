suppressPackageStartupMessages({
  library('tidyverse')
  library('data.table')
})

# introductory folder to load and save -------
todate <- format(Sys.time(), "%d_%b_%Y")
path_samples <- "R_virtual_shared_folder"
my_exp <- "mRNA_COLO205_pub"
dir.create(str_glue("/home/0/R_virtual_shared_folder/{my_exp}_analysis"))
dat_path <- str_glue("{my_exp}_analysis")

# load rds for hg38 feature counts -------
mRNA_feature <- dir(path_samples, full.names = TRUE,
            pattern = "mRNA_feature_counts_fastp_May_19_2019",
            recursive = TRUE)
dt2 <- read_rds(mRNA_feature)
hg38counts <- dt2$counts %>% 
  as_tibble(rownames = "smallRNA") %>% 
  rename(counts = 2) %>% 
  inner_join(dt2$annotation, by = c("smallRNA" = "GeneID")) 
# load mRNA expression HG19 -------
mRNA <- dir(path_samples, full.names = TRUE,
            pattern = "mRNA_raw_counts_COLO205_May_19_2019",
            recursive = TRUE)
dt1 <- read_tsv(mRNA)

mRNA <- dir(path_samples, full.names = TRUE,
                      pattern = "COLO205.txt",
                      recursive = TRUE)
dt <- fread(mRNA, select = 1:2) %>% 
  as_tibble() %>% 
  filter(!V1 %in% c("pRNA", "__no_feature",
                    "__ambiguous", "__too_low_aQual",
                    "__not_aligned", "__alignment_not_unique"),
         V2 > 0)

# load predicted targets IPP ------
targets_path <- dir(path_samples, full.names = TRUE,
                    pattern = "^all_IPPgenes_",
                    recursive = TRUE)

targets <- read_tsv(targets_path)
dt <- dt %>% mutate(logV2 = 
                log2(V2+1),
                cpm = V2 / sum(V2),
                lcpm = log2(cpm),
                color_expr = ifelse(V1 %in% targets$external_gene_name, "red", "black"))

DF <- dt %>% filter(V1 %in% targets$external_gene_name)

plot(dt$logV2,pch = 20,col = dt$color_expr)

# load Cell Lines file with expression --------
cell_lines_file <- dir(path_samples, full.names = TRUE,
                      pattern = "^all_mat_",
                      recursive = TRUE)
cpm_cell_lines <- read_tsv(cell_lines_file[1])
lcpm_cell_lines <- read_tsv(cell_lines_file[1])
cpm_cell_lines %>% separate(smallRNA, c("db","snRNA"), "_match_") %>% distinct(db)

cpm_cell_lines %>% filter(str_detect(smallRNA,"ensembl_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"miRNA_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"noAnnot_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"piRNA_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"rRNA_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"rfam_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"tRNA_match")) %>% summary
cpm_cell_lines %>% filter(str_detect(smallRNA,"tRNA_mature_match")) %>% summary

bds <- cpm_cell_lines %>% separate(smallRNA, c("db","snRNA"), "_match_")
bds %>% group_by(db) %>% 
  select(SW1417,RKO,Human_Testis) %>% 
  summarise_all(c(min = min, 
                      q25 = partial(quantile, probs = 0.25), 
                      median = median, 
                      q75 = partial(quantile, probs = 0.75), 
                      max = max,
                      mean = mean, 
                      sd = sd))
# load Yin et all piRNAs found in tissues -------

Yin_file <- dir(path_samples, full.names = TRUE,
                       pattern = "^data_of_Yin.xlsx",
                       recursive = TRUE)
Yin_table <- readxl::read_xlsx(Yin_file)

# load 317 IPP piRNAs ------
piRNA_IPP_317_file <- dir(path_samples, full.names = TRUE,
                pattern = "piRNA317_raw_reads.Input.txt",
                recursive = TRUE)

piRNA_IPP_317_file_table <- read_tsv(piRNA_IPP_317_file) %>% select(sRNA)
piRNA_IPP_317_file_table %>% inner_join(Yin_table, c("sRNA" = "piRBase"))
piRNA_IPP_317_file_table %>% filter(sRNA %in% Yin_table$piRBase)

# load IPP enriched smallRNAs ------
sRNA_IPP_all_file <- dir(path_samples, full.names = TRUE,
                          pattern = "LFCs_tables_comparisons",
                          recursive = TRUE)

sRNA_IPP_all_table <- read_tsv(sRNA_IPP_all_file) %>% 
  separate(smallRNA, c("db","snRNA"), "_match_")
  
seq_IPP <- DTRNA  %>% select(V2,V4,V6)

are_they_enriched <- sRNA_IPP_all_table %>% 
  inner_join(seq_IPP, by = c("snRNA" = "V4")) %>% 
  filter(V6 %in% Yin_table$sequence)

are_they_enriched %>% distinct(snRNA, .keep_all = TRUE)

# load methylation piRNAs -------
piRNA_methyl_file <- dir(path_samples, full.names = TRUE,
                          pattern = "Relative_ratio_Total_Apr_03_2019",
                          recursive = TRUE)[1]
                          
piRNA_methyl_file_table <- read_tsv(piRNA_methyl_file) 
piRNA_methyl_file_table <- read_tsv(piRNA_methyl_file) %>% separate(smallRNA, c("DB", "sRNA"), "_match_")

piRNA_methyl_file_table %>% group_by(DB,Status) %>% summarise(n())

methyl <- piRNA_methyl_file_table %>% 
  filter(Status %in% c("methylated","partially-methylated")) %>% 
  select(DB,sRNA, Tot_mean, Tot_NAIO4_mean)
methyl %>% filter(sRNA %in% Yin_table$piRBase)


piRNAs_com <- piRNA_methyl_file_table %>% 
  filter(sRNA %in% Yin_table$piRBase) %>% 
  select(sRNA,Tot_mean,Tot_median,Tot_NAIO4_mean,Tot_NAIO4_median,Ratio_mean,Ratio_median, Status)

# load the ambiguous piRNAs ------
piRNA_ambiguous_file <- dir(path_samples, full.names = TRUE,
                         pattern = "ambiguous",
                         recursive = TRUE)
piRNA_ambiguous_table <- readxl::read_xlsx(piRNA_ambiguous_file, col_names = TRUE) 
names(piRNA_ambiguous_table) <- names(piRNA_ambiguous_table) %>%
  str_replace_all(" ","_") 


# load the complete list of piRNAs in IPP ------
piRNA_IPP_file <- dir(path_samples, full.names = TRUE,
                            pattern = "IPP_journal",
                            recursive = TRUE)
piRNA_IPP <- read_tsv(piRNA_IPP_file)

# load the list of the 103 genes found in IPP -------
piRNA_genes_IPP_file <- dir(path_samples, full.names = TRUE,
                      pattern = "PIWIL1_piRNAs complex IP_May_2019",
                      recursive = TRUE)
piRNA_genes_IPP <- readxl::read_xlsx(piRNA_genes_IPP_file)
names(piRNA_genes_IPP) <- names(piRNA_genes_IPP) %>%
  str_replace_all(" ","_") 

piRNA_genes_IPP <- piRNA_genes_IPP %>% 
  select(List_5.2) %>% filter(!is.na(List_5.2))

dt %>% filter(V1 %in%  piRNA_genes_IPP$List_5.2)

dt <- dt %>% mutate(logV2 = log2(V2+1),
                    cpm = V2 / sum(V2),
                    lcpm = log2(cpm),
              color_expr = ifelse(V1 %in% piRNA_genes_IPP$List_5.2, "red", "black"))

dt %>%  select(logV2) %>% 
  summarise_all(c(min = min, 
                  q25 = partial(quantile, probs = 0.25), 
                  median = median, 
                  q75 = partial(quantile, probs = 0.75), 
                  max = max,
                  mean = mean, 
                  sd = sd))

dt %>% filter(V1 %in%  piRNA_genes_IPP$List_5.2) %>%  
  select(logV2) %>% 
  summarise_all(c(min = min, 
                  q25 = partial(quantile, probs = 0.25), 
                  median = median, 
                  q75 = partial(quantile, probs = 0.75), 
                  max = max,
                  mean = mean, 
                  sd = sd))

# cross all of them with the ambiguous ----
piRNA_methyl_file_table %>% filter(sRNA %in% piRNA_ambiguous_table$piRBase_ID) %>% select(sRNA)
piRNA_IPP_317_file_table %>% filter(sRNA %in% piRNA_ambiguous_table$piRBase_ID) 
Yin_table%>% filter(piRBase %in% piRNA_ambiguous_table$piRBase_ID)
bds %>% filter(snRNA %in% piRNA_ambiguous_table$piRBase_ID)
piRNA_IPP %>% filter(smallRNA %in% piRNA_ambiguous_table$piRBase_ID) %>% view

# make plots boxplot ------
dt %>% 
  ggplot(aes(x = "expression", y = lcpm)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), col = dt$color_expr)

# find which DE smallRNA are piRNAs in the study -------
LFCs_C_N_IPP_file <- dir(path_samples, full.names = TRUE,
                            pattern = "LFCs_C_N_2_input_2_results.txt",
                            recursive = TRUE)
LFCs_C_N <- read_tsv(LFCs_C_N_IPP_file) %>% 
  separate(smallRNA,c("DB","sRNA"),"_match_")
LFCs_C_N %>% group_by(DB) %>% summarise(n())
seq_of_573_DE_IPP <- DTRNA %>% filter(V4 %in% LFCs_C_N$sRNA) %>% select(FileName,V2,V4,V6)
LFCs_C_N %>% 
  inner_join(seq_of_573_DE_IPP, c("sRNA" = "V4")) %>% 
  arrange(desc(V2)) %>%  
  distinct(V6, .keep_all = TRUE) %>% 
  inner_join(Yin_table, c("V6" = "sequence")) %>% write_tsv("enriched_smallRNAs_yin.txt")
