library(tidyverse)
library(ggrepel)
library(phyloseq)
library(ANCOMBC)
library(qiime2R)  
library(limma)
library(Glimma)
library(edgeR)
library(ComplexHeatmap)
library(RColorBrewer)
library(vegan)
library(pairwiseAdonis)
library(rstatix)

# Place to save output results
dir.create(file.path(".", "sessioninfo"), recursive = TRUE, showWarnings = FALSE)
outdir <- file.path(".", "output_abundance")
dir.create(file.path(outdir, "toptables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "toptreat"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "MA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "mds"), recursive = TRUE, showWarnings = FALSE)

# Note that I started with run 2 (T282) so run 2 is note named and run 1 (T275) is named run -1

# qiime2 docs:
# https://amplicon-docs.readthedocs.io/en/latest/

# Convert metadata to tsv and make it match fastq file naming
md <- read_csv(file.path(".", "Metadata", "Metadata_overview_T282.csv"))%>%
  mutate(`sample-id` = gsub(" ","", `sample-id`))%>%
  dplyr::select(-`barcode-sequence`)

md_updated <- md[2:nrow(md),] %>%
  rename(light = "Constant light or typical light cycle")%>%
  rename(Days_since_infection = "days-since-experiment-start (where day 0 start means infection day)")%>%
  rename(Infection_status = "Infected or uninfected")%>%
  rename(Antibiotic = "reported-antibiotic-usage")%>%
  mutate(light = gsub(" ", "", light))%>%
  mutate(light = replace(light, light == "Baselinebeforeconstantlight", "Baseline"))%>%
  mutate(light = replace(light, light == "Baselinebeforetypicallight", "Baseline"))%>%
  mutate(all_group = paste0(Infection_status, "_", light, "_", Days_since_infection, "_",Antibiotic))%>%
  mutate(subject = gsub("-", "_", subject))%>%
  mutate(Run = "Run2")

# Add back on qiime2 classes
md_updated_save <- rbind(c(as.character(md[1,]), "categorical", "categorical"), md_updated)%>%
  write_tsv(file.path(".", "Qiime2", "sample-metadata.tsv"))

# Read and fix the metadata from the first run as well
run_1_md <-  read_tsv(file.path(".", "Metadata", "Metadata_overview_T275.csv"))%>%
  mutate(`sample-id` = gsub(" ","", `sample-id`))%>%
  # Looks like the barcode sequnce is the same as run 2 as well
  dplyr::select(-`barcode-sequence`)

run_1_md_updated <- run_1_md[2:nrow(run_1_md),] %>%
  rename(light = "Constant light or typical light cycle")%>%
  rename(Days_since_infection = "days-since-experiment-start (where start means infection day)")%>%
  rename(Infection_status = "Infected or uninfected")%>%
  rename(Antibiotic = "reported-antibiotic-usage")%>%
  mutate(light = gsub(" ", "", light))%>%
  mutate(light = replace(light, light == "Baselinebeforeconstantlight", "Baseline"))%>%
  mutate(light = replace(light, light == "Baselinebeforetypicallight", "Baseline"))%>%
  mutate(all_group = paste0(Infection_status, "_", light, "_", Days_since_infection, "_",Antibiotic))%>%
  # Fix baseline status naming consistency
  mutate(all_group = replace(all_group, all_group == c("Uninfected_Typicallight_-14_No"), "Uninfected_Baseline_neg_14_No"))%>%
  mutate(subject = gsub("-", "_", subject))%>%
  mutate(Run = "Run1")

# Add on qiime2 classes
run_1_md_updated_save <- rbind(c(as.character(run_1_md[1,]), "categorical", "categorical"), run_1_md_updated)%>%
  write_tsv(file.path(".", "Qiime2", "sample-metadata-run-1.tsv"))

# Save merged metadata from qiime2
merged_md_save <- rbind(c(as.character(run_1_md[1,]), "categorical", "categorical"), run_1_md_updated,md_updated)%>%
  write_tsv(file.path(".", "Qiime2", "sample-metadata-merged.tsv"))

# Make a manifest file will only need to be redone for if rerunning qiime2
# files <- list.files(file.path(".", "raw"), full.names = TRUE, pattern = "*.fastq.gz")
# 
# files_df <- data.frame(`absolute-filepath` = files,check.names = F)%>%
#   mutate(`sample-id` = basename(`absolute-filepath`))%>%
#   mutate(direction = ifelse(grepl("1.fastq.gz", `sample-id`), "forward", "reverse"))%>%
#   mutate(`sample-id` = gsub("_1.fastq.gz|_2.fastq.gz|", "", `sample-id`))%>%
#   dplyr::select(`sample-id`, `absolute-filepath`, direction)%>%
#   write_csv(file.path(".", "Qiime2", "manifest.csv"))

# Qiime2 output analysis starts here ----

# Read in the merged feature table
table_qza <- read_qza(file.path(".", "Qiime2", "merged_table.qza"))

taxonomy <- read_qza(file.path(".", "Qiime2", "merged-taxonomy.qza"))$data %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";",remove = F)

dt <- table_qza$data
mat <- as.matrix(dt)

anno_df <- data.frame(Feature.ID = rownames(mat))%>%
  left_join(taxonomy, by = "Feature.ID")

# Make an annotation object for limma using the qiime2 metadata
bulk_anno <- bind_rows(run_1_md_updated, md_updated)%>%
  mutate(all_group = gsub("-","neg_", all_group))%>%
  mutate(subject = gsub("-", "_", subject))%>%
  rename(sample.id = `sample-id`)%>%
  # Make sure subject is unique to the run
  mutate(subject = paste0(Run,"_", subject))%>%
  mutate(all_group = gsub("Uninfected", "Uni", all_group))%>%
  mutate(all_group = gsub("Infected", "Inf", all_group))%>%
  mutate(all_group = gsub("_neg_14_", "_n14_", all_group))%>%
  mutate(all_group = gsub("Constantlight", "Const", all_group))%>%
  mutate(all_group = gsub("Typicallight", "Typic", all_group))%>%
  mutate(all_group = gsub("_Yes", "_Anti", all_group))%>%
  mutate(all_group = gsub("_No", "_No_anti", all_group))%>%
  mutate(all_group = gsub("_Baseline", "", all_group))%>%
  mutate(group_run = paste0(all_group, "_", Run))%>%
  # Drop one sample that did not show cdiff infection T275 E5
  filter(sample.id != "T275E5_")

# Get a sensible factor order for the groups
order <- c("Uni_n14_No_anti", "Uni_Const_0_No_anti", "Uni_Typic_0_Anti", 
           "Uni_Const_0_Anti", "Uni_Typic_1.5_Anti", "Uni_Const_1.5_Anti", 
           "Inf_Typic_1.5_Anti", "Inf_Const_1.5_Anti", "Uni_Const_2_Anti",
           "Inf_Const_2_Anti", "Uni_Const_4_Anti", "Inf_Const_4_Anti")

# Check all the groups are there
order %in% bulk_anno$all_group

bulk_anno$all_group <- factor(bulk_anno$all_group, levels = order)

write_csv(bulk_anno, file.path(".", "Metadata", "Combined_sample_metadata_neat.csv"))

# Check everything lines up
sum(anno_df$Feature.ID == rownames(mat)) == nrow(mat)

mat <- mat[,bulk_anno$sample.id]
sum(colnames(mat) == bulk_anno$sample.id) == nrow(bulk_anno)

counts <- DGEList(mat, genes = anno_df[,2:ncol(taxonomy)])

# Individual mice are repeated
design <- model.matrix(~0 + all_group + Run, data = bulk_anno)
colnames(design) <- gsub("all_group", "", colnames(design))
rownames(design) <- rownames(counts$samples)

hist(rowSums(mat), breaks = 100000, xlim = c(0,100))

table(bulk_anno$all_group)

# You could use these counts for a limma-style analysis
#keep <- rowSums(mat>5)>3
keep <- filterByExpr(mat, design)
table(keep)
counts <- counts[keep,, keep.lib.sizes=FALSE]
dim(counts)

# Apply TMM normalisation to the DGElist
# TMMwsp works better when counts have a high proportion of zeros
counts <- calcNormFactors(counts, method = "TMMwsp")

lcpm <- edgeR::cpm(counts, log = T, normalized.lib.sizes = TRUE)

lcpm_df <- lcpm%>%
  data.frame()%>%
  rownames_to_column("Bacterium")

write_csv(lcpm_df, file.path(outdir, "TMM_normalised_log2_CPMs.csv"))

# A LOT of species are lost after antibiotic treatment
# bulk_anno_test <- bulk_anno
# bulk_anno_test$f936920a3a442ae373b5d3832c0f5a35 <- mat["f936920a3a442ae373b5d3832c0f5a35",]
# bulk_anno_test$lcpm_f936920a3a442ae373b5d3832c0f5a35 <- lcpm["f936920a3a442ae373b5d3832c0f5a35",]
mds_save <- file.path(outdir, "glimma", "mds", "MDS.html")

htmlwidgets::saveWidget(glimmaMDS(counts, groups = bulk_anno$all_group,
                                  labels = bulk_anno), mds_save)

colnames(design)

# Make a contrast matrix
cont.matrix <- makeContrasts(
  # Determine if constant light alters the microbiome at day 0
  Uni_Const_0_Anti_vs_Uni_Typic_0_Anti = Uni_Const_0_Anti - Uni_Typic_0_Anti,
  # Determine if constant light alters the microbiome at day 1.5
  Uni_Const_1.5_Anti_vs_Uni_Typic_1.5_Anti = Uni_Const_1.5_Anti - Uni_Typic_1.5_Anti,
  # Determine if constant light alters the microbiome at typical peak C. difficile infection
  Inf_Const_1.5_Anti_vs_Inf_Typic_1.5_Anti = Inf_Const_1.5_Anti - Inf_Typic_1.5_Anti,
  # Light average over both days
  light_both_days = (Uni_Const_0_Anti + Inf_Const_1.5_Anti)/2 - 
    (Uni_Typic_0_Anti + Inf_Typic_1.5_Anti)/2,
  # The effect of antibiotics
  Uni_Const_0_Anti_vs_Uni_Const_0_No_anti = Uni_Const_0_Anti - Uni_Const_0_No_anti,
  # The effect of time and light (but not antis) from d-14 to d0 
  Uni_Const_0_No_anti_vs_Uni_n14_No_anti = Uni_Const_0_No_anti - Uni_n14_No_anti,
  # The effect of diff
  Inf_Typic_1.5_Anti_vs_Uni_Typic_1.5_Anti = Inf_Typic_1.5_Anti - Uni_Typic_1.5_Anti,
  levels = colnames(design))

# Check my contrasts are equal to 0
colSums(cont.matrix)

fit <- voomLmFit(counts, design, plot=TRUE, block = bulk_anno$subject)
fit <- contrasts.fit(fit, contrasts=cont.matrix)

fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
summa.fit <- decideTests(fit)
summary(summa.fit)

# Make a TREAT fit object
tfit <- treat(fit, fc=1.5)

for(contrast in colnames(cont.matrix)){
  
  toptable <- topTable(fit,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(file.path(outdir, "toptables", paste0(contrast, "_toptable.csv")))
  
  toptreat <- topTreat(tfit,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(file.path(outdir, "toptreat", paste0(contrast, "_toptreat.csv")))
  
  vol_save <- file.path(outdir, "glimma", "volcano", paste0(contrast, "_Volcano.html"))
  
  htmlwidgets::saveWidget(glimmaVolcano(fit, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = lcpm,transform.counts	= "none",
                                        dge = counts, groups = bulk_anno$all_group), vol_save)
  
  ma_save <- file.path(outdir, "glimma", "MA", paste0(contrast, "_MA.html"))
  
  htmlwidgets::saveWidget(glimmaMA(fit, coef = contrast,main = gsub("_"," ",contrast),
                                   counts = lcpm,
                                   dge = counts, groups = bulk_anno$all_group), ma_save)
  
}

# Remove glimma intermediate files safely
files_to_remove <- list.files(file.path(outdir, "glimma"), pattern = "_files$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
unlink(files_to_remove, recursive = TRUE)

# Compile the toptables 
all_toptables <- list.files(file.path(outdir, "toptables"), full.names = TRUE)

tt_list <- list()
for(i in seq_along(all_toptables)){
  
  contrast <- gsub("_toptable\\.csv$", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
}

# Compile toptables and save the significant results
if(length(tt_list) > 0) {
  toptables_compiled <- bind_rows(tt_list)
  
  toptables_signif <- toptables_compiled %>%
    filter(adj.P.Val < 0.05)%>%
    arrange(adj.P.Val)%>%
    write_csv(file.path(outdir, "Compiled_toptables_significant_bacteria.csv"))
}

# Compile the toptreat 
all_toptreat <- list.files(file.path(outdir, "toptreat"), full.names = TRUE)

tt_list <- list()
for(i in seq_along(all_toptreat)){
  
  contrast <- gsub("_toptreat\\.csv$", "", basename(all_toptreat[i]))
  
  tt <- read_csv(all_toptreat[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
}

# Compile toptables and save the significant results
if(length(tt_list) > 0) {
  toptreat_compiled <- bind_rows(tt_list)
  
  toptreat_signif <- toptreat_compiled %>%
    filter(adj.P.Val < 0.05)%>%
    arrange(adj.P.Val)%>%
    write_csv(file.path(outdir, "Compiled_toptreat_significant_bacteria.csv"))
}

# Remove the run level batch effect from the data
# design2 <- model.matrix(~0 + all_group, data = bulk_anno)
# lcpm_batch_removed <- removeBatchEffect(lcpm, batch = bulk_anno$Run, design = design2)

# Generate the 12-color "Paired" palette
palette_paired <- brewer.pal(12, "Paired")

palette_paired[11] <- "pink"

# Create a named vector mapping group names to colors
cb_palette_final <- setNames(palette_paired, unique(bulk_anno$all_group))

newpal <- c("Uni_Typic_0_Anti" = "#A2DF8D",
            "Uni_Const_0_Anti" = "#D9BBF7", 
            "Uni_Typic_1.5_Anti" = "#008000",
            "Uni_Const_1.5_Anti" = "#C000C0",
            "Inf_Typic_1.5_Anti" = "#0396FE",
            "Inf_Const_1.5_Anti" = "#FFA040"
)

cb_palette_final[names(newpal)] <- newpal

# Run Bray Shannon indexes for both runs but remove the groups without antis
for(run in unique(bulk_anno$Run)){
  
  bulk_run <- bulk_anno%>%
    filter(Run %in% run)%>%
    filter(Antibiotic == "Yes")
  
  # Create phyloseq object
  otu_table <- otu_table(mat[,bulk_run$sample.id], taxa_are_rows = TRUE)
  
  # Make a sample data object
  sample_data <- sample_data(bulk_run %>% column_to_rownames("sample.id"))
  
  # Create a phyloseq object WITH RAW counts
  physeq_raw <- phyloseq(otu_table(mat[,bulk_run$sample.id], taxa_are_rows = TRUE),
                         sample_data(bulk_run %>% column_to_rownames("sample.id")))
  
  # Calculate Bray-Curtis distance using the raw counts phyloseq object
  bray_curtis_raw <- phyloseq::distance(physeq_raw, method = "bray")
  
  # Also get the metadata correctly formatted for vegan functions
  sample_df_raw <- data.frame(sample_data(physeq_raw))
  
  # Ensure your grouping variable is a factor
  sample_df_raw$all_group <- factor(sample_df_raw$all_group, levels = order) %>% droplevels()
  
  # --- Step 2: Global PERMANOVA Test ---
  set.seed(42) # for reproducible permutations
  permanova_res_raw <- adonis2(bray_curtis_raw ~ all_group,
                               data = sample_df_raw,
                               permutations = 999)
  
  # --- Step 3: Check Group Dispersion ---
  beta_disp_raw <- betadisper(bray_curtis_raw, sample_df_raw$all_group)
  permutest_disp_raw <- permutest(beta_disp_raw, permutations = 999)
  
  print(permutest_disp_raw)
  
  # Only run if global PERMANOVA was significant
  if (permanova_res_raw$`Pr(>F)`[1] < 0.05) {
    pairwise_permanova_raw <- pairwise.adonis(bray_curtis_raw,
                                              factors = sample_df_raw$all_group,
                                              p.adjust.m = "BH") # Benjamini-Hochberg correction
  } else {
    print("Global PERMANOVA not significant, skipping pairwise tests.")
    pairwise_permanova_raw <- NULL
  }
  
  write_csv(pairwise_permanova_raw, file.path(outdir, "tables", paste0(run, "-bray-pairwise_permanova-antis.csv")))
  
  
  # Run PCoA using the raw counts phyloseq object and raw counts distance matrix
  pcoa_bray_raw <- ordinate(physeq_raw, method = "PCoA", distance = bray_curtis_raw)
  
  # Make the plot
  pcoa_plot_raw <- plot_ordination(physeq_raw, pcoa_bray_raw, color = "all_group") +
    labs(title=paste0("PCoA on Bray-Curtis (raw counts, ", run, ")"), colour = "All conditions") +
    geom_point(size = 2.5, alpha = 0.7) + # Adjust point appearance
    stat_ellipse(aes(group = all_group, color = all_group), # Use group aesthetic here too
                 type = "t", level = 0.95, linewidth=0.7) + # Ellipses for groups
    scale_color_manual(values = cb_palette_final) + # Use your defined colors
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          aspect.ratio = 1)
  
  # Extract global p-value and R2 from results in Step 2
  global_permanova_p <- permanova_res_raw$`Pr(>F)`[rownames(permanova_res_raw) == "Model"]
  global_permanova_r2 <- permanova_res_raw$R2[rownames(permanova_res_raw) == "Model"]
  permanova_label <- paste0("PERMANOVA (all_group):\nR2 = ", round(global_permanova_r2, 3),
                            "; p = ", format.pval(global_permanova_p, digits = 3))
  
  # Add annotation text
  pcoa_plot_annotated <- pcoa_plot_raw +
    annotate("text", x = Inf, y = -Inf, # Bottom-right corner (adjust if needed)
             label = permanova_label,
             hjust = 1.05, vjust = -0.5, # Adjust position
             size = 3.5, col = "black")
  
  ggsave(filename = file.path(outdir, "plots", paste0(run, "_bray_pcoa-antis.pdf")), 
         plot = pcoa_plot_annotated,width = 7, height = 5.5)
  
  # Convert distance object to matrix
  bc_matrix <- as.matrix(bray_curtis_raw)
  
  # Identify samples belonging to the baseline group 
  if(run == "Run1"){
    baseline_group_name <- "Uni_Typic_1.5_Anti" 
  }else{
    baseline_group_name <- "Uni_Typic_0_Anti"
  }
  
  baseline_samples <- rownames(sample_df_raw)[sample_df_raw$all_group == baseline_group_name]
  
  # --- Calculate average distance to baseline for each sample ---
  all_samples <- rownames(bc_matrix)
  dist_to_baseline <- numeric(length(all_samples))
  names(dist_to_baseline) <- all_samples
  
  for (sample_id in all_samples) {
    distances_to_baseline_group <- bc_matrix[sample_id, setdiff(baseline_samples, sample_id)]
    dist_to_baseline[sample_id] <- mean(distances_to_baseline_group, na.rm = TRUE)
  }
  
  # Create a data frame
  dist_to_baseline_df <- data.frame(
    sample.id = names(dist_to_baseline),
    Avg_Dist_To_Baseline = dist_to_baseline
  )
  
  metadata_with_dist <- sample_df_raw %>%
    tibble::rownames_to_column("sample.id") %>% 
    left_join(dist_to_baseline_df, by = "sample.id")
  
  metadata_with_dist$all_group <- factor(metadata_with_dist$all_group, levels = order)%>%
    droplevels()
  
  alpha_group_var <- "all_group" 
  y_var <- "Avg_Dist_To_Baseline" 
  facet_var <- NULL 
  
  dunn_test_dist <- NULL 
  
  # Kruskal-Wallis test
  stat_test_dist <- kruskal_test(data = metadata_with_dist, formula = reformulate(alpha_group_var, y_var))
  print(paste("Kruskal-Wallis test for", y_var, ":"))
  print(stat_test_dist)
  
  # Dunn's post-hoc test if Kruskal-Wallis is significant
  if (stat_test_dist$p < 0.05) {
    if (!is.null(facet_var)) { 
      dunn_test_dist <- dunn_test(data = metadata_with_dist, formula = reformulate(alpha_group_var, y_var), p.adjust.method = "BH") %>%
        add_xy_position(x = alpha_group_var, facet.by = facet_var) %>%
        add_significance()
    } else { 
      dunn_test_dist <- dunn_test(data = metadata_with_dist, formula = reformulate(alpha_group_var, y_var), p.adjust.method = "BH") %>%
        add_xy_position(x = alpha_group_var) %>%
        add_significance()
    }
  }
  
  dist_boxplot <- ggplot(metadata_with_dist, aes(x = !!sym(alpha_group_var), y = !!sym(y_var))) +
    geom_boxplot(outlier.shape = NA, aes(fill = !!sym(alpha_group_var))) +
    geom_jitter(width = 0.25, height = 0, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cb_palette_final) +
    labs(title = "Average Bray-Curtis distance to baseline group",
         subtitle = get_test_label(stat_test_dist, detailed = TRUE),
         x = "Group",
         y = "Avg. Bray-Curtis dist. to baseline") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  if (!is.null(facet_var)) {
    dist_boxplot <- dist_boxplot + facet_wrap(vars(!!sym(facet_var)))
  }
  
  if (exists("dunn_test_dist") && stat_test_dist$p < 0.05) {
    dist_boxplot <- dist_boxplot +
      ggpubr::stat_pvalue_manual(
        data = dunn_test_dist,
        label = "p.adj.signif",
        tip.length = 0.01,
        hide.ns = TRUE
      )
  }
  
  print(dist_boxplot)
  
  ggsave(filename = file.path(outdir, "plots", paste0(run, "_bray_boxplot-antis.pdf")), 
         plot = dist_boxplot,width = 7, height = 7)
  
  if(!is.null(dunn_test_dist)){
    write_csv(dunn_test_dist, 
              file.path(outdir, "tables", paste0(run, "_bray-dunn-test-dist-all-relative-to-baseline-antis.csv")))
  }
  
  # Make a Shannon index plot 
  mat_for_vegan <- t(mat[,bulk_run$sample.id])
  identical(bulk_run$sample.id, rownames(mat_for_vegan))
  
  shannon_indices <- vegan::diversity(mat_for_vegan, index = "shannon")
  
  shannon_df <- data.frame(
    sample.id = names(shannon_indices),
    Shannon = shannon_indices,
    row.names = names(shannon_indices) 
  )
  
  metadata_with_shannon <- left_join(bulk_run, shannon_df, by = "sample.id")
  metadata_with_shannon$all_group <- factor(metadata_with_shannon$all_group, levels = order)%>%
    droplevels()
  
  dunn_test_shannon <- NULL 
  
  stat_test_shannon <- kruskal_test(data = metadata_with_shannon, Shannon ~ all_group)
  
  if (stat_test_shannon$p < 0.05) {
    dunn_test_shannon <- dunn_test(data = metadata_with_shannon, Shannon ~ all_group, p.adjust.method = "BH") %>%
      add_xy_position(x = "all_group") %>%
      add_significance()
    print("Shannon Dunn's Post-Hoc Test Results (with facet positions):")
    
    write_csv(dunn_test_shannon, file.path(outdir, "tables", paste0(run, "_dunn-test-table-shannon_index-antis.csv")))
  }
  
  shannon_plot <- ggplot(metadata_with_shannon, aes(x = !!sym(alpha_group_var), y = Shannon)) +
    geom_boxplot(outlier.shape = NA, aes(fill = !!sym(alpha_group_var))) +
    geom_jitter(width = 0.25, height = 0, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cb_palette_final) +
    labs(title = "Shannon diversity comparison",
         x = "Group",
         y = "Shannon index (ln)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  if (!is.null(stat_test_shannon)) {
    shannon_plot <- shannon_plot + labs(subtitle = get_test_label(stat_test_shannon, detailed = TRUE))
    
    if (exists("dunn_test_shannon") && stat_test_shannon$p < 0.05) {
      shannon_plot <- shannon_plot +
        ggpubr::stat_pvalue_manual(
          dunn_test_shannon,
          label = "p.adj.signif", 
          tip.length = 0.01,
          hide.ns = TRUE 
        )
    }
  }
  
  ggsave(filename = file.path(outdir, "plots", paste0(run, "_shannon_plot-antis.pdf")), 
         plot = shannon_plot, width = 7, height = 7)
  
}

# Compare groups of interest

# Final filter from your script: Groups of interest from Run2
bulk_run <- bulk_anno %>%
  filter(Run == "Run2") %>%
  filter(all_group %in% c(
    "Uni_Const_0_Anti", "Uni_Typic_0_Anti",     # Pair 1
    "Uni_Const_1.5_Anti", "Uni_Typic_1.5_Anti", # Pair 2
    "Inf_Const_1.5_Anti", "Inf_Typic_1.5_Anti"  # Pair 3
  )) %>%
  mutate(all_group = factor(all_group))

my_comparisons_list <- list(
  c("Uni_Const_0_Anti", "Uni_Typic_0_Anti"),
  c("Uni_Const_1.5_Anti", "Uni_Typic_1.5_Anti"),
  c("Inf_Const_1.5_Anti", "Inf_Typic_1.5_Anti")
)

groups_in_analysis <- levels(bulk_run$all_group)
analysis_palette <- cb_palette_final[groups_in_analysis]

bray_subset <- vegan::vegdist(t(mat[, bulk_run$sample.id]), method = "bray")
group_counts_bray <- table(bulk_run$all_group)

set.seed(123)
all_pairwise_permanova <- pairwiseAdonis::pairwise.adonis2(
  bray_subset ~ all_group,
  data = bulk_run,
  permutations = 999,
  p.adjust.m = "BH"
)

comparison_names <- names(all_pairwise_permanova)[-1] 

results_list <- list()

for (comp_name in comparison_names) {
  comp_result_table <- all_pairwise_permanova[[comp_name]]
  f_stat <- comp_result_table$F[1]
  r_sq   <- comp_result_table$R2[1]
  p_val  <- comp_result_table$`Pr(>F)`[1]
  groups <- strsplit(comp_name, "_vs_")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  results_list[[comp_name]] <- data.frame(
    Group1 = group1, Group2 = group2, F.Model = f_stat,
    R2 = r_sq, p.value = p_val, comparison_name = comp_name
  )
}

target_pairs_sorted <- sapply(my_comparisons_list, function(p) paste(sort(p), collapse = "_vs_"))

all_results_df <- do.call(rbind, results_list)%>%
  filter(comparison_name %in% target_pairs_sorted)%>%
  select(Group1, Group2, comparison_name, F.Model, R2, p.value)

all_results_df$p.adjusted <- p.adjust(all_results_df$p.value, method = "BH")

write_csv(all_results_df, file.path(outdir, "tables", "Groups_of_interest_permanova.csv"))

# Shannon Index: Pairwise Wilcoxon for Specific Contrasts ---
mat_subset_t <- t(mat[, bulk_run$sample.id])
shannon_subset <- vegan::diversity(mat_subset_t, index = "shannon")

shannon_df_subset <- data.frame(sample.id = names(shannon_subset), Shannon = shannon_subset)
bulk_run_shannon <- bulk_run %>% left_join(shannon_df_subset, by = "sample.id")

all_pairwise_shannon <- rstatix::wilcox_test(
  data = bulk_run_shannon,
  formula = Shannon ~ all_group,
  p.adjust.method = "BH"
)

pairwise_shannon_select <- all_pairwise_shannon %>%
  mutate(current_pair_sorted = pmap_chr(list(group1, group2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  filter(current_pair_sorted %in% target_pairs_sorted) %>%
  select(-current_pair_sorted)

pairwise_shannon_select$p.adj <- p.adjust(pairwise_shannon_select$p, method = "BH")

pairwise_shannon_plotting <- pairwise_shannon_select %>%
  rstatix::add_xy_position(x = "all_group", data = bulk_run_shannon) 

write_csv(pairwise_shannon_plotting, file.path(outdir, "tables", "Groups_of_interest_wilcox_shannon.csv"))

shannon_boxplot_select <- ggplot(bulk_run_shannon, aes(x = all_group, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, aes(fill = all_group)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  scale_fill_manual(values = analysis_palette) + 
  ggpubr::stat_pvalue_manual(
    y.position = 3.5,
    pairwise_shannon_plotting, 
    label = "p.adj.signif",   
    tip.length = 0.01,
    hide.ns = T           
  ) +
  labs(
    title = "Shannon diversity (run 2 specific groups)",
    subtitle = "Pairwise comparisons: Wilcoxon test, BH adjustment",
    x = "Group",
    y = "Shannon index (ln)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_blank(),
    legend.position = "none" 
  ) +
  coord_cartesian(ylim = c(NA, max(bulk_run_shannon$Shannon, na.rm = TRUE) * 1.15)) 

print(shannon_boxplot_select)

ggsave(filename = file.path(outdir, "plots", "Run2_shannon_boxplot_specific_contrasts.pdf"),
       plot = shannon_boxplot_select, width = 8, height = 6)


# --- 4. Bray-Curtis Derived Boxplot (Avg. Distance to Baseline) ---
baseline_group_name_subset <- "Uni_Typic_0_Anti"
baseline_samples_subset <- bulk_run$sample.id[bulk_run$all_group == baseline_group_name_subset]

bc_matrix_subset <- as.matrix(bray_subset)
all_samples_subset <- bulk_run$sample.id
dist_to_baseline_subset <- numeric(length(all_samples_subset))
names(dist_to_baseline_subset) <- all_samples_subset

for (sample_id in all_samples_subset) {
  distances_to_baseline_group <- bc_matrix_subset[sample_id, setdiff(baseline_samples_subset, sample_id)]
  if(length(distances_to_baseline_group) > 0){
    dist_to_baseline_subset[sample_id] <- mean(distances_to_baseline_group, na.rm = TRUE)
  } else {
    dist_to_baseline_subset[sample_id] <- NA 
    warning(paste("Could not calculate distance to baseline for sample", sample_id))
  }
}

dist_to_baseline_df_subset <- data.frame(
  sample.id = names(dist_to_baseline_subset),
  Avg_Dist_To_Baseline = dist_to_baseline_subset
)

bulk_run_dist <- bulk_run %>%
  left_join(dist_to_baseline_df_subset, by = "sample.id")

stat_test_dist_subset <- kruskal_test(data = bulk_run_dist, Avg_Dist_To_Baseline ~ all_group)

if (stat_test_dist_subset$p < 0.05) {
  dunn_test_dist_subset_all <- dunn_test(data = bulk_run_dist, Avg_Dist_To_Baseline ~ all_group, p.adjust.method = "BH")
  
  dunn_test_dist_plotting <- dunn_test_dist_subset_all %>%
    mutate(current_pair_sorted = pmap_chr(list(group1, group2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
    filter(current_pair_sorted %in% target_pairs_sorted) %>%
    select(-current_pair_sorted) %>%
    rstatix::add_xy_position(x = "all_group", data = bulk_run_dist) 
  
  print("Filtered Dunn's Test (Avg Dist to Baseline) Results:")
  print(dunn_test_dist_plotting)
  
  dunn_test_dist_plotting$p.adj <-  p.adjust(dunn_test_dist_plotting$p, method = "BH")
  dunn_test_dist_plotting <- dunn_test_dist_plotting%>%
    mutate(p.adj.signif = replace(p.adj.signif, p.adj < 0.05, "*"))
  
  write_csv(dunn_test_dist_plotting, file.path(outdir, "tables", "Groups_of_interest_Dunn_test_Bray_dist.csv"))
  
} else {
  message("Overall Kruskal-Wallis test for Avg_Dist_To_Baseline not significant (p=", stat_test_dist_subset$p, "). Skipping Dunn's test and annotations.")
}

dist_boxplot_select <- ggplot(bulk_run_dist, aes(x = all_group, y = Avg_Dist_To_Baseline)) +
  geom_boxplot(outlier.shape = NA, aes(fill = all_group)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  scale_fill_manual(values = analysis_palette) + 
  labs(
    title = "Average Bray-Curtis distance to baseline (Uni_Typic_0_Anti)",
    subtitle = get_test_label(stat_test_dist_subset, detailed = TRUE), 
    x = "Group",
    y = "Avg. Bray-Curtis dist. to baseline"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

if (exists("dunn_test_dist_plotting") && stat_test_dist_subset$p < 0.05) {
  dist_boxplot_select <- dist_boxplot_select +
    ggpubr::stat_pvalue_manual(
      y.position = 1,
      dunn_test_dist_plotting,
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = TRUE
    ) +
    coord_cartesian(ylim = c(NA, max(bulk_run_dist$Avg_Dist_To_Baseline, na.rm = TRUE) * 1.3))
}

print(dist_boxplot_select)

ggsave(filename = file.path(outdir, "plots", "Run2_bray_dist_to_baseline_boxplot_specific_contrasts_dunn.pdf"),
       plot = dist_boxplot_select, width = 8, height = 6)


# Keep only the bacteria that were DE
DE <- lcpm[rownames(lcpm)%in% unique(toptables_signif$Gene),]

cpm_scaled <- DE[,bulk_anno$sample.id]%>%
  t()%>%
  scale()%>%
  t()

row_df <- data.frame(Feature.ID = rownames(cpm_scaled))%>%
  left_join(taxonomy, by = "Feature.ID")

rownames(cpm_scaled) <- row_df$Taxon

ha = HeatmapAnnotation(Treat = bulk_anno$all_group,
                       Days_since_infection = bulk_anno$Days_since_infection,
                       Infection_status = bulk_anno$Infection_status,
                       Light = bulk_anno$light,
                       Antibiotic = bulk_anno$Antibiotic,
                       col = list(Treat=cb_palette_final))

hm <- Heatmap(cpm_scaled, top_annotation = ha, name = "Row Z score",
              show_column_names = F, 
              cluster_columns = F,
              cluster_rows = T,
              column_split = bulk_anno$all_group,
              column_title_rot = 45,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(30, "cm"))

pdf(file.path(outdir, "plots", "DE_bacteria_heatmap.pdf"), width = 40, height = 20) 
draw(hm) 
dev.off()  

# Try a family level analysis ----
outdir2 <- file.path(".", "output_family_level_abundance")
dir.create(file.path(outdir2, "toptables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "toptreat"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "glimma", "volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "glimma", "MA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir2, "glimma", "mds"), recursive = TRUE, showWarnings = FALSE)

dt <- table_qza$data
mat <- as.matrix(dt)

anno_df <- data.frame(Feature.ID = rownames(mat))%>%
  left_join(taxonomy, by = "Feature.ID")%>%
  filter(!is.na(Family))

mat <- mat[anno_df$Feature.ID,]

counts_df <- mat %>%
  data.frame()%>%
  rownames_to_column(var = "Feature.ID")

tax_family <- anno_df %>%
  select(Feature.ID, Family)%>%
  left_join(counts_df, by = "Feature.ID")

family_counts <- tax_family %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

family_mat <- family_counts[,2:ncol(family_counts)]%>%
  as.matrix()

rownames(family_mat) <- family_counts$Family

family_mat <- family_mat[,bulk_anno$sample.id]
sum(colnames(family_mat) == bulk_anno$sample.id) == nrow(bulk_anno)

family_mat

table(bulk_anno$Run, bulk_anno$all_group)

bulk_anno$Acholeplasmataceae <- family_mat["f__Acholeplasmataceae",]
bulk_anno$RF39 <- family_mat["f__RF39",]

bulk_test <- bulk_anno%>%
  dplyr::select(RF39, Acholeplasmataceae, Run, subject, Days_since_infection)%>%
  filter(Acholeplasmataceae > 0)%>%
  filter(Days_since_infection == -14)

ggplot(data = bulk_test, aes(x = Run, y = Acholeplasmataceae))+
  geom_boxplot(outlier.alpha = 0)+
  facet_wrap(~Days_since_infection)+
  geom_jitter()

bulk_test <- bulk_anno%>%
  dplyr::select(RF39, Acholeplasmataceae, Run, subject, Days_since_infection)%>%
  filter(RF39 > 0)%>%
  filter(Days_since_infection == -14)

ggplot(data = bulk_test, aes(x = Run, y = RF39))+
  geom_boxplot(outlier.alpha = 0)+
  facet_wrap(~Days_since_infection)+
  geom_jitter()


counts <- DGEList(family_mat)

design <- model.matrix(~0 + all_group + Run, data = bulk_anno)
colnames(design) <- gsub("all_group", "", colnames(design))
rownames(design) <- rownames(counts$samples)

hist(rowSums(family_mat), breaks = 100000, xlim = c(0,100))

table(bulk_anno$all_group)

keep <- rowSums(family_mat>5)>3
table(keep)
rownames(counts)[!keep]
counts <- counts[keep,, keep.lib.sizes=FALSE]
dim(counts)

counts <- calcNormFactors(counts, method = "TMMwsp")

lcpm_family <- edgeR::cpm(counts, log = T, normalized.lib.sizes = TRUE)

lcpm_df_family <- lcpm_family%>%
  data.frame()%>%
  rownames_to_column("Family")

write_csv(lcpm_df_family, file.path(outdir2, "TMM_normalised_log2_CPMs_family.csv"))

mds_save <- file.path(outdir2, "glimma", "mds", "MDS.html")

htmlwidgets::saveWidget(glimmaMDS(counts, groups = bulk_anno$all_group,
                                  labels = bulk_anno), mds_save)

design2 <- model.matrix(~0 + all_group, data = bulk_anno)

lcpm_batch_removed_family <- removeBatchEffect(lcpm_family, 
                                               batch = bulk_anno$Run, 
                                               design = design2)

mds_save <- file.path(outdir2, "glimma", "mds", "MDS_run_batch_removed.html")

htmlwidgets::saveWidget(glimmaMDS(lcpm_batch_removed_family, groups = bulk_anno$all_group,
                                  labels = bulk_anno), mds_save)

colnames(design)

cont.matrix <- makeContrasts(
  Uni_Const_0_Anti_vs_Uni_Typic_0_Anti = Uni_Const_0_Anti - Uni_Typic_0_Anti,
  Uni_Const_1.5_Anti_vs_Uni_Typic_1.5_Anti = Uni_Const_1.5_Anti - Uni_Typic_1.5_Anti,
  Inf_Const_1.5_Anti_vs_Inf_Typic_1.5_Anti = Inf_Const_1.5_Anti - Inf_Typic_1.5_Anti,
  light_both_days = (Uni_Const_0_Anti + Inf_Const_1.5_Anti)/2 - 
    (Uni_Typic_0_Anti + Inf_Typic_1.5_Anti)/2,
  Uni_Const_0_Anti_vs_Uni_Const_0_No_anti = Uni_Const_0_Anti - Uni_Const_0_No_anti,
  Uni_Const_0_No_anti_vs_Uni_n14_No_anti = Uni_Const_0_No_anti - Uni_n14_No_anti,
  Inf_Typic_1.5_Anti_vs_Uni_Typic_1.5_Anti = Inf_Typic_1.5_Anti - Uni_Typic_1.5_Anti,
  levels = colnames(design))

colSums(cont.matrix)

fit <- voomLmFit(counts, design, plot=TRUE, block = bulk_anno$subject)

fit <- eBayes(fit)

toptable_run <- topTable(fit,coef="RunRun2",sort.by="p",number = Inf)%>%
  rownames_to_column("Gene")%>%
  write_csv(file.path(outdir2, "toptables", "Run_toptable.csv"))

vol_save <- file.path(outdir2, "glimma", "volcano", "Run_Volcano.html")

htmlwidgets::saveWidget(glimmaVolcano(fit,coef="RunRun2", main = "Run2",
                                      counts = cpm(counts, log = T),transform.counts	= "none",
                                      dge = counts, groups = bulk_anno$group_run), vol_save)

fit2 <- contrasts.fit(fit, contrasts=cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit, main="Final model: Mean-variance trend")
summa.fit <- decideTests(fit2)
summary(summa.fit)

tfit <- treat(fit2, fc=1.5)

for(contrast in colnames(cont.matrix)){
  
  toptable <- topTable(fit2,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(file.path(outdir2, "toptables", paste0(contrast, "_toptable.csv")))
  
  toptreat <- topTreat(tfit,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(file.path(outdir2, "toptreat", paste0(contrast, "_toptreat.csv")))
  
  vol_save <- file.path(outdir2, "glimma", "volcano", paste0(contrast, "_Volcano.html"))
  
  htmlwidgets::saveWidget(glimmaVolcano(fit2, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = cpm(counts, log = T),transform.counts	= "none",
                                        dge = counts, groups = bulk_anno$all_group), vol_save)
  
  ma_save <- file.path(outdir2, "glimma", "MA", paste0(contrast, "_MA.html"))
  
  htmlwidgets::saveWidget(glimmaMA(fit2, coef = contrast,main = gsub("_"," ",contrast),
                                   counts = cpm(counts, log = T),transform.counts	= "none",
                                   dge = counts, groups = bulk_anno$all_group), ma_save)
  
}

files_to_remove2 <- list.files(file.path(outdir2, "glimma"), pattern = "_files$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
unlink(files_to_remove2, recursive = TRUE)

all_toptables <- list.files(file.path(outdir2, "toptables"), full.names = TRUE)

tt_list <- list()
for(i in seq_along(all_toptables)){
  
  contrast <- gsub("_toptable\\.csv$", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
}

if(length(tt_list) > 0) {
  toptables_compiled <- bind_rows(tt_list)
  
  toptables_signif <- toptables_compiled %>%
    filter(adj.P.Val < 0.05)%>%
    arrange(adj.P.Val)%>%
    write_csv(file.path(outdir2, "Compiled_toptables_significant_families.csv"))
}

all_toptreat <- list.files(file.path(outdir2, "toptreat"), full.names = TRUE)

tt_list <- list()
for(i in seq_along(all_toptreat)){
  
  contrast <- gsub("_toptreat\\.csv$", "", basename(all_toptreat[i]))
  
  tt <- read_csv(all_toptreat[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
}

if(length(tt_list) > 0) {
  toptreat_compiled <- bind_rows(tt_list)
  
  toptreat_signif <- toptreat_compiled %>%
    filter(adj.P.Val < 0.05)%>%
    arrange(adj.P.Val)%>%
    write_csv(file.path(outdir2, "Compiled_toptreat_significant_families.csv"))
}

conts <- c("Inf_Const_1.5_Anti_vs_Inf_Typic_1.5_Anti","Uni_Const_1.5_Anti_vs_Uni_Typic_1.5_Anti")

for (cont in conts){
  
  toptable_plot <- toptables_compiled %>%
    filter(contrast == cont)
  
  if(nrow(toptable_plot) == 0) {
    message(paste("No data found for contrast:", cont, "- Skipping plot."))
    next
  }
  
  plot_data <- toptable_plot %>%
    mutate(Significance = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Increase",
      adj.P.Val < 0.05 & logFC < 0 ~ "Decrease",
      TRUE ~ "Not differentially abundant"
    )) %>%
    mutate(Significance = factor(Significance, levels = c("Decrease", "Not differentially abundant", "Increase"))) %>%
    mutate(Clean_Family = gsub("^f__", "", Gene))
  
  volcano_colors <- c(
    "Increase" = "red", 
    "Decrease" = "blue", 
    "Not differentially abundant" = "grey80"
  )
  
  custom_volcano <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(fill = Significance), alpha = 0.7, size = 2, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(name = "Differential abundance", values = volcano_colors) +
    labs(
      x = expression(Log[2]*"FC constant light relative to 12-hour light/12-hour dark"),
      y = expression(-Log[10](P[adj]))
    ) +
    theme_bw(base_size = 14) 
  
  sig_data <- filter(plot_data, Significance != "Not differentially abundant")
  
  if(nrow(sig_data) > 0) {
    custom_volcano <- custom_volcano +
      geom_text_repel(
        data = sig_data,
        aes(label = Clean_Family),
        size = 3.5,
        max.overlaps = 20,
        show.legend = FALSE 
      )
  }
  ggsave(filename = file.path(outdir2, "plots", paste0(cont, "_custom_Volcano.pdf")),
         plot = custom_volcano, width = 8, height = 6)
  
}

# Family level composition plot
family_counts_plot <- family_counts%>%
  pivot_longer(cols = -Family, names_to = "sample.id", values_to = "Count")%>%
  group_by(sample.id)%>%
  mutate(Total = sum(Count))%>%
  ungroup()%>%
  mutate(pct_of_total = Count/Total*100)%>%
  left_join(bulk_anno, by = "sample.id")%>%
  mutate(all_group = factor(all_group, levels = order))%>%
  mutate(Family_plot = replace(Family, pct_of_total < 1, "Other"))%>%
  write_csv(file.path(outdir2, "tables", "Family level proportions table.csv"))

cbp32 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", 
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", 
  "#AA4499", "#661100", "#6699CC", "#888888", "#FAD0C3", "#79A5C2", "#A1EDB8", "#FECA81", 
  "#B496CF", "#796A52", "#E54F6D", "#56E3B8", "#2E86AB", "#F5B82E", "#721817", "#1E3F20",
  "grey"
)

family_names <- c(unique(family_counts_plot$Family),"Other")

family_colors <- setNames(cbp32[1:length(unique(family_names))], unique(family_names)) 

family_counts_plot_r1 <- family_counts_plot%>%
  filter(Run == "Run1")

family_counts_plot_r2 <- family_counts_plot%>%
  filter(Run == "Run2")

family_plot_combined <- ggplot(data = family_counts_plot, aes(x = sample.id, y = pct_of_total, fill = Family))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total")+
  scale_fill_manual(values = family_colors)

family_plot_combined

ggsave(filename = file.path(outdir2, "plots", "family-abundance-plot-relative_all.pdf"), 
       plot = family_plot_combined,
       width = 18, height = 7)

family_plot_r1 <- ggplot(data = family_counts_plot_r1, aes(x = sample.id, y = pct_of_total, fill = Family))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total")+
  scale_fill_manual(values = family_colors)

ggsave(filename = file.path(outdir2, "plots", "family-abundance-plot-relative-r1_all.pdf"), 
       plot = family_plot_r1,
       width = 18, height = 7)

family_plot_r1 <- ggplot(data = family_counts_plot_r1, aes(x = sample.id, y = pct_of_total, fill = Family_plot))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total", fill = "Family")+
  scale_fill_manual(values = family_colors)

ggsave(filename = file.path(outdir2, "plots", "family-abundance-plot-relative-r1.pdf"), 
       plot = family_plot_r1,
       width = 18, height = 7)

family_plot_r2 <- ggplot(data = family_counts_plot_r2, aes(x = sample.id, y = pct_of_total, fill = Family))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total")+
  scale_fill_manual(values = family_colors)

ggsave(filename = file.path(outdir2, "plots", "family-abundance-plot-relative-r2_all.pdf"), 
       plot = family_plot_r2,
       width = 18, height = 7)

family_plot_r2 <- ggplot(data = family_counts_plot_r2, aes(x = sample.id, y = pct_of_total, fill = Family_plot))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total", fill = "Family")+
  scale_fill_manual(values = family_colors)

ggsave(filename = file.path(outdir2, "plots", "family-abundance-plot-relative-r2.pdf"), 
       plot = family_plot_r2,
       width = 18, height = 7)

# Repeat at the phylum level
tax_phylum <- anno_df %>%
  select(Feature.ID, Phylum)%>%
  left_join(counts_df, by = "Feature.ID")

phylum_counts <- tax_phylum %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

phylum_counts_plot <- phylum_counts%>%
  pivot_longer(cols = -Phylum, names_to = "sample.id", values_to = "Count")%>%
  group_by(sample.id)%>%
  mutate(Total = sum(Count))%>%
  ungroup()%>%
  mutate(pct_of_total = Count/Total*100)%>%
  left_join(bulk_anno, by = "sample.id")%>%
  mutate(all_group = factor(all_group, levels = order))

phylum_counts_plot_r1 <- phylum_counts_plot%>%
  filter(Run == "Run1")

phylum_counts_plot_r2 <- phylum_counts_plot%>%
  filter(Run == "Run2")

phylum_plot_r1 <- ggplot(data = phylum_counts_plot_r1, 
                         aes(x = sample.id, y = pct_of_total, fill = Phylum))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total")

phylum_plot_r2 <- ggplot(data = phylum_counts_plot_r2, 
                         aes(x = sample.id, y = pct_of_total, fill = Phylum))+
  geom_bar(stat = "identity")+
  facet_wrap(~all_group, scales = "free", nrow = 1,
             strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text.x = element_text(
          angle = 90))+
  labs(x = "Condition", y = "% of total")

ggsave(filename = file.path(outdir2, "plots", "phylum-abundance-plot_run_1.pdf"), 
       plot = phylum_plot_r1,
       width = 16, height = 5)

ggsave(filename = file.path(outdir2, "plots", "phylum-abundance-plot_run_2.pdf"), 
       plot = phylum_plot_r2,
       width = 16, height = 5)

# Family counts plot over all conditions
if (exists("toptables_signif")) {
  family_counts_plot_signif <- family_counts_plot%>%
    filter(Family %in% toptables_signif$Gene)%>%
    mutate(Family = gsub("^[dfkpcogs]__", "", Family)) 
  
  if(nrow(family_counts_plot_signif) > 0) {
    family_abundance_plot <- ggplot(family_counts_plot_signif, aes(x = all_group, y = pct_of_total)) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "grey60") +
      stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
      stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 1, color = "red") +
      facet_wrap(~ Family, scales = "free_y", ncol = 3) + 
      labs(
        title = "Relative abundance of all significant families across groups",
        x = "Condition",
        y = "Relative abundance (%)" 
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        strip.text = element_text(size = 8, face = "bold"), 
        strip.background = element_rect(fill = "grey90", color = NA)
      )
    
    family_abundance_plot
    
    ggsave(filename = file.path(outdir2, "plots", "family-abundance-line-plot.pdf"), plot = family_abundance_plot,
           width = 11,height = 10)
  }
}

# Keep only the bacteria that were DE
design2 <- model.matrix(~0 + all_group, data = bulk_anno)
lcpm_batch_removed_family <- removeBatchEffect(lcpm_family, 
                                               batch = bulk_anno$Run, 
                                               design = design2)

if(exists("toptables_signif") && nrow(toptables_signif) > 0) {
  DE <- lcpm_batch_removed_family[rownames(lcpm_batch_removed_family)%in% unique(toptables_signif$Gene),]
  
  cpm_scaled <- DE[,bulk_anno$sample.id]%>%
    t()%>%
    scale()%>%
    t()
  
  bulk_anno$all_group <- factor(bulk_anno$all_group, levels = order)
  
  ha = HeatmapAnnotation(Treat = bulk_anno$all_group,
                         Days_since_infection = bulk_anno$Days_since_infection,
                         Infection_status = bulk_anno$Infection_status,
                         Light = bulk_anno$light,
                         Antibiotic = bulk_anno$Antibiotic,
                         col = list(Treat=cb_palette_final))
  
  hm <- Heatmap(cpm_scaled, top_annotation = ha, name = "Row Z score",
                show_column_names = F, 
                cluster_columns = F,
                cluster_rows = T,
                column_split = bulk_anno$all_group,
                column_title_rot = 45,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 10),
                row_names_max_width = unit(30, "cm"))
  
  pdf(file.path(outdir2, "plots", "DE_bacteria_heatmap_family.pdf"), width = 20, height = 10) 
  draw(hm) 
  dev.off()  
}

# Make a plotting object.
family_long <- lcpm_df_family%>%
  pivot_longer(cols = -Family, names_to = "sample.id", values_to = "log2CPM")%>%
  left_join(bulk_anno, by = "sample.id")%>%
  write_csv(file.path(outdir2, "Family_log2CPMs_ggplot2_format.csv"))

# ---------------------------------------------------------------------------
# Session Info Export
# ---------------------------------------------------------------------------
session_dir <- file.path(".", "Code_for_submission", "sessioninfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(session_dir, "Session_info_abundance.txt"))
