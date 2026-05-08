library(tidyverse)
library(Rsubread)
library(limma)
library(edgeR)
library(Glimma)
library(biomaRt)
library(GSA)
library(parallel)
library(GSVA)

# Colours for groups
contrast_levs <- c("Uninfected_Typicallight", "Uninfected_Constantlight","Infected_Typicallight", "Infected_Constantlight")
group_colors <- c("Uninfected_Typicallight" = "#0066CC", "Uninfected_Constantlight" = "#99CCFF","Infected_Typicallight" = "#CC5200", "Infected_Constantlight" = "#FFB366")

# FIX: Use file.path for robust base directories. Remove trailing slashes.
outdir <- file.path(".", "output_RNASeq")
indir <- file.path(".", "Transcriptomics_inputs")

# ---------------------------------------------------------------------------
# Alignment code is commented out but was run for the paper. 
# ---------------------------------------------------------------------------

# FIX: Use dir.create instead of system("mkdir ...")
# dir.create(file.path(outdir, "bams"), recursive = TRUE, showWarnings = FALSE)
# 
# # Find the fastq files
# files <- list.files(file.path(".", "raw_data"), pattern = "*.fq.gz", full.names = TRUE, recursive = TRUE)
# 
# files_df <- data.frame(files) %>%
#   mutate(Sample = basename(files)) %>%
#   mutate(Sample = gsub("_1.fq.gz|_2.fq.gz", "", Sample)) %>%
#   mutate(Outfile = file.path(outdir, "bams", paste0(Sample, "_aligned.bam"))) %>%
#   write_csv(file.path(outdir, "Metadata.csv"))
# 
# # ... (rest of commented alignment code) ...

# ---------------------------------------------------------------------------
# Load existing data
# ---------------------------------------------------------------------------

fc_bam <- readRDS(file.path(indir, "featurecounts_object.rds"))
stats <- fc_bam$stat

# Get the counts out of the object
counts <- fc_bam$counts

annot_df_mouse <- readRDS(file.path(indir, "gene_annotation.rds"))
counts2 <- readRDS(file.path(indir, "counts.rds"))

sample_md <- read_csv(file.path(indir, "T282 RNA seq analysis.csv"))

metadata <- data.frame(Bam_name = colnames(counts2)) %>%
  mutate(`Sample name on tube` = gsub("-.*", "", Bam_name)) %>%
  left_join(sample_md, by = "Sample name on tube") %>%
  mutate(light_inf = paste0(`Infected or uninfected`, "_", `Constant light or typical light cycle`)) %>%
  mutate(light_inf = gsub(" ", "", light_inf)) %>%
  mutate(light_inf = factor(light_inf, levels = contrast_levs))

# Check everything is ordered at the sample 
sum(colnames(counts2) == metadata$Bam_name) == ncol(counts2)
# Check everything is ordered at the gene level
sum(rownames(counts2) == annot_df_mouse$entrezgene_id) == nrow(counts2)

# Make our limma object
rnaseq <- DGEList(counts2, genes = annot_df_mouse)
design <- model.matrix(~0 + light_inf, data = metadata)
colnames(design) <- gsub("light_inf", "", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

# Limma filter by expression
keep <- filterByExpr(rnaseq, design = design)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

# Apply TMM normalisation
rnaseq <- calcNormFactors(rnaseq, method = "TMM")
lcpm <- cpm(rnaseq, log = TRUE)

# FIX: Create directories natively
dir.create(file.path(outdir, "toptables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "MA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "mds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots", "barcode_plots"), recursive = TRUE, showWarnings = FALSE)

mds_save <- file.path(outdir, "glimma", "mds", "MDS.html")
htmlwidgets::saveWidget(glimmaMDS(rnaseq, groups = metadata$light_inf, labels = metadata), mds_save)

# Normalise and fit linear model
v <- voomWithQualityWeights(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(
  Light_Effect_Infected     = Infected_Constantlight - Infected_Typicallight,
  Light_Effect_Uninfected   = Uninfected_Constantlight - Uninfected_Typicallight,
  Infection_Effect_Constant = Infected_Constantlight - Uninfected_Constantlight,
  Infection_Effect_Typical  = Infected_Typicallight - Uninfected_Typicallight,
  Light_both = (Infected_Constantlight + Uninfected_Constantlight)/2 - (Infected_Typicallight + Uninfected_Typicallight)/2,
  Infection_Light_Interaction = (Infected_Constantlight - Uninfected_Constantlight) - (Infected_Typicallight - Uninfected_Typicallight),
  levels = design
)

fit <- contrasts.fit(fit, contrasts=cont.matrix)
fit <- eBayes(fit)
summa.fit <- decideTests(fit)

de_summary <- data.frame(summary(summa.fit), check.names = FALSE) %>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3) %>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig"))) %>%
  arrange(`Direction if significant`, `Number of genes`) %>%
  write_csv(file.path(outdir, "Significant_genes_summary.csv"))

signif <- de_summary %>% filter(`Direction if significant` != "NotSig")

plot_summary <- de_summary %>%
  filter(`Direction if significant` != "NotSig") %>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] * -1)) %>%
  arrange(`Direction if significant`, -`Number of genes`) %>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`)) +
  geom_bar(stat = "identity")

ggsave(file.path(outdir, "plots", "Signficant genes.pdf"), width = 7, height = 5)

for(contrast in colnames(cont.matrix)){
  toptable <- topTable(fit, coef=contrast, sort.by="p", number = Inf) %>%
    rownames_to_column("Gene") %>%
    write_csv(file.path(outdir, "toptables", paste0(contrast, "_toptable.csv")))
  
  vol_save <- file.path(outdir, "glimma", "volcano", paste0(contrast, "_Volcano.html"))
  htmlwidgets::saveWidget(glimmaVolcano(fit, coef = contrast, main = gsub("_"," ",contrast),
                                        counts = lcpm, transform.counts = "none",
                                        dge = rnaseq, groups = metadata$light_inf), vol_save)
  
  ma_save <- file.path(outdir, "glimma", "MA", paste0(contrast, "_MA.html"))
  htmlwidgets::saveWidget(glimmaMA(fit, coef = contrast, main = gsub("_"," ",contrast),
                                   counts = lcpm, transform.counts = "none",
                                   dge = rnaseq, groups = metadata$light_inf), ma_save)
}

# FIX: Use R's native unlink instead of system("rm") to remove directories
files_to_remove <- list.files(file.path(outdir, "glimma"), pattern = "_files$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
unlink(files_to_remove, recursive = TRUE)

get_geneset_list <- function(collection){
  gene_set <- quiet(GSA.read.gmt(collection))
  gene_set_formatted <- gene_set$genesets
  gene_set_formatted <- lapply(gene_set_formatted, function(x) x[nzchar(x)])
  names(gene_set_formatted) <- gene_set$geneset.names
  return(gene_set_formatted)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Gene set collections
gsea_gmt_dir <- file.path(indir, "GSEA")
collections <- list.files(gsea_gmt_dir, full.names = TRUE, pattern = "*.symbols.gmt")
collections <- c(collections, 
                 file.path(gsea_gmt_dir, "HC_Select Hs.Signalling_Geneset_combined_2024.Hs.gmt"),
                 file.path(gsea_gmt_dir, "c2.cp.kegg.v2023.1.Hs.symbols.gmt"),
                 file.path(gsea_gmt_dir, "HC_Select_Hs.CellType_Geneset_combined_2024.Hs.gmt"),
                 file.path(gsea_gmt_dir, "HC_Select_Hs.Oncogenic_Geneset_combined_2024.Hs.gmt"))

gene_mappings <- read_csv(file.path(gsea_gmt_dir, "Human_mouse_gene_mappings.csv"))

collections_all <- list()
for(i in 1:length(collections)){
  collection <- collections[i]
  if(grepl("HC_Select|c2.cp.kegg.v2023.1.Hs.symbols", collection)){
    human <- get_geneset_list(collection)
    for(i2 in 1:length(human)){
      gene_set <- human[[i2]]
      gene_set2 <- gene_mappings$Symbol.x[gene_mappings$Symbol.y %in% gene_set]
      collections_all[[names(human)[i2]]] <- gene_set2
    }
  } else {
    collections_all <- c(collections_all, get_geneset_list(collection))
  }
}

dir.create(file.path(outdir, "gsea", "camera"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "gsea", "fry"), recursive = TRUE, showWarnings = FALSE)

names_df <- data.frame(names(collections_all))
indexed <- ids2indices(gene.sets = collections_all, identifiers = rnaseq$genes$mgi_symbol, remove.empty=TRUE)

for(contrast in colnames(cont.matrix)){
  camera_result <- camera(y = v, index = indexed, design = design, contrast = cont.matrix[,contrast]) %>%
    rownames_to_column("Gene set") %>%
    dplyr::select(`Gene set`, "NGenes" , "Direction", "PValue", "FDR") %>%
    filter(FDR <= 0.05) %>%
    mutate(Contrast= contrast) %>%
    write_csv(file.path(outdir, "gsea", "camera", paste0(contrast, ".csv")))
  
  fry_result <- fry(y = v, index = indexed, design = design, contrast = cont.matrix[,contrast]) %>%
    rownames_to_column("Gene set") %>%
    dplyr::select(`Gene set`, "NGenes" , "Direction", "PValue", "FDR") %>%
    filter(FDR <= 0.05) %>%
    mutate(Contrast= contrast) %>%
    write_csv(file.path(outdir, "gsea", "fry", paste0(contrast, ".csv")))
}

# Compile the camera results
all_camera <- list.files(file.path(outdir, "gsea", "camera"), full.names = TRUE)
clist <- list()
for(i in seq_along(all_camera)){
  contrast_name <- gsub("\\.csv$", "", basename(all_camera[i]))
  tt <- read_csv(all_camera[i], col_types = cols(.default = "c")) %>%
    mutate(contrast = contrast_name) %>%
    dplyr::select(-Contrast)
  clist[[i]] <- tt
}

if(length(clist) > 0) {
  camera_compiled <- bind_rows(clist) %>%
    mutate(FDR = as.numeric(FDR)) %>%
    arrange(FDR) %>%
    write_csv(file.path(outdir, "Compiled_significant_gene_sets_camera.csv"))
}

# Compile Fry results
all_fry <- list.files(file.path(outdir, "gsea", "fry"), full.names = TRUE)
clist <- list()
for(i in seq_along(all_fry)){
  contrast_name <- gsub("\\.csv$", "", basename(all_fry[i]))
  tt <- read_csv(all_fry[i], col_types = cols(.default = "c")) %>%
    mutate(contrast = contrast_name) %>%
    dplyr::select(-Contrast)
  clist[[i]] <- tt
}

if(length(clist) > 0) {
  fry_compiled <- bind_rows(clist) %>%
    mutate(FDR = as.numeric(FDR)) %>%
    arrange(FDR) %>%
    write_csv(file.path(outdir, "Compiled_significant_gene_sets_fry.csv"))
}

# Compile the toptables
all_toptables <- list.files(file.path(outdir, "toptables"), full.names = TRUE)
tt_list <- list()
for(i in seq_along(all_toptables)){
  # FIX: Properly remove the entire suffix so contrast names remain clean
  contrast_name <- gsub("_toptable\\.csv$", "", basename(all_toptables[i]))
  tt <- read_csv(all_toptables[i]) %>%
    mutate(contrast = contrast_name)
  tt_list[[i]] <- tt
}

if(length(tt_list) > 0) {
  toptables_compiled <- bind_rows(tt_list)
  toptables_signif <- toptables_compiled %>%
    filter(adj.P.Val <= 0.05) %>%
    arrange(adj.P.Val) %>%
    write_csv(file.path(outdir, "Compiled_toptables_significant_genes.csv"))
}

gene_sets <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_MYC_TARGETS_V2", 
               "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
               "GOBP_PROTEIN_REFOLDING", "GOBP_REGULATION_OF_TYPE_I_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
               "GOBP_CHROMOSOME_SEPARATION")

# Barcode plots
gene_set <- quiet(GSA.read.gmt(file.path(indir, "GSEA", "mh.all.v2023.2.Mm.symbols.gmt")))
gene_set_formatted_hm <- gene_set$genesets
names(gene_set_formatted_hm) <- gene_set$geneset.names

gene_set <- quiet(GSA.read.gmt(file.path(indir, "GSEA", "m5.go.v2023.2.Mm.symbols.gmt")))
gene_set_formatted_go <- gene_set$genesets
names(gene_set_formatted_go) <- gene_set$geneset.names

gene_set_formatted <- c(gene_set_formatted_hm, gene_set_formatted_go)
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rnaseq$genes$mgi_symbol, remove.empty=TRUE)

for(contrast in colnames(cont.matrix)){
  for(gene_set_name in gene_sets){
    save_name <- file.path(outdir, "plots", "barcode_plots", paste0(contrast, "_", gene_set_name, ".pdf"))
    index <- indexed[[gene_set_name]]
    
    pdf(save_name, width = 7, height = 4)
    barcodeplot(fit$t[,contrast],
                index=index,
                labels = c("Down","Up"),
                main= paste0(gene_set_name, "\nN genes=", length(index)),
                xlab = paste0("Limma t statistic: ", contrast))
    dev.off()
  }
}

# GSVA analysis and plotting
lcpm_gene <- cpm(rnaseq, log = TRUE)
rownames(lcpm_gene) <- rnaseq$genes$mgi_symbol
lcpm_gene <- lcpm_gene[!is.na(rownames(lcpm_gene)), ]
lcpm_gene <- lcpm_gene[!duplicated(rownames(lcpm_gene)), ]

gene_sets_GSVA <- gene_set_formatted[names(gene_set_formatted) %in% gene_sets]
param <- gsvaParam(exprData = lcpm_gene, geneSets = gene_sets_GSVA)
gsva.es <- gsva(param, verbose=FALSE)

# FIX: Replace deprecated `gather()` with `pivot_longer()`
scores <- data.frame(gsva.es, check.rows = FALSE, check.names = FALSE) %>%
  rownames_to_column("Gene_set") %>%
  pivot_longer(cols = -Gene_set, names_to = "Bam_name", values_to = "GSVA_score") %>%
  mutate(Gene_set = gsub("_", " ", Gene_set)) %>%
  mutate(Gene_set = gsub("HALLMARK ", "", Gene_set)) %>%
  left_join(metadata, by = "Bam_name") %>%
  mutate(light_inf = factor(light_inf, levels = contrast_levs)) 

GSVA <- ggplot(data = scores, aes(x = light_inf, y = GSVA_score, colour = light_inf))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(height = 0, width = 0.3)+
  labs(x = "Light/infection", y = "GSVA score", title = "GSVA scores key gene sets", colour = "Light/infection")+
  facet_wrap(~Gene_set, ncol = 3)+
  scale_colour_manual(values = group_colors)+
  guides(colour = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = GSVA, filename = file.path(outdir, "plots", "GSVA_boxplot_key_gene_sets.pdf"), 
       width = 200, height = 170, units = "mm")

writeLines(capture.output(sessionInfo()), file.path(".", "sessioninfo", "Session_info_transcriptomics.txt"))
