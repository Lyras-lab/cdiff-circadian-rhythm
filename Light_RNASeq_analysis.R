library(tidyverse)
library(Rsubread)
library(qs)
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

# FIX: Use file.path for robust base directories
outdir <- file.path(".", "output_RNASeq")
indir <- file.path(".", "Transcriptomics_inputs")

# ---------------------------------------------------------------------------
# Alignment code is commented out but was run for the paper. 
# Uncomment and modify for GEO counts to rerun.
# ---------------------------------------------------------------------------

# Grab the raw data from the SRA

# Where the reference data is from
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.primary_assembly.genome.fa.gz

# Check the files (all good)
# md5sum --c md5.md5 > md5_check_result_on_labserver.txt

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
# # All 50 samples
# sum(table(files_df$Sample) == 2) == length(unique(files_df$Sample))
# 
# # Basename of the index
# index <- "./RCm38.primary_assembly.genome"
# 
# f1 <- files_df %>%
#   filter(grepl("_1.fq.gz",files))
# 
# f2 <- files_df %>%
#   filter(grepl("_2.fq.gz",files))
# 
# # Sanity check
# f2$Sample == f1$Sample

# Align all at once so I don't need to reload the index
# align.stat <- Rsubread::align(index = index,
#                                 readfile1 = f1$files,
#                                 readfile2 = f2$files,
#                                 output_file = f1$Outfile,
#                                 nthreads = 40)

# List the output files for the organoid data
# bams <- list.files(outdir, pattern = "*.bam$", recursive = TRUE, full.names = TRUE)
# 
# # Run featurecounts over the aligned bams
# fc_bam <- featureCounts(files = bams,
#                           annot.inbuilt = "mm10",
#                           annot.ext = NULL,
#                           isPairedEnd = TRUE,
#                           # Looks like the library is unstranded based on count assignment using stranded options
#                           strandSpecific = 0,
#                           nthreads = 20)

# saveRDS(fc_bam, file.path(indir, "featurecounts_object.rds"))

# ---------------------------------------------------------------------------
# Load existing data
# ---------------------------------------------------------------------------

fc_bam <- readRDS(file.path(indir, "featurecounts_object.rds"))
stats <- fc_bam$stat

# Grab the annotation for the entrezgene IDs do not rerun unless you want updated annotations
# ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# attrs <- listAttributes(ensembl) %>%
#   data.frame()
# annot <- getBM(c("ensembl_gene_id","entrezgene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=ensembl)
# 
# mouse_mart <- annot %>%
#   mutate(entrezgene_id = as.character(entrezgene_id)) %>%
#   filter(!is.na(entrezgene_id)) %>%
#   filter(!duplicated(entrezgene_id))

# Get the counts out of the object
counts <- fc_bam$counts

# counts[1:5,1:5]
# 
# annot_df_mouse <- data.frame(entrezgene_id = rownames(counts)) %>%
#   left_join(mouse_mart)

# saveRDS(annot_df_mouse, file.path(indir, "gene_annotation.rds"))

annot_df_mouse <- readRDS(file.path(indir, "gene_annotation.rds"))

# counts2 <- counts
# colnames(counts2) <- gsub("_aligned.bam", "", colnames(counts2))

# counts_and_anno_csv_for_degust <- counts2 %>%
#   data.frame() %>%
#   rownames_to_column("entrezgene_id") %>%
#   left_join(annot_df_mouse) %>%
#   write_csv(file.path(outdir, "Counts_and_anno_for_degust.csv"))

# saveRDS(counts2, file.path(indir, "counts.rds"))

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

# https://www.biostars.org/p/9567892/
# mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
# 
# # separate human and mouse 
# mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
# human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
# 
# mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
# 
# # remove some columns
# mouse <- mouse[,c(1,4)]
# human <- human[,c(1,4)]
# 
# # merge the 2 dataset  (note that the human list is longer than the mouse one)
# mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) %>%
#   filter(!duplicated(Symbol.x))
# 
# sum(duplicated(mh_data$Symbol.x))
# 
# mouse_genes <- rnaseq$genes$HUGOSymbol %>%
#   data.frame()%>%
#   dplyr::rename(Symbol.x = 1)%>%
#   left_join(mh_data)%>%
#   write_csv(file.path(indir, "GSEA", "Human_mouse_gene_mappings.csv"))

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
    rownames_
