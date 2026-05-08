library(tidyverse)
library(limma)
library(edgeR)
library(Glimma)
library(KEGGREST)
library(stringr)
library(ComplexHeatmap)
library(uwot)

# Set ggplot2 themes for the paper
blank_theme <- theme_bw(base_size = 7)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7))

# Looking at this for analysis
# https://support.bioconductor.org/p/64484/

condition_colors <- c(
  "Uni_Typical_light" = "#56B4E9",   # Light Blue
  "Uni_Const_light" = "#0072B2",   # Dark Blue
  "Inf_Typical_light" = "#D55E00",   # Reddish Orange
  "Inf_Const_light" = "#E69F00"   # Orange
)

# FIX: Use file.path for robust base directories
outdir <- file.path(".", "output_bile_acids")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE) # Ensure plots dir exists

# Read in the protein data
acids <- read_csv(file.path(".", "Bile_acids_input", "Bile_acids.csv"))

acids_mat <- as.matrix(acids[,2:ncol(acids)])
rownames(acids_mat) <- acids$Filename
colnames(acids_mat) <- gsub(",| ", "-", colnames(acids_mat))

acids_mat <- t(acids_mat)

md <- data.frame(Sample = colnames(acids_mat))%>%
  mutate(Infection = ifelse(grepl("Cdifficile",Sample), "Inf", "Uni"))%>%
  mutate(Light = ifelse(grepl("Typical_light",Sample), "Typical_light", "Const_light"))%>%
  mutate(Infection_light = paste0(Infection, "_", Light))%>%
  mutate(Infection_light = replace(Infection_light, grepl("QC",Sample), "QC"))%>%
  mutate(Infection_light = replace(Infection_light, grepl("Blank",Sample), "Blank"))%>%
  filter(!Infection_light %in% c("QC", "Blank"))%>%
  mutate(Infection_light = factor(Infection_light, levels = names(condition_colors)))

# Drop the peaks that don't appear consistently in at least one condition (or 5 samples)
table(md$Infection_light)

acids_mat <- acids_mat[,md$Sample]

colnames(acids_mat) == md$Sample

# Log transfrom the counts
log2_spectra <- log2(acids_mat + 1)

# 19 samples
dim(log2_spectra)

# Set up a design matrix
design <- model.matrix(~0+ Infection_light, data = md)
colnames(design) <- gsub("Infection_light", "", colnames(design))
rownames(design) <- colnames(acids_mat)

hist(log2_spectra, breaks = 300)
# "Normalizes expression intensities so that the intensities 
# or log-ratios have similar distributions across a set of arrays."
normalized <- normalizeBetweenArrays(log2_spectra, method="quantile")

um <- umap(t(normalized))

# Colours for umap
umdf <- um%>%
  data.frame()%>%
  rownames_to_column("Sample")%>%
  left_join(md, by = "Sample")

um_plot <- ggplot(data = umdf, aes(x = X1, y = X2, colour = Infection_light)) + 
  geom_point(size = 2) +
  stat_ellipse(aes(fill = Infection_light), geom = "polygon", alpha = 0.2, colour = NA) +
  labs(x = "UMAP 1", y = "UMAP 2") + 
  blank_theme +
  theme(
    legend.key.size = unit(3, 'mm'),
    aspect.ratio = 1
  ) +
  scale_colour_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  labs(colour = "Infection/light", fill = "Infection/light", title = "All metabolites")

ggsave(plot = um_plot, filename = file.path(outdir, "plots", "Annotated UMAP.pdf"), 
       width = 100, height = 85, units = "mm")

pdf(file.path(outdir, "plots", "Unnormalised_spectra.pdf"), width = 7, height = 5)
par(mar = c(7, 4.1, 4.1, 2.1))
# Check distributions of samples using boxplots
boxplot(log2_spectra, xlab="", ylab="Log2 counts",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log2_spectra),col="blue")
title("Boxplots of log intensities (unnormalised)")
dev.off()

pdf(file.path(outdir, "plots", "Normalised_spectra.pdf"), width = 7, height = 5)
par(mar = c(7, 4.1, 4.1, 2.1))
# Check distributions of samples using boxplots
boxplot(normalized, xlab="", ylab="Log2 counts",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(normalized),col="blue")
title("Boxplots of log intensities (normalised)")
dev.off()

# Fit the model
fit <- lmFit(normalized, design = design)

colnames(design)

# Make a contrast matrix
contrast.matrix <- makeContrasts(
  Inf_Const_light_vs_Inf_Typical_light  = Inf_Const_light - Inf_Typical_light,
  Uni_Const_light_vs_Uni_Typical_light = Uni_Const_light - Uni_Typical_light,
  Inf_vs_uni = (Inf_Const_light+Inf_Typical_light)/2 - (Uni_Const_light + Uni_Typical_light)/2,
  Const_light = (Inf_Const_light+Uni_Const_light)/2 - (Inf_Typical_light + Uni_Typical_light)/2,
  levels = colnames(design))

fit2 <- contrasts.fit(fit, contrast.matrix)

# Use limma-trend to allow an "intensity-dependent trend for the prior variance"
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# Plot residual standard deviation versus average log 
# expression for a fitted microarray linear model.

pdf(file.path(outdir, "plots", "Limma SA.pdf"), width = 7, height = 5)
# Check distributions of samples using boxplots
plotSA(fit3)
# Let's add a blue horizontal line that corresponds to the median logCPM
title("Limma SA (mean vs residual variance)")
dev.off()

summa.fit <- decideTests(fit3)
summary(summa.fit)

de_summary <- data.frame(summary(summa.fit), check.names = FALSE)%>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3)%>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig")))%>%
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(file.path(outdir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(`Direction if significant`,-`Number of genes`)%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")

# FIX: Create directories natively
dir.create(file.path(outdir, "glimma", "mds"), recursive = TRUE, showWarnings = FALSE)
mds_save <- file.path(outdir, "glimma", "mds", "MDS.html")

htmlwidgets::saveWidget(glimmaMDS(normalized, groups = md$Infection_light,
                                  labels = md), mds_save)

# FIX: Create directories natively
dir.create(file.path(outdir, "toptables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "glimma", "MA"), recursive = TRUE, showWarnings = FALSE)

for(contrast in colnames(contrast.matrix)){
  
  toptable <- topTable(fit3,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(file.path(outdir, "toptables", paste0(contrast, "_toptable.csv")))
  
  vol_save <- file.path(outdir, "glimma", "volcano", paste0(contrast, "_Volcano.html"))
  
  htmlwidgets::saveWidget(glimmaVolcano(fit3, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = normalized,transform.counts	= "none",
                                        dge = fit3, groups = md$Infection_light,
                                        xlab = "log2_normalised_intensity_FC"), vol_save)
  
  ma_save <- file.path(outdir, "glimma", "MA", paste0(contrast, "_MA.html"))
  
  htmlwidgets::saveWidget(glimmaMA(fit3, coef = contrast,main = gsub("_"," ",contrast),
                                   counts = normalized,transform.counts	= "none",
                                   dge = fit3, groups = md$Infection_light,
                                   xlab = "log2_normalised_intensity"), ma_save)
  
}

# FIX: Use native unlink instead of system("rm -r ...")
files_to_remove <- list.files(file.path(outdir, "glimma"), pattern = "_files$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
unlink(files_to_remove, recursive = TRUE)

# Compile the toptables
all_toptables <- list.files(file.path(outdir, "toptables"), full.names = TRUE)

tt_list <- list()
for(i in seq_along(all_toptables)){
  
  # FIX: Properly remove the entire suffix so contrast names remain clean
  contrast_name <- gsub("_toptable\\.csv$", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast_name)
  
  tt_list[[i]] <- tt
  
}

# Compile toptables and save the significant results
if(length(tt_list) > 0) {
  toptables_compiled <- bind_rows(tt_list)
  
  toptables_signif <- toptables_compiled %>%
    #filter(adj.P.Val <= 0.05)%>%
    arrange(adj.P.Val)%>%
    write_csv(file.path(outdir, "Compiled_toptables_all_bile_acids.csv"))
}

# ---------------------------------------------------------------------------
# Session Info Export
# ---------------------------------------------------------------------------
dir.create(file.path(".", "Code_for_submission", "sessioninfo"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(".", "sessioninfo", "Session_info_bile_acids.txt"))
