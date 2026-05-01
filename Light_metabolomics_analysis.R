library(tidyverse)
library(limma)
library(edgeR)
library(Glimma)
library(KEGGREST)
library(stringr)
library(ComplexHeatmap)
library(uwot)
library(limpa)
library(fgsea)

# Set ggplot2 themes for the paper
blank_theme <- theme_bw(base_size = 7)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7))

# Looking at this for analysis
# https://support.bioconductor.org/p/64484/

condition_colors <- c(
  "Uni_Typical_light" = "#008000",   
  "Uni_Const_light" = "#C000C0",   
  "Inf_Typical_light" = "#0396FE",   
  "Inf_Const_light" = "#FFA040"   
)

outdir <- "./output_metabolomics/"
system(paste0("mkdir -p ",outdir, "/plots/barcode_plots/"))

# Read in the protein data
protein <- read_csv("./Metabolomics_inputs/P21_0421_Exp5_StatsTable.csv")

duplicated_row <- protein[,grepl("Allysine", colnames(protein))]

protein_mat <- as.matrix(protein[,3:ncol(protein)])
rownames(protein_mat) <- protein$sample

protein_mat <- t(protein_mat)

min(protein_mat)

protein_mat_na <- protein_mat
protein_mat_na[protein_mat_na == 0] <- NA

md <- data.frame(Sample = colnames(protein_mat))%>%
  # Drop a sample where diff was not properly detected
  filter(Sample != "Typical_light_Cdifficile_5")%>%
  mutate(Infection = ifelse(grepl("Cdifficile",Sample), "Inf", "Uni"))%>%
  mutate(Light = ifelse(grepl("Typical_light",Sample), "Typical_light", "Const_light"))%>%
  mutate(Infection_light = paste0(Infection, "_", Light))%>%
  mutate(Infection_light = replace(Infection_light, grepl("QC",Sample), "QC"))%>%
  filter(Infection_light != "QC")%>%
  mutate(Infection_light = factor(Infection_light, levels = names(condition_colors)))

# Drop the peaks that don't appear consistently in at least one condition (or 5 samples)
table(md$Infection_light)

protein_mat <- protein_mat[,md$Sample]

colnames(protein_mat) == md$Sample

# Log transfrom the counts
log2_spectra <- log2(protein_mat)
log2_spectra[log2_spectra == -Inf] <- NA

protein.id <- rownames(log2_spectra)

dpcest <- dpcCN(log2_spectra)
y.protein <- dpcQuant(log2_spectra, protein.id, dpc=dpcest)

plotMDSUsingSEs(y.protein)

# Set up a design matrix
design <- model.matrix(~0+ Infection_light, data = md)
colnames(design) <- gsub("Infection_light", "", colnames(design))
rownames(design) <- colnames(protein_mat)

hist(y.protein$E, breaks = 300)

um <- umap(t(y.protein$E))

# Colours for umap
umdf <- um%>%
  data.frame()%>%
  rownames_to_column("Sample")%>%
  left_join(md)

um_plot <- ggplot(data = umdf, aes(x = X1, y = X2, colour = Infection_light)) + 
  geom_point(size = 2) +
  stat_ellipse(aes(colour = Infection_light), geom = "polygon", alpha = 0.2, fill = NA) +
  labs(x = "UMAP 1", y = "UMAP 2") + 
  blank_theme +
  theme(
    legend.key.size = unit(3, 'mm'),
    aspect.ratio = 1
  ) +
  scale_colour_manual(values = condition_colors) +
  labs(colour = "Infection/light", title = "All metabolites")

ggsave(plot = um_plot,paste0(outdir, "/plots/Annotated UMAP.pdf"), 
       width = 100, height = 85, units = "mm")

pdf(paste0(outdir, "/plots/Unnormalised_spectra.pdf"), width = 7, height = 5)
par(mar = c(7, 4.1, 4.1, 2.1))
# Check distributions of samples using boxplots
boxplot(log2_spectra, xlab="", ylab="Log2 counts",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log2_spectra),col="blue")
title("Boxplots of log intensities (unnormalised)")
dev.off()

# Fit the model
fit <- dpcDE(y.protein, design)

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
fit3 <- eBayes(fit2)

# Plot residual standard deviation versus average log 
# expression for a fitted microarray linear model.

pdf(paste0(outdir, "/plots/Limma SA.pdf"), width = 7, height = 5)
# Check distributions of samples using boxplots
plotSA(fit3)
# Let's add a blue horizontal line that corresponds to the median logCPM
title("Limma SA (mean vs residual variance)")
dev.off()

summa.fit <- decideTests(fit3)
summary(summa.fit)

de_summary <- data.frame(summary(summa.fit), check.names = F)%>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3)%>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig")))%>%
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(paste0(outdir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  #filter(`Number of genes` >10)%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(`Direction if significant`,-`Number of genes`)%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")

system(paste0("mkdir -p ", outdir, "/glimma/mds/"))
mds_save <-paste0(outdir,"/glimma/mds/", "MDS.html")

htmlwidgets::saveWidget(glimmaMDS(y.protein$E, groups = md$Infection_light,
                                  labels = md), mds_save)

system(paste0("mkdir -p ", outdir, "/toptables"))
system(paste0("mkdir -p ", outdir, "/glimma/volcano/"))
system(paste0("mkdir -p ", outdir, "/glimma/MA/"))

for(contrast in colnames(contrast.matrix)){
  
  toptable <- topTable(fit3,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("Gene")%>%
    write_csv(paste0(outdir,"toptables/", contrast, "_toptable.csv"))
  
  vol_save <- paste0(outdir,"/glimma/volcano/",contrast, "_Volcano.html")
  
  htmlwidgets::saveWidget(glimmaVolcano(fit3, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = y.protein$E,transform.counts	= "none",
                                        dge = fit3, groups = md$Infection_light,
                                        xlab = "log2_normalised_intensity_FC"), vol_save)
  
  ma_save <- paste0(outdir,"/glimma/MA/",contrast, "_MA.html")
  
  htmlwidgets::saveWidget(glimmaMA(fit3, coef = contrast,main = gsub("_"," ",contrast),
                                   counts = y.protein$E,transform.counts	= "none",
                                   dge = fit3, groups = md$Infection_light,
                                   xlab = "log2_normalised_intensity"), ma_save)
  
}

system(paste0("rm -r ", outdir, "glimma/*/*_files"))

# Compile the toptables
all_toptables <- list.files(paste0(outdir, "toptables/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptables)){
  
  contrast <- gsub("_toptable.csv", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
  
}

# Compile toptables and save the significant results
toptables_compiled <- bind_rows(tt_list)

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val <= 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(outdir, "Compiled_toptables_significant_genes.csv"))

all_metabolites <- rownames(protein_mat)

duplicated_names <- gsub("\\.\\.\\..*", "", rownames(protein_mat))
duplicated_names[duplicated(duplicated_names)]

unique_names <- unique(duplicated_names)

sum(duplicated_names == rownames(y.protein$E))

length(unique_names)

unique_df <- data.frame(unique_names)

# Save the metabolite sets
df_metabolites <- data.frame(Metabolite_Name = unique_names)

# listDatabases()
# org <- keggList("organism")
# org_df <- data.frame(org)
# 
# # Read in the kegg compound database
# compound_all <- read_tsv("/oldvol/apattison/Data/Lyras_lab/Desirel/cdiff_light_protein/ref/compound.txt", col_names = c("kegg_id","compound"))
# 
# comp_list <- list()
# # Loop over each compound and get the kegg database entry
# for(i in 1:length(unique_names)){
#   
#   query <- unique_names[i]
#   
#   if(nchar(query) == 0){
#     
#     next()
#     
#   }
#   
#   if(grepl("\\[", query)){
#     next()
#   }
#   
#   print(query)
#   
#   # Looks like everything is in the compound database
#   compound <- compound_all%>%
#     filter(grepl(query, compound))
#   
#   # Skip this interation if nothing is found
#   if(length(compound) == 0){
#     
#     next()
#     
#   }
#   
#   cdf <- data.frame(compound)%>%
#     filter (grepl(paste0("^",query,"\\;"), compound) | 
#               grepl(paste0("\\; ", query, "\\;"), compound) |
#               grepl(paste0("\\; ", query, "$"), compound)|
#               grepl(paste0("^",query,"$"), compound) )%>%
#     mutate(query = query)
#   
#   if(nrow(cdf) == 0){
#     
#     cdf <- data.frame(compound)%>%
#       filter (grepl(paste0("^D-",query,"\\;"), compound) | 
#                 grepl(paste0("\\; D-", query, "\\;"), compound) |
#                 grepl(paste0("\\; D-", query, "$"), compound)|
#                 grepl(paste0("^D-",query,"$"), compound)|
#                 grepl(paste0("^L-",query,"\\;"), compound) | 
#                 grepl(paste0("\\; L-", query, "\\;"), compound) |
#                 grepl(paste0("\\; L-", query, "$"), compound)|
#                 grepl(paste0("^L-",query,"$"), compound))%>%
#       mutate(query = query)
#     
#   }
#   
#   if(nrow(cdf) == 0){
#     
#     next()
#     
#   }
#   
#   
#   print(nrow(cdf))
#   
#   comp_list[[i]] <- cdf
#   
# }
# 
# all_comps <- bind_rows(comp_list)%>%
#   group_by(query)%>%
#   mutate(Count = n())%>%
#   ungroup()
# 
# # Remove any remaining duplicates
# comps_dups_removed <- all_comps%>%
#   filter(!duplicated(query))
# 
# brite_list <- list() 
# # Loop oer the kegg IDs and get the brite pathays
# for(i in 1:nrow(comps_dups_removed)){
#   
#   kegg_id <- comps_dups_removed$kegg_id[i]
#   
#   print(kegg_id)
#   
#   kegg_results <- tryCatch(
#     {keggGet(kegg_id)},
#     error = function(cond) {
#       message(conditionMessage(cond))
#       # Choose a return value in case of error
#       NULL
#     })
#   if(is.null(kegg_results)){
#     next
#   }
#   
#   if(length(kegg_results)>0){
#     
#     brite <- kegg_results[[1]]$BRITE
#     
#     if(!is.null(brite)){
#     
#     brite_df <- data.frame(brite, kegg_id,compound =  comps_dups_removed$compound[i], query = comps_dups_removed$query[i])
#     
#     brite_list[[i]] <- brite_df
#     }
#   }
#   
#   # Wait 0.2 seconds as the server didn't like all my requests
#   Sys.sleep(0.2)
# }
# 
# all_brite <- bind_rows(brite_list)%>%
#   mutate(brite = str_trim(brite))%>%
#   group_by(brite)%>%
#   mutate(count = n())
# 
# # Make gene set lists
# metabolite_list <- list()
# for(i in 1:length(unique(all_brite$brite))){
#   
#   metabolite_set <- unique(all_brite$brite)[i]
#   
#   metabolites <- all_brite%>%
#     filter(brite == metabolite_set)
#   
#   metabolite_list[[metabolite_set]] <- metabolites$query
#   
# }
# 
# head(metabolite_list)
# 
# saveRDS(metabolite_list, paste0(outdir,"/ref/metabolites_for_GSEA.rds"))

metabolite_list <- readRDS("./Metabolomics_inputs/metabolites_for_GSEA.rds")
# 
# metabolo_list_df <- do.call(rbind, lapply(names(metabolite_list), function(n) {
#   data.frame(Metabolite_set = n, Metabolite = metabolite_list[[n]], stringsAsFactors = FALSE)
# }))
# 
# write_csv(metabolo_list_df, paste0(outdir,"/ref/metabolites_for_GSEA_data_frame.csv"))

metabolo_list_df <- read_csv("./Metabolomics_inputs/metabolites_for_GSEA_data_frame.csv")

# Use the duplicated names for identification
indexed <- ids2indices(gene.sets = metabolite_list, identifiers = duplicated_names, remove.empty=TRUE)
# Make a directory
system(paste0("mkdir -p ", outdir,"gsea/camera/"))
system(paste0("mkdir -p ", outdir,"gsea/fry/"))

for(contrast in colnames(contrast.matrix)){
  
  camera_result <- camera(y = y.protein ,index = indexed, design = design,
                          contrast = contrast.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    filter(FDR < 0.1)%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    mutate(Contrast= contrast)
  
  fry_result <- fry(y = y.protein ,index = indexed, design = design,
                    contrast = contrast.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    filter(FDR < 0.05)%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    mutate(Contrast= contrast)
  
  write_csv(camera_result, paste0(outdir,"gsea/camera/","Metabolites", "_", contrast,".csv"))
  write_csv(fry_result, paste0(outdir,"gsea/fry/","Metabolites", "_", contrast,".csv"))
  
}

# Compile the camera results
all_camera <- list.files(paste0(outdir,"gsea/camera/"), full.names = T)

clist <- list()
for(i in 1:length(all_camera)){
  
  contrast <- gsub("\\.csv", "", basename(all_camera[i]))
  
  tt <- read_csv(all_camera[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    dplyr::select(-Contrast)
  
  clist[[i]] <- tt
  
}

# Compile toptables and save the significant results
camera_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  write_csv(paste0(outdir, "Compiled_significant_gene_sets_camera.csv"))

# Compile the fry results
all_fry <- list.files(paste0(outdir,"gsea/fry/"), full.names = T)

clist <- list()
for(i in 1:length(all_fry)){
  
  contrast <- gsub("\\.csv", "", basename(all_fry[i]))
  
  tt <- read_csv(all_fry[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    dplyr::select(-Contrast)
  
  clist[[i]] <- tt
  
}

# Compile toptables and save the significant results
fry_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(-FDR)%>%
  mutate(NGenes = as.numeric(NGenes))%>%
  filter(NGenes>1)%>%
  mutate(`Gene set` = factor(`Gene set`, levels = unique(`Gene set`)))%>%
  write_csv(paste0(outdir, "Compiled_significant_gene_sets_fry.csv"))

fry_compiled_plot <- fry_compiled%>%
  filter(contrast == "Metabolites_Inf_Const_light_vs_Inf_Typical_light")

inf_light <- toptables_compiled%>%
  filter(contrast == "Inf_Const_light_vs_Inf_Typical_light")

DE <- inf_light%>%
  arrange(t)

ranks <- as.numeric(DE$t)
names(ranks) <- DE$Gene

head(ranks)

# Run FGSEA with our data
# Gene set list is made up from the eggonog data
fgseaRes <- fgsea(pathways = metabolite_list, 
                  stats    = ranks,
                  minSize  = 2,
                  maxSize  = 2000)

# Keep just the signifcant pathways
fgseaRes_ordered <- fgseaRes%>%
  arrange(padj)
  
fry_gene_sets <- ggplot(fry_compiled_plot, aes(x = -log10(FDR), y = `Gene set`, size = NGenes, color = Direction)) +
  geom_point(alpha = 0.8) +
  geom_segment(aes(x = 0, xend = -log10(FDR), y = `Gene set`, yend = `Gene set`, color = Direction), size = 0.8) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +  # Red for up, blue for down
  scale_radius(breaks = c(2,4, 8, max(fry_compiled_plot$NGenes)))+
  labs(
    x = expression(-log[10]("FDR")),
    y = "Gene Set",
    size = "Gene Count",
    color = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank()
  )+
  ggtitle("Inf const light vs\ninf typical light (FDR < 0.05)")

fry_gene_sets

ggsave(plot = fry_gene_sets, filename = paste0(outdir, "/plots/", "Fry_gene_sets_Inf_Const_light_vs_Inf_Typical_light.pdf"),width = 10, height = 7)

# Make some barcode plots of some key gene sets
colnames(contrast.matrix)

# Make a barcodeplot for all gene sets that were significant for Inf_Const_light_vs_Inf_Typical_light
for(gene_set in unique(c(as.character(fry_compiled$`Gene set`), as.character(camera_compiled$`Gene set`)))){
  
  print(gene_set)
  
  save_name <- paste0(outdir, "plots/barcode_plots/Inf_Const_light_vs_Inf_Typical_light_", 
                      gene_set ,".pdf")
  
  index <- indexed[[gene_set]]
  
  pdf(save_name, width = 7, height = 4)
  barcodeplot(fit3$t[,"Inf_Const_light_vs_Inf_Typical_light"],
              index=index,
              labels = c("Down","Up"),
              main= paste0(gene_set, "\nN metabolites=", length(index)),
              xlab = "Limma t statistic: Inf const light vs inf typical light")
  dev.off()
  
}

gene_set_heatmap_metabolo <- function(index, desc, small = F){
  
  md_orderd <- md %>%
    arrange(Infection_light)
  
  # Grab the genes from the CPM
  genes <- y.protein$E[index,md_orderd$Sample]%>%
    t()%>%
    scale()%>%
    t()
  
  ha = HeatmapAnnotation(Condition = md_orderd$Infection_light,
                         col = list(Condition=condition_colors))
 
  if(small == T){
    
    Heatmap(genes,column_title = desc, name = "Z score", row_names_side = "left",top_annotation = ha,
            row_names_gp = gpar(fontsize = 5),
            column_title_gp = gpar(fontsize = 7),
            row_names_max_width = unit(30, "cm"), cluster_columns = F, column_split = md_orderd$Infection_light)
    
  }else{
    Heatmap(genes,column_title = desc, name = "Z score", row_names_side = "left",top_annotation = ha,
            cluster_columns = F, column_split = md_orderd$Infection_light)
  }

}

gene_set_heatmap_metabolo_inf <- function(index, desc, small = F){
  
  md_orderd <- md %>%
    arrange(Infection_light)%>%
    filter(Infection =="Inf")%>%
    mutate(Infection_light = droplevels(Infection_light))
  
  # Grab the genes from the CPM
  genes <- y.protein$E[index,md_orderd$Sample]%>%
    t()%>%
    scale()%>%
    t()
  
  ha = HeatmapAnnotation(Condition = md_orderd$Infection_light,
                         col = list(Condition=condition_colors))
  
  if(small == T){
    
    Heatmap(genes,column_title = desc, name = "Z score", row_names_side = "left",top_annotation = ha,
            row_names_gp = gpar(fontsize = 5),
            column_title_gp = gpar(fontsize = 7),
            row_names_max_width = unit(30, "cm"), cluster_columns = F, column_split = md_orderd$Infection_light)
    
  }else{
    Heatmap(genes,column_title = desc, name = "Z score", row_names_side = "left",top_annotation = ha,
            cluster_columns = F, column_split = md_orderd$Infection_light)
  }
  
}


system(paste0("mkdir -p ", outdir,"plots/heatmaps/"))
system(paste0("mkdir -p ", outdir,"plots/heatmaps_inf/"))

# Make for loop to run the heatmap function for all signifcant gene sets
for(gene_set in unique(fry_compiled$`Gene set`)){
  
  print(gene_set)
  
  save_name <- paste0(outdir, "plots/heatmaps/", gene_set ,".pdf")
  
  pdf(save_name, width = 14, height = 4 + (0.5*length(indexed[[gene_set]])))
  hm <- gene_set_heatmap_metabolo(indexed[[gene_set]] , desc =  gene_set)
  draw(hm) 
  dev.off()
  
}

for(gene_set in unique(fry_compiled$`Gene set`)){
  
  print(gene_set)
  
  save_name <- paste0(outdir, "plots/heatmaps_inf/", gene_set ,".pdf")
  
  pdf(save_name, width = 14, height = 4 + (0.5*length(indexed[[gene_set]])))
  hm <- gene_set_heatmap_metabolo_inf(indexed[[gene_set]] , desc =  gene_set)
  draw(hm) 
  dev.off()
  
}

# Save the R session info for methods section
writeLines(capture.output(sessionInfo()), paste0("./sessioninfo/Session_info_metabolomics.txt"))
