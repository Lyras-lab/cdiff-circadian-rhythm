#!/bin/bash
# Logging the commands that I ran

# All runs in current working directory

# On the command line
# conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
# Looks like amplicon should have the plugins I need
conda activate qiime2-amplicon-2024.10

# Tabuluate run 1
qiime metadata tabulate \
  --m-input-file sample-metadata-run-1.tsv \
  --o-visualization tabulated-sample-metadata-run-1.qzv

# Tabulate run 2
qiime metadata tabulate \
  --m-input-file sample-metadata.tsv \
  --o-visualization tabulated-sample-metadata.qzv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33
  
# I have paired-end-demux_run_1.qza
qiime demux summarize \
  --i-data paired-end-demux_run_1.qza \
  --o-visualization paired-end-demux_run_1.qzv

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv
  
# Check for consistent start sequences
zcat T282B2d0_1.fastq.gz | awk 'NR%4==2' | cut -c 1-20 | sort | uniq -c | sort -nr | head -20
zcat T282B2d0_2.fastq.gz | awk 'NR%4==2' | cut -c 1-20 | sort | uniq -c | sort -nr | head -20  

# Trim off the adapters for run 1
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux_run_1.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences trimmed_demux_run_1.qza \
  --verbose
  
# Trim off the adapters for run 2
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences trimmed_demux.qza \
  --verbose
  
qiime demux summarize \
  --i-data trimmed_demux_run_1.qza \
  --o-visualization trimmed_demux_run_1.qzv
  
qiime demux summarize \
  --i-data trimmed_demux.qza \
  --o-visualization trimmed_demux.qzv
  
# Run with threads to speed up and selected quailty scores
# Based on inspection of trimmed_demux_run_1.qzv quaity boxplots
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trimmed_demux_run_1.qza  \
    --p-trunc-len-f 256 \
    --p-trunc-len-r 186 \
    --p-n-threads 20 \
    --o-representative-sequences representative-sequences_run_1.qza \
    --o-table feature-table-0_run_1.qza \
    --o-denoising-stats denoising-stats_run_1.qza  

# Run with threads to speed up and selected quailty scores
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trimmed_demux.qza  \
    --p-trunc-len-f 253 \
    --p-trunc-len-r 183 \
    --p-n-threads 20 \
    --o-representative-sequences representative-sequences.qza \
    --o-table feature-table-0.qza \
    --o-denoising-stats denoising-stats.qza  
    
qiime feature-table merge \
  --i-tables feature-table-0_run_1.qza \
  --i-tables feature-table-0.qza  \
  --o-merged-table merged_table.qza
  
qiime feature-table merge-seqs \
  --i-data representative-sequences_run_1.qza \
  --i-data representative-sequences.qza \
  --o-merged-data merged-representative-sequences.qza

# Summarise the merged feature table  
qiime feature-table summarize \
  --i-table merged_table.qza \
  --o-visualization merged_table.qzv \
  --m-sample-metadata-file sample-metadata-merged.tsv

# Tabulate the merged seqs
qiime feature-table tabulate-seqs \
  --i-data merged-representative-sequences.qza \
  --o-visualization merged-representative-sequences.qzv

# Look at the denoising stats for run 1
qiime metadata tabulate \
  --m-input-file denoising-stats_run_1.qza   \
  --o-visualization dada2-stats-summ_run_1.qzv    

# Look at the denoising stats for run 2
qiime metadata tabulate \
  --m-input-file denoising-stats.qza  \
  --o-visualization dada2-stats-summ.qzv
  
# Assign taxonomy to merged reads
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads merged-representative-sequences.qza \
  --p-n-jobs 20 \
  --o-classification merged-taxonomy.qza

# Assign taxonomy to reads
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads representative-sequences.qza \
  --p-n-jobs 20 \
  --o-classification taxonomy.qza

# Assign taxonomy to reads
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads representative-sequences_run_1.qza \
  --p-n-jobs 20 \
  --o-classification taxonomy_run_1.qza
  
# Make a plot of taxonomy assignment for run1
qiime metadata tabulate \
  --m-input-file taxonomy_run_1.qza \
  --o-visualization taxonomy_run_1.qzv

# Make a plot of taxonomy assignment
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Make a per-sample feature assignment visualisation
qiime feature-table summarize \
  --i-table feature-table-0.qza \
  --o-visualization feature-assignment-table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

# De novo multiple sequence alignment    
qiime alignment mafft \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza
  
# 'Mask (i.e., filter) unconserved and highly gapped columns from an alignment.'
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza  
  
# Build and unrooted tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza

# Midpoint root an unrooted phylogenetic tree. Roughly in the middle of the data  
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Export the tree to a nwk file
qiime tools export \
    --input-path rooted-tree.qza \
    --output-path exported-rooted-tree

# Changed sampling depth to 
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table feature-table-0.qza \
  --p-sampling-depth 8500 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity-core-metrics-phylogenetic

# Taxonomy barplot
qiime taxa barplot \
  --i-table feature-table-0.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
  
# Run ancomBC
# ancombc multi formula with reference levels
qiime composition ancombc \
  --i-table feature-table-0.qza \
  --m-metadata-file sample-metadata.tsv \
    --p-formula 'all_group + subject' \
    --p-reference-levels all_group::Uninfected_Typicallight_0_Yes subject::subject_1 \
    --o-differentials dataloaf.qza

qiime composition da-barplot \
  --i-data dataloaf.qza \
  --p-significance-threshold 0.001 \
  --o-visualization da-barplot-subject.qzv

