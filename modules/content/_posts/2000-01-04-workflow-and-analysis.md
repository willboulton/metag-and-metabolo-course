---
title: "Workflow and Analysis: Running the Complete Pipeline"
---

## Overview

This module walks through the complete analysis pipeline from raw data to biological insights. We'll use the DIME dataset to demonstrate each step, but the workflow is generalizable to other studies.

---

## Part 1: Metagenomics Analysis Pipeline

### Step 1: Quality Control and Preprocessing

Before assembly, we need to ensure our reads are high quality.

**Why this matters:**
- Low-quality reads introduce errors into assemblies
- Adapters from sequencing can confuse assembly algorithms
- Contamination should be detected and removed

**Tools:** FastQC, MultiQC for visualization; Trimmomatic or fastp for trimming

**Output:** Clean FASTQ files ready for assembly

### Step 2: Assembly

Assembly reconstructs the genome(s) in your sample from overlapping short reads.

**MATAFILER4 assembler:** Typically uses metaSPAdes or MEGAHIT
- **Input:** High-quality reads
- **Output:** Contigs (assembled sequences, typically 500 bp to several Mb)
- **Key parameter:** k-mer size (affects trade-off between coverage and accuracy)

**What to expect:**
- Assembly stats report: N50, max contig length, number of contigs
- Typically 50-80% of reads assemble into contigs for well-characterized environments
- Assembly quality varies by sample complexity and sequencing depth

### Step 3: Abundance Estimation

Determine how many reads from each sample map to each contig.

**Process:**
1. Map all reads (from all samples) back to assembled contigs
2. Count overlapping reads to estimate coverage
3. Normalize by contig length to get abundance per sample

**Output:** Abundance matrix (contigs × samples with coverage/RPKM values)

### Step 4: Binning and MAG Recovery

Binning groups contigs into genomic bins (potential genomes) based on:
- Coverage patterns across samples (similar contigs should have similar coverage)
- GC content (organism-specific)
- Tetranucleotide frequency (4bp composition bias by organism)

**Binning methods:** MaxBin2, MetaBAT2, CONCOCT (often run in ensemble)

**Output:** Individual genome bins (.fa files), one per MAG

**Quality assessment:**
- Use tools like CheckM to evaluate **completeness** (% of essential genes present)
- Assess **contamination** (% of genes suggesting multiple organisms)
- Filter: typically keep MAGs with >50% completeness and <5-10% contamination

### Step 5: Dereplication

When analyzing multiple samples, you may recover the same genome multiple times (as different bins). Dereplication creates one representative per unique genome.

**Tools:** dRep or similar
- **Input:** All quality MAGs from all samples
- **Output:** Non-redundant set of genomes + mapping of original bins to representatives

### Step 6: Taxonomy Assignment

Assign taxonomic identity to each MAG.

**Methods:**
- **Gene-based:** Use gene annotations to infer taxonomy (if >50% of genes match a taxon)
- **Genome-based:** Compare MAG to reference genomes using tools like GTDB-tk

**Output:** Taxonomy for each MAG (Phylum, Class, Order, Family, Genus, Species level when possible)

### Step 7: Functional Annotation

Identify genes in MAGs and annotate their functions.

**Process:**
1. Gene prediction: Use prodigal or similar to identify open reading frames
2. Functional assignment: Compare genes to databases (KEGG, COG, Pfam)
3. Pathway mapping: Assign genes to metabolic pathways

**Outputs:**
- Gene inventory for each MAG
- Annotation to functional categories
- Predicted metabolic pathways (can reveal "what this organism can do")

### Step 8: Create Analysis-Ready Tables

Convert MATAFILER4 outputs into tables suitable for R analysis.

**Metagenomics tables needed:**
1. **Abundance table**: MAGs × samples (relative or absolute abundance)
2. **Taxonomy table**: MAG IDs with full taxonomy (Phylum - Species)
3. **Functional table** (optional): MAG × KEGG pathways (presence/absence or gene count)

**Example format (abundance_table.csv):**
```
MAG_ID,Sample_1,Sample_2,Sample_3,...
MAG_000001,0.05,0.03,0.08,...
MAG_000002,0.12,0.15,0.10,...
...
```

---

## Part 2: Metabolomics Analysis Pipeline

### Step 1: Feature Detection and Preprocessing

Starting point: Raw mass spectrometry files (.raw, .mzML, etc.)

**Preprocessing includes:**
- Peak picking (identify discrete m/z signals)
- Baseline correction
- Noise filtering

**Tools:** XCMS, MSConvert, or vendor-specific software

**Output:** Feature matrix (m/z × retention time × samples) with intensity values

### Step 2: Normalization

Adjust for technical variation in peak intensities.

**Common approaches:**
- **Total intensity normalization:** Divide each sample by total peak intensity
- **Log-transformation:** Account for variance being proportional to mean
- **Batch correction:** Account for technical variation between instrument runs

**Why normalize:**
- Accounts for differences in injection volume, instrument sensitivity, etc.
- Makes features comparable across samples
- Improves statistical power

### Step 3: Feature Annotation

Attempt to identify which metabolites your peaks represent.

**Approaches:**
- **Database matching:** Compare m/z and retention time to compound databases
- **In-house standards:** Compare to known compounds you've measured
- **Elemental composition prediction:** Narrow down possibilities based on accurate mass

**Reality check:**
- ~10-30% of features are typically identified in untargeted analysis
- Many features remain as "unknown features" (Feature_m/z_12.3456_RT_4.5)
- Can use KEGG ID if available to connect to metabolic pathways

**Example metadata (features_metadata.csv):**
```
Feature_ID,m/z,RT_min,Compound_ID,KEGG_ID,Annotation_Source
Feature_001,104.107,2.3,L-Aspartate,C00049,KEGG_database
Feature_002,118.088,3.1,Unknown,NA,Unknown
...
```

### Step 4: Quality Filtering

Remove low-quality features.

**Criteria:**
- Remove features with >50% missing values across samples
- Remove features with very low abundance (below detection limit in most samples)
- Remove likely artifacts/contaminants (features with extreme m/z or chemistry)

**Output:** Filtered feature matrix ready for analysis

### Step 5: Create Analysis-Ready Tables

Convert preprocessed data into R-compatible format.

**Metabolomics tables needed:**
1. **Features table**: Features × samples (peak intensities, normalized)
2. **Features metadata**: Feature ID, m/z, RT, compound name, KEGG ID, etc.
3. **Samples metadata**: (same as metagenomics)

---

## Part 3: Integrated Analysis in R

### Load Data into R

```r
# Load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(igraph)
library(SpiecEasi)

# Load metagenomics data
mag_abundance <- read.csv("data/metagenomics/abundance_tables.csv", 
                          row.names = 1)
mag_taxonomy <- read.csv("data/metagenomics/taxonomy.csv", row.names = 1)

# Load metabolomics data
metabolite_abundance <- read.csv("data/metabolomics/features.csv", 
                                 row.names = 1)
metabolite_metadata <- read.csv("data/metabolomics/features_metadata.csv", 
                                row.names = 1)

# Load sample metadata (shared between both)
sample_metadata <- read.csv("data/metadata/sample_info.csv", row.names = 1)
```

### Step 1: Separate Analyses

First, analyse each data type independently to understand the structure of each dataset before combining them.

#### **Metagenomics diversity analysis:**

```r
library(vegan)

# Alpha diversity (within-sample diversity)
# Operates on sample rows — transpose so samples are rows
alpha_div <- data.frame(
  Shannon = diversity(t(mag_abundance), index = "shannon"),
  Simpson = diversity(t(mag_abundance), index = "simpson"),
  Richness = specnumber(t(mag_abundance))
)
alpha_div$treatment <- sample_metadata[rownames(alpha_div), "treatment"]

# Test for treatment differences in Shannon diversity
kruskal.test(Shannon ~ treatment, data = alpha_div)

# Beta diversity — Bray-Curtis distance for ordination
bray_curtis_dist <- vegdist(t(mag_abundance), method = "bray")
pcoa_result <- cmdscale(bray_curtis_dist, k = 2, eig = TRUE)

# Statistical test: PERMANOVA (note: use adonis2, not the deprecated adonis)
permanova_mag <- adonis2(bray_curtis_dist ~ treatment,
                         data = sample_metadata,
                         permutations = 999)
print(permanova_mag)

# Visualise PCoA
plot(pcoa_result$points,
     col = as.numeric(factor(sample_metadata$treatment)),
     pch = 19,
     xlab = paste0("PCoA1 (", round(pcoa_result$eig[1]/sum(pcoa_result$eig)*100, 1), "%)"),
     ylab = paste0("PCoA2 (", round(pcoa_result$eig[2]/sum(pcoa_result$eig)*100, 1), "%)"),
     main = "MAG community composition (Bray-Curtis PCoA)")
legend("topright", legend = levels(factor(sample_metadata$treatment)),
       col = 1:nlevels(factor(sample_metadata$treatment)), pch = 19)
```

> **Note on `adonis2` vs `adonis`:** The original `adonis()` function in vegan is deprecated as of vegan 2.6. Use `adonis2()`, which has a slightly different interface (formula on the left-hand side) but equivalent functionality.

#### **Metabolomics exploration:**

```r
# CLR transformation before PCA (addresses compositionality)
# Add pseudocount to handle zeros
clr_transform <- function(mat, pseudocount = 0.5) {
  mat_pseudo <- mat + pseudocount
  log(mat_pseudo / apply(mat_pseudo, 2, function(x) exp(mean(log(x)))))
}

metabolite_clr <- clr_transform(metabolite_abundance)

# PCA of CLR-transformed metabolite profiles
pca_result <- prcomp(t(metabolite_clr), scale. = FALSE)  # CLR already standardises
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

plot(pca_result$x[, 1:2],
     col = as.numeric(factor(sample_metadata$treatment)),
     pch = 19,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "Metabolite profiles (CLR-PCA)")

# PERMANOVA on metabolite Aitchison distance (= Euclidean distance in CLR space)
aitchison_dist <- dist(t(metabolite_clr), method = "euclidean")
permanova_met <- adonis2(aitchison_dist ~ treatment,
                         data = sample_metadata,
                         permutations = 999)
print(permanova_met)
```

### Step 2: Network-Based Co-analysis with SpiecEasi

SpiecEasi builds a network where nodes are MAGs and metabolites, and edges represent significant partial correlations — indicating direct associations after accounting for all other variables. This is more robust than pairwise Pearson correlations for sparse, compositional data (see the Theory section for the conceptual background).

#### **Prepare combined data:**

SpiecEasi expects a **samples × features** matrix (samples as rows, features as columns) with **count or relative abundance** data — it applies its own log-ratio transformation internally. Filter to features present in a minimum number of samples to reduce noise.

```r
library(SpiecEasi)

# Transpose so rows = samples, columns = features
mag_t <- t(mag_abundance)          # samples × MAGs
met_t <- t(metabolite_abundance)   # samples × metabolites

# Filter features: keep those detected in at least 20% of samples
min_prevalence <- 0.2
mag_keep <- colMeans(mag_t > 0) >= min_prevalence
met_keep <- colMeans(met_t > 0) >= min_prevalence
mag_filt <- mag_t[, mag_keep]
met_filt <- met_t[, met_keep]

cat(sprintf("Retained %d MAGs and %d metabolite features for network analysis\n",
            ncol(mag_filt), ncol(met_filt)))

# Ensure same sample order
stopifnot(all(rownames(mag_filt) == rownames(met_filt)))
```

#### **Run SpiecEasi (multi-domain mode):**

For cross-domain inference (MAGs and metabolites together), use `spiec.easi()` with the `data2` argument — this fits a single network across both tables while respecting that they come from different compositional domains.

```r
# Run SpiecEasi in MB (Meinshausen-Bühlmann) mode for cross-domain inference
# lambda.min.ratio and nlambda control regularisation — smaller ratio = denser network
# StARS selects lambda based on network stability across subsamples
set.seed(42)
se_result <- spiec.easi(
  list(mag_filt, met_filt),
  method      = "mb",
  lambda.min.ratio = 1e-2,
  nlambda     = 20,
  pulsar.params = list(rep.num = 20, seed = 42)
)

# Extract the adjacency matrix (binary: edge present/absent)
net_matrix <- as.matrix(getRefit(se_result))
rownames(net_matrix) <- c(colnames(mag_filt), colnames(met_filt))
colnames(net_matrix) <- c(colnames(mag_filt), colnames(met_filt))

cat(sprintf("Network has %d nodes and %d edges\n",
            nrow(net_matrix), sum(net_matrix) / 2))
```

#### **Visualise network:**

```r
library(igraph)

# Build igraph object
g <- graph_from_adjacency_matrix(net_matrix, mode = "undirected", diag = FALSE)

# Label nodes by type
n_mags <- ncol(mag_filt)
node_type <- c(rep("MAG", n_mags), rep("Metabolite", ncol(met_filt)))
V(g)$type   <- node_type
V(g)$color  <- ifelse(node_type == "MAG", "#6BAED6", "#FC8D59")
V(g)$size   <- ifelse(node_type == "MAG", 6, 4)
V(g)$label  <- NA  # hide labels for dense networks

# Use a force-directed layout
set.seed(42)
layout_fr <- layout_with_fr(g)

plot(g,
     layout     = layout_fr,
     vertex.color = V(g)$color,
     vertex.size  = V(g)$size,
     edge.color   = "grey70",
     edge.width   = 0.5,
     main = "MAG–Metabolite co-occurrence network (SpiecEasi MB)")
legend("bottomleft",
       legend = c("MAG", "Metabolite"),
       fill   = c("#6BAED6", "#FC8D59"),
       bty    = "n")
```

#### **Extract associations of interest:**

```r
# Find all metabolites connected to a specific taxon
# (e.g., a Faecalibacterium MAG of interest)
target_mag <- "MAG_000042"   # replace with your MAG ID
neighbours <- names(neighbors(g, v = target_mag))
metabolite_neighbours <- neighbours[neighbours %in% colnames(met_filt)]
cat("Metabolites correlated with", target_mag, ":\n")
print(metabolite_neighbours)
```

### Step 3: Statistical Association Testing

Test which organisms and metabolites are significantly associated with treatment.

#### **Differential abundance of MAGs:**

```r
library(DESeq2)

# Round abundance to integers for DESeq2 (it expects count-like data)
mag_counts <- round(mag_abundance * 1e6)  # scale to pseudo-counts

# Build DESeq2 object
coldata <- data.frame(
  treatment = factor(sample_metadata$treatment),
  row.names = colnames(mag_counts)
)
dds <- DESeqDataSetFromMatrix(countData = mag_counts,
                              colData   = coldata,
                              design    = ~ treatment)
dds <- DESeq(dds)
res_mag <- results(dds, contrast = c("treatment", "intervention", "control"))

# Summarise significant MAGs
sig_mags <- subset(res_mag, padj < 0.05 & abs(log2FoldChange) > 1)
cat(sprintf("%d MAGs significantly changed (FDR < 0.05, |log2FC| > 1)\n",
            nrow(sig_mags)))
```

#### **Differential abundance of metabolites:**

```r
# Log-transform and test with linear models (paired design if longitudinal)
metabolite_log <- log2(metabolite_abundance + 0.5)

# Simple t-test (or use limma for better variance estimation)
library(limma)
design_mat <- model.matrix(~ treatment, data = sample_metadata)
fit <- lmFit(metabolite_log, design = design_mat)
fit <- eBayes(fit)
top_metabolites <- topTable(fit, coef = "treatmentintervention",
                            number = Inf, adjust.method = "BH")
sig_mets <- subset(top_metabolites, adj.P.Val < 0.05 & abs(logFC) > 0.5)
cat(sprintf("%d metabolite features significantly changed\n", nrow(sig_mets)))
```

#### **Using DDT for automated analysis:**

DDT runs the above analyses (and more) in a single function call and generates a standardised set of output figures and tables. This is particularly useful for exploring the full dataset systematically before diving into specific hypotheses.

```r
library(DDT)

# Build a phyloseq object for DDT input
# (see phyloseq documentation for full details)
OTU  <- otu_table(mag_abundance, taxa_are_rows = TRUE)
TAX  <- tax_table(as.matrix(mag_taxonomy))
META <- sample_data(sample_metadata)
physeq <- phyloseq(OTU, TAX, META)

# Run DDT — generates alpha diversity plots, PCoA, differential abundance
ddt_results <- run_ddt(
  physeq     = physeq,
  group_var  = "treatment",
  output_dir = "results/DDT/metagenomics/",
  contrast   = c("intervention", "control")
)

# DDT returns a list of result tables and figure paths
str(ddt_results, max.level = 1)
```

### Step 4: Biological Interpretation

Connect analysis results back to biological hypotheses.

**For DIME specifically:**
- Which MAGs correlate with polyphenol metabolites? (potential polyphenol degraders)
- Do SCFA-producing taxa correlate with SCFA metabolites?
- How do microbial community shifts predict metabolic changes?

**Interpretation examples:**
```r
# Find MAGs positively correlated with butyrate (if annotated)
butyrate_idx <- which(rownames(metabolite_metadata) == "Butyrate")
butyrate_correlations <- net[butyrate_idx, ]
butyrate_partners <- names(which(butyrate_correlations != 0))

# Are these known SCFA producers? Check their annotations...
```

---

## Part 4: Functional Integration

Link metabolites back to organism capabilities.

**Approach:**
1. Identify metabolites of interest
2. Propose biochemical pathways that produce/consume them
3. Check which MAGs encode genes for these pathways
4. Verify: Do MAGs with pathway genes correlate with metabolite abundance?

**Example: Butyrate production**

```r
# Find MAGs annotated with butyrate synthesis genes
butyrate_synthase_kegg <- "K00929"  # Example KEGG ID
mag_has_butyrate_gene <- sapply(mag_annotations, 
                                function(x) butyrate_synthase_kegg %in% x)

# Correlate these MAGs with butyrate metabolite abundance
# Do they correlate? This suggests functional relationship
```

---

## Part 5: Visualization Best Practices

### Network Visualization

**Effective network plots show:**
- Clear node separation by type (MAG vs metabolite)
- Edge colors/thickness reflecting correlation strength
- Node position reflecting network topology (central vs peripheral)
- Legends explaining all visual encoding

### Heatmaps

```r
# Co-abundance heatmap (MAGs × Metabolites)
correlation_matrix <- cor(t(mag_abundance), t(metabolite_abundance))
pheatmap::pheatmap(correlation_matrix,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean")
```

### Ordination Plots

Show how samples cluster based on composition/metabolism

```r
# Combined PCA plot
combined_pca <- prcomp(rbind(scale(t(mag_abundance)), 
                              scale(t(metabolite_abundance))),
                       scale. = FALSE)
```

---

## Expected Results from DIME Analysis

By end of this section, you should have:

1. ✓ Characterized microbial community composition across treatment groups
2. ✓ Recovered and identified microbial genomes (MAGs)
3. ✓ Measured metabolite profiles across samples
4. ✓ Inferred co-occurrence networks between organisms and metabolites
5. ✓ Identified which organisms potentially produce key metabolites
6. ✓ Quantified the statistical significance of community/metabolite shifts

---

## Troubleshooting Common Issues

**Issue:** SpiecEasi takes very long to run
- Solution: Subset data to highly variable features only; reduce lambda.min.ratio

**Issue:** Network too dense to visualize
- Solution: Filter edges by correlation strength; use network layout algorithms

**Issue:** Few features annotated
- Solution: Use KEGG IDs; focus on known compounds relevant to your hypothesis

---

## Next Steps

Proceed to the Conclusion section to synthesize these results and discuss broader implications!
