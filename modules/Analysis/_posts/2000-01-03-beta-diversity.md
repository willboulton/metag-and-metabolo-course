---
title: "Analysis: Beta Diversity Ordination and Distance-Based Redundancy Analysis"
---

## Overview

This section reproduces panels **A** (PCoA) and **B** (dbRDA) from Figure 2 of the DIME
paper. These plots show how the gut microbiome community structure changes in response to
the High Bioactive and Low Bioactive intervention diets, using ordination methods.

The analysis uses **species-level** metagenomic abundance data.

---

## 1  Load packages and data-loading helpers

We load the standard packages directly. The only project-specific code we source is
`data.R` (for `read_taxa()` and `read_metadata()`, which are quite simple functions you 
should take a look at; they read and format the data better) and `theme.R` 
(for the shared colour palette and ggplot2 theme).

```r
library(tidyverse)
library(vegan)    # vegdist, capscale, envfit, scores, ordiArrowMul
library(ggrepel)  # geom_text_repel for non-overlapping labels

source("../R/data.R")    # read_metadata(), read_taxa(), order_table()
source("../R/theme.R")   # BIOACTIVE_COLORS, BIOACTIVE_COLORS3, THEME_DIME
```

---

## 2  Load data

The sample metadata; we've seen this table before. 

```r
tbl_sample_metadata <- read_metadata("../data/source/sample_metadata.csv")

tbl_sample_metadata |> head()
```

```
## # A tibble: 6 × 8
##   sample_id participant time_point sample_group sample_arm  sequence diet    
##   <chr>     <chr>       <chr>      <chr>        <chr>       <chr>    <chr>   
## 1 01-MP     01          MP         midpoint     midpoint    low_high midpoint
## 2 01-V1     01          V1         before       before_low  low_high low     
## 3 01-V2     01          V2         after        after_low   low_high low     
## 4 01-V3     01          V3         before       before_high low_high high    
## 5 01-V4     01          V4         after        after_high  low_high high    
## 6 02-MP     02          MP         midpoint     midpoint    high_low midpoint
## # ℹ 1 more variable: baseline <lgl>
```

- `sample_id` - matches the column names in the abundance matrix
- `participant` - participant identifier (crossover design, so each person appears multiple
  times)
- `time_point` - V1 = Visit 1 = baseline, V2 = after 1st diet, MP = midpoint, V4 = after 2nd diet
- `sample_arm` - `before_high`, `after_high`, `before_low`, `after_low`
- `baseline` - `TRUE` for the Visit 1 sample; this is reflective of pre-intervention baseline
- `sequence` - Did this participant crossover from low to high, or high to low?


### 2.2 Species-level taxonomic abundance

`read_taxa()` reads a tab-delimited MATAFILER file and does a small amount of
housekeeping (uppercases sample IDs, drops a negative control column). The result is a
tibble with species on rows and samples on columns. We saw a similar table in the previous
session, for genera rather than species. 

```r
# L6 = species-level assignment
tbl_species <- read_taxa("../data/source/microbiome/taxa/MGS.matL6.txt") |>
  order_table()

dim(tbl_species)   # rows = species, cols = 1 name column + samples
```

```
## [1] 815 101
```

### 2.3 Relative abundance (total-sum scaling)

Sequencing depth varies between samples, so raw read counts are not directly comparable.
We convert each column to proportions that sum to 1 using `prop.table(margin = 2)`.
This is called total-sum scaling (TSS), which we also saw last time.

```r
# Pull the species-name column aside, scale, put it back

count_matrix  <- tbl_species |>
  column_to_rownames("L6") |>
  as.matrix()

rel_matrix <- prop.table(count_matrix, margin = 2)   # column proportions

# Confirm: every sample column should now sum to 1
range(colSums(rel_matrix))
```

```
## [1] 1 1
```

---

## 3  PCoA — Panel A

**Principal Coordinates Analysis (PCoA)** is an ordination technique. We start from a
pairwise-distance matrix and project the samples into a low-dimensional space so that
distances on the plot approximate the original distances as closely as possible.

### 3.1 Compute Bray-Curtis distances

**Bray-Curtis dissimilarity** is quite standard for microbiome data. It ranges from 0
(communities identical in composition) to 1 (no shared species). `vegdist()` from `vegan`
computes it. It expects a matrix with **samples on rows**, so we transpose.

```r
# We only want baseline + post-intervention samples (drop before_high / before_low)
tbl_md_after <- tbl_sample_metadata |>
  filter(baseline | sample_arm %in% c("after_high", "after_low")) |>
  mutate(Diet = case_when(
    baseline             ~ "Baseline",
    sample_arm == "after_high" ~ "High Bioactive",
    sample_arm == "after_low"  ~ "Low Bioactive"
  ))

# Subset the abundance matrix to those samples, then transpose
rel_after <- rel_matrix[, tbl_md_after$sample_id]

bc_dist <- vegan::vegdist(t(rel_after), method = "bray")

# Quick look: distance matrix dimensions
as.matrix(bc_dist)[1:4, 1:4] |> round(3)
```

```
##       01-V1 01-V2 01-V4 02-V1
## 01-V1 0.000 0.198 0.270 0.751
## 01-V2 0.198 0.000 0.272 0.769
## 01-V4 0.270 0.272 0.000 0.713
## 02-V1 0.751 0.769 0.713 0.000
```

**Exercise:** Try using some other distance functions in `vegdist` and see what the differences are. Which are suitable for microbiome data? 

### 3.2 Run PCoA with `cmdscale()`

`cmdscale()` performs classical metric multidimensional scaling on the distance matrix.
Requesting `eig = TRUE` also returns eigenvalues, which tell us how much variance each
axis captures.

```r
pcoa_result <- cmdscale(bc_dist, k = 2, eig = TRUE)  

# % variance explained by each axis
eig <- pcoa_result$eig
pct_var <- round(eig[1:2] / sum(eig[eig > 0]) * 100, 1)
cat("PC1:", pct_var[1], "%\nPC2:", pct_var[2], "%\n")
```

```
## PC1: 16 %
## PC2: 8.7 %
```

### 3.3 Build the PCoA plot

`cmdscale()` returns coordinates in `$points`. We join these with the metadata to get
Diet and participant labels, then plot with `ggplot2`.

The lines connect each participant's baseline sample to their two post-intervention
samples, showing the individual trajectories.

```r
# Extract coordinates and join metadata
tbl_pcoa <- pcoa_result$points |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  rename(PC1 = V1, PC2 = V2) |>
  left_join(tbl_md_after, by = "sample_id")

# Axis labels with variance explained
x_lab <- paste0("PC1 (", pct_var[1], "% variance)")
y_lab <- paste0("PC2 (", pct_var[2], "% variance)")

plt_pcoa <- ggplot(tbl_pcoa, aes(x = PC1, y = PC2, colour = Diet)) +
  # Lines: baseline - High Bioactive 
  geom_line(
    data = filter(tbl_pcoa, Diet %in% c("Baseline", "High Bioactive")),
    aes(group = participant),
    colour   = BIOACTIVE_COLORS["High Bioactive"],
    linewidth = 0.15
  ) +
  # Lines: baseline - Low Bioactive
  geom_line(
    data = filter(tbl_pcoa, Diet %in% c("Baseline", "Low Bioactive")),
    aes(group = participant),
    colour   = BIOACTIVE_COLORS["Low Bioactive"],
    linewidth = 0.15
  ) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = BIOACTIVE_COLORS3) +
  labs(
    title = "A  Species PCoA",
    subtitle = "Bray-Curtis dissimilarity · baseline & post-intervention samples",
    x = x_lab, y = y_lab
  ) +
  THEME_DIME

plt_pcoa
```

![PCoA of species Bray-Curtis dissimilarity, coloured by diet arm with participant trajectories]({{ site.baseurl }}/img/analysis/beta-diversity-pcoa.png)

This has plotted the PCoA and linked points from the same participant. 
You can see quite clearly that there is a large per-participant effect. 
It's harder to spot a pattern in the interventions. 

The next section, distance-based redundancy analysis (dbRDA), is a method to
deal with this. In dbRDA, we constrain the ordination
so that it can only represent changes that are attributable to the diet. 

---

## 4  dbRDA — Panel B

**Distance-based Redundancy Analysis (dbRDA)** is a *constrained* ordination implemented
in `vegan` as `capscale()`. Unlike PCoA, it does not show all sources of variation —
instead it finds the ordination axes best explained by predictor variables, here
**diet arm**.

Because every participant ate both diets (crossover design), each person acts as their own
control. The `Condition(participant)` term removes between-person variation before
fitting the diet effect, making this equivalent to a paired analysis.

**Formula:** `abundance ~ Diet + Condition(participant)`

### 4.1 Prepare a filtered metadata table and matched abundance matrix

Once again, we restrict to three timepoints: V1 (baseline), V2 (after first diet), V4 (after second
diet). Together these give one baseline, one after_high, and one after_low sample per
participant.

```r

# Subset the relative abundance matrix and transpose (capscale needs samples × species)
rel_dbrda <- rel_matrix[, tbl_md_after$sample_id]
mat_dbrda  <- t(rel_dbrda)   # samples x species

dim(mat_dbrda)
```

```
## [1]  60 815
```

### 4.2 Run dbRDA with `capscale()`

`capscale()` accepts a matrix (or data frame) of community data on the left side of the
formula. The `distance = "bray"` argument tells it to compute Bray-Curtis distances
internally before running the constrained ordination.

```r
dbrda_result <- vegan::capscale(
  mat_dbrda ~ Diet + Condition(participant),
  data     = tbl_md_after,
  distance = "bray"
)

# How much variance is explained?
s <- summary(dbrda_result)
cat(sprintf(
  "Conditioned (participant effect removed): %.1f%%\n",
  s$partial.chi / s$tot.chi * 100
))
cat(sprintf(
  "Constrained by diet:                      %.1f%%\n",
  s$constr.chi / s$tot.chi * 100
))
cat(sprintf(
  "Residual (unexplained):                   %.1f%%\n",
  s$unconst.chi / s$tot.chi * 100
))
```

```
## Conditioned (participant effect removed): 83.4%
## Constrained by diet:                      1.0%
## Residual (unexplained):                   15.6%
```

### 4.3 Extract scores with `scores()`

`scores()` retrieves the coordinates for sites (samples), centroids (group means), and
biplot vectors from the fitted model.

```r
dbrda_scores <- vegan::scores(
  dbrda_result, display = c("sites", "cn"), scaling = 2
)

# Sites: one row per sample
tbl_sites <- dbrda_scores$sites |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  left_join(tbl_md_after, by = "sample_id")

# Centroids: one row per diet group
tbl_centroids <- dbrda_scores$centroids |>
  as.data.frame() |>
  rownames_to_column("label") |>
  mutate(label = gsub("^Diet", "", label))   # strip "Diet" prefix added by vegan

head(tbl_sites, 4)
```

```
##   sample_id        CAP1       CAP2 participant time_point sample_group
## 1     01-V1 -0.05230122  0.3020460          01         V1       before
## 2     01-V2 -0.90025028  0.4364091          01         V2        after
## 3     01-V4  0.95255150 -0.7384551          01         V4        after
## 4     02-V1 -0.57816531  1.1147867          02         V1       before
##    sample_arm sequence diet baseline           Diet
## 1  before_low low_high  low     TRUE       Baseline
## 2   after_low low_high  low    FALSE  Low Bioactive
## 3  after_high low_high high    FALSE High Bioactive
## 4 before_high high_low high     TRUE       Baseline
```

### 4.4 Build the base dbRDA plot

```r
# Axis labels: % constrained variance captured by each CAP axis
eig_cap  <- dbrda_result$CCA$eig
pct_cap  <- round(eig_cap / sum(eig_cap) * s$constr.chi / s$tot.chi * 100, 1)

plt_dbrda <- ggplot(tbl_sites, aes(x = CAP1, y = CAP2, colour = Diet)) +
  geom_point(size = 1.5) +
  stat_ellipse(linewidth = 0.5) +
  geom_text(
    data   = tbl_centroids,
    aes(x = CAP1, y = CAP2, label = label),
    colour = "black", size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = BIOACTIVE_COLORS3) +
  labs(
    title    = "B  Species dbRDA",
    subtitle = "community ~ Diet + Condition(participant)  ·  Bray-Curtis",
    x = paste0("CAP1 (", pct_cap[1], "% of total variance)"),
    y = paste0("CAP2 (", pct_cap[2], "% of total variance)")
  ) +
  THEME_DIME

plt_dbrda
```

![dbRDA of species composition constrained by diet, with participant effect removed]({{ site.baseurl }}/img/analysis/beta-diversity-dbrda-base.png)

### 4.5 Add species vectors with `envfit()`

`envfit()` fits each species vector onto the ordination by regression and tests
significance by permutation. We overlay arrows for species that are significantly
associated with the constrained axes (p < 0.05), scaled to fill the plot space using
`ordiArrowMul()`. You might find this takes a little time. 

```r
# Fit species onto the dbRDA ordination
species_fit <- vegan::envfit(dbrda_result, mat_dbrda, permutations = 999)

# Scale arrows to fill ~half the ordination space
arrow_scale  <- vegan::ordiArrowMul(species_fit, fill = 2)
arrow_coords <- vegan::scores(species_fit, display = "vectors") * arrow_scale

# Keep only significant species (p ≤ 0.05)
sig_pvals   <- species_fit$vectors$pvals
sig_arrows  <- arrow_coords[sig_pvals <= 0.05, , drop = FALSE]

cat(sum(sig_pvals <= 0.05), "of", length(sig_pvals), "species significantly associated\n")
```

```
## 5 of 815 species significantly associated
```

Exercise: What are the six significantly associated species? 

The species names in this dataset are full lineage strings separated by `;`.  We trim
each to its deepest known level for readability:

```r
trim_taxon_name <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  # Last non-'?' level
  known <- parts[parts != "?"]
  tail(known, 1)
}

tbl_arrows <- sig_arrows |>
  as.data.frame() |>
  rownames_to_column("raw_name") |>
  mutate(taxon = map_chr(raw_name, trim_taxon_name))

head(tbl_arrows, 4)
```

```
##                                                                                                raw_name
## 1           Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;CAG-127;CAG-127 sp900319515
## 2 Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;Butyribacter;Butyribacter sp000436755
## 3      Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides eggerthii
## 4   Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;Butyribacter;Butyribacter intestini
##        CAP1       CAP2                    taxon
## 1  1.377960 -1.4842113      CAG-127 sp900319515
## 2  2.000000 -1.0840485 Butyribacter sp000436755
## 3 -1.998979  0.1680821    Bacteroides eggerthii
## 4  1.673647 -1.4716780   Butyribacter intestini
```

```r
plt_dbrda_annotated <- plt_dbrda +
  geom_segment(
    data = tbl_arrows,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    arrow     = arrow(length = unit(0.2, "cm")),
    colour    = "grey40",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data        = tbl_arrows,
    aes(x = CAP1, y = CAP2, label = taxon),
    size        = 2,
    colour      = "grey20",
    inherit.aes = FALSE
  )

plt_dbrda_annotated
```

![dbRDA plot annotated with significant species vectors]({{ site.baseurl }}/img/analysis/beta-diversity-dbrda-annotated.png)

We see our friends the Butyribacter are two of the six species significantly associated with the constrained axes, and point towards the High Bioactive arm. 

