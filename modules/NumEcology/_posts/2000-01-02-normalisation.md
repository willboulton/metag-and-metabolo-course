---
title: "Numerical Ecology: Building, Cleaning & Normalising the Feature Table"
description: "Build a phyloseq object from LotuS2 output, strip contaminants, and normalise counts by rarefaction or composition."
---

## Overview

Build a phyloseq object from LotuS2 output, strip unwanted taxa and contaminants, then normalise by rarefaction or composition.

*R / phyloseq · adapted from materials by Ezgi Ozkurt · expanded with Kaden Muffett workflow tools*

---

## Load the Packages

```r
library(vegan)
library(phyloseq)
library(decontam)
library(ggplot2)
library(ape)
library(dplyr)
```

## Build a phyloseq Object from LotuS2 Output

LotuS2 writes the four tables you need. Read them in, then merge them into one object. There are different protocols for object creation from different sources, always validate your phyloseq objects!

```r
#abundances
OTU <- read.delim("~/shared-team/datasets/day3/Arad_input/OTU.txt", row.names=1)
#taxonomy
hiera_BLAST <- read.delim("~/shared-team/datasets/day3/Arad_input/hiera_BLAST.txt",row.names = 1)
#sample data
sd <- read.delim("~/shared-team/datasets/day3/Arad_input/in.map",row.names = 1)
#tree
tree <- read.tree("~/shared-team/datasets/day3/Arad_input/OTUphylo.nwk")


tax      <- as.matrix(hiera_BLAST)
otu_mat  <- otu_table(OTU, taxa_are_rows = TRUE)
taxa_mat <- tax_table(tax)
sd       <- sample_data(sd)
physeq   <- phyloseq(otu_mat, taxa_mat, sd, tree)
physeq
```

This object has 36 samples! It is great for some things, but will not let us practice other tools (like decontam). Load this alternative into the background, we will use it later.

```r
classviz<-readRDS("~/shared-team/datasets/day3/class_viz.rds")
classviz
```

> *How many samples, taxa and metadata variables do you have? Look at each component:*

```r
head(otu_table(physeq))
head(tax_table(physeq))
head(sample_data(physeq))
summary(as.factor(sample_data(physeq)$Population))
```

## Subset Samples

`subset_samples()` keeps samples matching a metadata condition.

```r
phy_root <- subset_samples(physeq, Tissue == "Root")
phy_root
```

**Exercise:** Keep every sample that is *not* root.

```r
phy_nonroot <- subset_samples(physeq, Tissue != "Root")
phy_nonroot
```

*Note:* See that the taxa numbers are still the same? No taxa will be dropped until you subset_taxa()

Or drop shallow samples by read count:

```r
physeq1 <- prune_samples(sample_sums(physeq) > 10000, physeq)
```

## Filter Unwanted Taxa

Remove host organelle reads (if they've made it through) — chloroplast and mitochondria.

> **Keep the `is.na()`**
>
> Filtering on a rank also silently drops taxa that are *unassigned* at that rank. Add `| is.na(rank)` so you only remove confirmed hits, not everything unclassified.

```r
physeq <- subset_taxa(physeq, Order  != "Chloroplast"  | is.na(Order))
physeq <- subset_taxa(physeq, Class  != "Chloroplast"  | is.na(Class))
physeq <- subset_taxa(physeq, Family != "Mitochondria" | is.na(Family))
physeq

classviz <- subset_taxa(classviz, Order  != "Chloroplast"  | is.na(Order))
classviz <- subset_taxa(classviz, Class  != "Chloroplast"  | is.na(Class))
classviz <- subset_taxa(classviz, Family != "Mitochondria" | is.na(Family))
classviz
```

## Remove Contaminants with decontam

> **Negative controls**
>
> If your run included blanks, identify and drop likely contaminants before any diversity analysis. Mark controls in a logical column (here `Is.neg`), then score taxa by their prevalence in controls versus real samples.

First, see what DNA made it into the blanks:

```r
classviz@sam_data$Neg<-classviz@sam_data$Body_cav == "BLNK"
plot_bar(subset_samples(classviz, Neg == TRUE), fill = "Phylum")
```

Identify contaminants by the prevalence method:

```r
contamdf <- isContaminant(classviz, method = "prevalence", neg = "Neg")
table(contamdf$contaminant)
```

Plot prevalence in controls against true samples to sanity-check the call:

```r
ps.pa     <- transform_sample_counts(classviz, function(a) 1 * (a > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Neg == TRUE,  ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Neg == FALSE, ps.pa)
df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos),
                    pa.neg = taxa_sums(ps.pa.neg),
                    contaminant = contamdf$contaminant)
ggplot(df.pa, aes(pa.neg, pa.pos, color = contaminant)) + geom_point() + xlim(c(0,5))+
  xlab("Prevalence (negative controls)") + ylab("Prevalence (true samples)")

#physeq <- prune_taxa(!contamdf$contaminant, physeq)
```

We are going to be more conservative than the prevalence based approach. Prune the contaminants, drop the blanks, and remove now-empty taxa:

```r
?isContaminant

contaminantlist<-taxa_names(subset_taxa(subset_samples(classviz, Neg == TRUE), taxa_sums(subset_samples(classviz, Neg == TRUE))>200))

cleanclassvis <- prune_taxa(!taxa_names(classviz) %in% contaminantlist, classviz)
cleanclassvis <- subset_samples(cleanclassvis, Neg == FALSE)
cleanclassvis <- prune_taxa(taxa_sums(cleanclassvis) > 0, cleanclassvis)
```

## Inspect Read Depth

```r
summary(sample_sums(physeq))
sort(phyloseq::sample_sums(physeq))
hist(sample_sums(physeq))

summary(sample_sums(cleanclassvis))
sort(phyloseq::sample_sums(cleanclassvis))
hist(sample_sums(cleanclassvis))
```

*Question:* How do you think these differences impact downstream analysis?

Plot ASV totals across samples to spot the spread:

```r
cASV <- data.frame(counts  = taxa_sums(physeq),
                   samples = names(taxa_sums(physeq)))
ggplot(cASV, aes(reorder(samples, -counts), counts)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("ASV") + ggtitle("ASV count across samples")
```

## Normalise I — Rarefaction

> **Diagnose first**
>
> Before rarefying, draw a rarefaction curve. A curve that flattens means the sample is saturated; pick a depth past the elbow so you keep most of the diversity.

```r
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


vclass <- psotu2veg(physeq)            # phyloseq -> vegan community matrix
rarecurve(vclass, step = 50, cex = 0.5)

vclass2 <- psotu2veg(cleanclassvis)            # phyloseq -> vegan community matrix
rarecurve(vclass2, step = 50, cex = 0.5)
```

`psotu2veg()` and `pssd2veg()` are small helpers that pull the OTU matrix and the sample table out of a phyloseq object in the orientation vegan expects.

The outputs of our rarecurve are important. We have data that has very uneven sampling effort.

Rarefy to the minimum depth, without replacement, with a fixed seed for reproducibility:

```r
phy.rarefied <- rarefy_even_depth(physeq, 
                  sample.size = min(sample_sums(physeq)), replace = FALSE)
ntaxa(physeq); ntaxa(phy.rarefied)
```

> *Why do we lose taxa when we rarefy? What happens to the deepest sample?*

**Exercise:** Rarefy to a depth *larger* than the smallest sample (say 30000). How many samples drop out?

```r
phy.rrf.notmin <- rarefy_even_depth(physeq, rngseed = 1, sample.size = 30000, replace = FALSE)
```

> **Large datasets**
>
> `rarefy_even_depth` is fine for small studies but slow at scale. `rtk` does the same job much faster:

```r
library(rtk)
data       <- matrix(sample(c(rep(0, 1500), rep(1:10, 500), 1:1000), 1.2e6, replace = TRUE), 1000)
samplesize <- min(colSums(data))
system.time(rtk(input = data, depth = samplesize, threads = 12))
```

## Normalise II — Compositional Transform

> **An alternative to rarefying**
>
> Rarefaction discards reads. Converting counts to proportions keeps every read and is what most ordination and bar-plot code expects. **This is NOT** appropriate with uneven sample depths!

```r
library(microbiome)
classvis_p <- microbiome::transform(cleanclassvis, "compositional")
```

Rarefy *or* transform depending on what comes next: richness metrics want even sampling depth (rarefy); ordinations and composition bars are happy with compositional data. Don't do both blindly.

---

## Next Steps

With a cleaned and normalised feature table in hand, move on to alpha diversity.
