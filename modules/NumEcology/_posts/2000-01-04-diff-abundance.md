---
title: "Numerical Ecology: Differential Abundance with ALDEx2"
description: "Find taxa that differ between treatments the compositionally-aware way: CLR-transform the counts, test, and rank hits by effect size."
---

## Overview

Find the individual taxa that differ between two treatments the compositionally-aware way: CLR-transform the counts, test, and rank hits by effect size rather than by p-value alone.

*R / phyloseq · adapted from materials by Ezgi Ozkurt · expanded with class workflow tools by KMM*

---

## Load the Packages

``` r
library(phyloseq)
library(ALDEx2)
library(dplyr)
library(ggplot2)
library(ggrepel)
```

## Load the data and set up the comparison

Diversity told us *whether* the communities differ. Differential abundance asks *which taxa* drive it.

``` r
rare_classvis
filter_low_prevalence <- function(ps, min_reads = 10, min_samples = 20) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  keep <- rowSums(otu >= min_reads) >= min_samples
  prune_taxa(keep, ps)
}
rare_classvis_filtered <- filter_low_prevalence(rare_classvis, min_reads = 10, min_samples = 10)
#Down to 231 taxa.
table(sample_data(rare_classvis_filtered)$Body_cav)      # check your two treatment levels
```

> **Proportions lie**
>
> A read count only means something relative to the rest of its sample. If one taxon blooms, everything else *looks* like it dropped even when nothing changed. ALDEx2 sidesteps this by working in centred log-ratio (CLR) space and by modelling the sampling noise with Monte-Carlo draws from the counts. *this choice is not perfect*.

> **Feed it counts, not proportions**
>
> ALDEx2 does its own normalization. Give it the raw integer table — **never** a compositional or rarefied-to-proportions object.

Keep the two treatments, drop everything else, and remove taxa that are now all-zero:

``` r
pf <- prune_taxa(taxa_sums(rare_classvis_filtered) > 0, rare_classvis_filtered)

conds <- as.character(sample_data(pf)$Body_cav)
```

> **Check your orientation**
>
> ALDEx2 wants **features in rows, samples in columns**. A dada2 table is built samples-in-rows, so it needs flipping — but let the object tell you rather than hard-coding `t()`.

``` r
reads <- as(otu_table(pf), "matrix")
if (!taxa_are_rows(pf)) reads <- t(reads)     
reads <- reads[rowSums(reads) > 0, ]
```

## CLR transform and test

Two levels means a Welch t-test plus an effect size. `aldex.kw` adds a Kruskal–Wallis and glm p-value to the same object.

``` r
x        <- aldex.clr(reads, conds, mc.samples = 128, denom = "all", verbose = FALSE)
#MC.samples value is based on sample size. Always look at the Aldex2 documentation!
x.tt     <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)
x.kw     <- aldex.kw(x)
x.effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
x.all    <- data.frame(x.tt, x.kw, x.effect)
```

> **denom = "all" vs "iqlr"**
>
> `"all"` centrs each sample on every feature. If many taxa are genuinely shifting (an asymmetric MA plot below), that reference drifts — switch to `denom = "iqlr"`, which centres only on the middle-variance features.

## Diagnostics

``` r
par(mfrow = c(1, 2))
aldex.plot(x.all, type = "MA", test = "welch", cutoff.pval = 0.05)   # abundance vs difference
aldex.plot(x.all, type = "MW", test = "welch", cutoff.pval = 0.05)   # dispersion vs difference
```

The effect distribution should sit around zero — a big skew is your cue to reconsider `denom`:

``` r
par(mfrow = c(1, 2))
hist(x.all$diff.btw)     # raw between-group difference
hist(x.all$effect)       # standardized effect - should be centered on 0
```

## Call the significant features

``` r
sig <- x.all %>%
  filter(effect > 1 | effect < -1) %>%
  filter(glm.eBH < 0.05) %>%
  filter(we.eBH  < 0.05)

nrow(sig)
```

> **Rank by effect, not p**
>
> ALDEx2's `effect` is the median CLR difference standardized by within-group scatter — far more reproducible across runs than the p-value. A common call is `|effect| > 1`; we recomend the stricter `|effect| > 2` *and* a BH-corrected p below 0.05. Loosening the threshold until a hit appears is p-hacking.

Effect direction follows the treatment levels alphabetically, so a **positive** effect means enriched in `GVC`, negative means enriched in `Bell`. Attach the taxonomy so the hits are readable:

``` r
sig_tax <- cbind(sig, as(tax_table(pf)[rownames(sig), ], "matrix"))
sig_tax[, c("effect", "we.eBH", "Genus")]
```

## Volcano plot

``` r
ggplot(x.all, aes(effect, -log10(we.eBH))) +
  geom_point(colour = "grey70") +
  geom_point(data = sig, colour = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey") +
  geom_text_repel(data = sig, aes(label = rownames(sig)),
                  min.segment.length = 0, max.overlaps = 15) +
  theme_bw() + ggtitle("Bell vs GVC")
```

> To label with names instead of ASV IDs, add `sig$label <- tax_table(pf)[rownames(sig), "Genus"]` and swap `label = label` into the `geom_text_repel` aesthetic.

**Exercise:** Re-run the CLR transform, diagnostics, and significant-feature steps above with `denom = "iqlr"`. Does the significant set change? Which reference frame do you trust here, and why?

> *Which taxa separate the two treatments — and would you have found them from the diversity analyses alone?*

---

## Next Steps

This concludes the numerical ecology detour. Return to the **Analysis** module to continue the DIME integrated analysis.
