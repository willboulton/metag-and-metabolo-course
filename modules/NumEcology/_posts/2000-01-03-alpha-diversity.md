---
title: "Numerical Ecology: Alpha Diversity — Within-Sample Richness"
---

## Overview

Visualise composition, estimate alpha diversity, and test group differences by both parametric and non-parametric routes.

*R / phyloseq · adapted from materials by Ezgi Ozkurt · expanded with class workflow tools by KMM*

---

## Load the Packages

``` r
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(cowplot)
```

## Composition bar plots

Agglomerate to a taxonomic rank, then plot the taxa per sample:

``` r
phyPhylum <- tax_glom(phy.rarefied, taxrank = "Phylum")
plot_bar(phyPhylum, fill = "Phylum")
```

> **Comparable bars + facets**
>
> On compositional data the bars become proportions you can compare across samples, and you can split them by metadata.

``` r
classvis_p<-subset_samples(classvis_p,Body_cav=="Bell" |Body_cav =="GVC")
plot_bar( classvis_p, fill = "Phylum") +
  facet_wrap(~ Body_cav, scale = "free") + theme_linedraw()
```

Yikes! Add `geom_bar(stat = "identity")` to collapse the thin per-ASV outlines into solid phylum blocks.

``` r
library(microViz)

classvis_p %>% tax_fix() %>%
  comp_barplot(tax_level = "Order") +
  coord_flip() + facet_wrap(~ Body_cav, scale = "free")

```

Wow wow! Some tidying goes a long way. Under the hood, its vasically this:

```r
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)

# collapse to Order level, convert to relative abundance
ps_order <- classvis_p %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x))

df <- psmelt(ps_order)

# keep top 12 orders, lump the rest as "Other"
top_taxa <- df %>%
  group_by(Order) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 12) %>%
  pull(Order)

df <- df %>%
  mutate(Order = ifelse(Order %in% top_taxa, Order, "Other")) %>%
  group_by(Sample, Body_cav, Order) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# order samples by Bray-Curtis clustering, like microViz's default sample_order
otu_mat <- as(otu_table(ps_order), "matrix")
if (taxa_are_rows(ps_order)) otu_mat <- t(otu_mat)
hc <- hclust(vegdist(otu_mat, method = "bray"), method = "average")
df$Sample <- factor(df$Sample, levels = rownames(otu_mat)[hc$order])

pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(df$Order)))

ggplot(df, aes(Sample, Abundance, fill = Order)) +
  geom_col(colour = "grey20", linewidth = 0.1, width = 0.9) +
  scale_fill_manual(values = pal) +
  coord_flip() +
  facet_wrap(~ Body_cav, scales = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(), legend.title = element_text(face = "bold"))
```

## Alpha diversity metrics

Plot observed richness and Shannon, before and after rarefaction:

``` r
plot_richness(physeq,       x = "Tissue", measures = c("Observed","Shannon")) + geom_boxplot() + theme_bw()
plot_richness(phy.rarefied, x = "Tissue", measures = c("Observed","Shannon")) + geom_boxplot() + theme_bw()
```

> *How do the two compare? What did rarefaction change, and why?*

**Exercise:** Plot Observed, Shannon, Chao1 and InvSimpson across Population (`?plot_richness` for the options).

``` r
p <- plot_richness(phy.rarefied, x = "Population",
       measures = c("Observed","Shannon","Chao1","InvSimpson")) + geom_boxplot() + theme_bw()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

Pull the raw numbers out for testing:

``` r
richness <- estimate_richness(phy.rarefied)
head(richness)
```

Remember, use a diversity scoring system that makes sense. Choosing one with the highest deviation between groups is p-hacking.

### Alpha with `rtk`

> **Rarefy *and* score in one step**
>
> `estimate_richness` scores a **single** rarefied draw, so the numbers wobble with `rngseed`. [`rtk`](https://github.com/hildebra/Rarefaction) (Rarefaction ToolKit) does the rarefying and the diversity math together in C++, repeats it many times, and hands back the **average** — so you lean less on one lucky draw and it stays fast on big tables.

``` r
library(rtk)
```

`rtk` wants a plain count matrix with **features in rows and samples in columns**, so pull the OTU table out of the phyloseq object and transpose if needed:

``` r
otu <- as(otu_table(class_viz), "matrix")
if (!taxa_are_rows(class_viz)) otu <- t(otu)   # rtk expects features (rows) x samples (cols)

depth <- min(colSums(otu))                  # same even depth phyloseq would rarefy to

set.seed(1)
rtk_res <- rtk(otu, repeats = 100, depth = depth, ReturnMatrix = 0, margin = 2)
```

Each sample now carries a *vector* of scores (one per repeat) in `rtk_res$divvs`. Collapse those to one row per sample by averaging across the repeats:

``` r
rtk_alpha <- do.call(rbind, lapply(rtk_res$divvs, function(s) data.frame(
  Sample     = s$samplename,
  Observed   = mean(s$richness),
  Shannon    = mean(s$shannon),
  InvSimpson = mean(s$invsimpson),
  Chao1      = mean(s$chao1)
)))
head(rtk_alpha)
```

Glue the metadata back on and box-plot it just like `plot_richness`:

``` r
library(dplyr)
meta <- data.frame(sample_data(class_viz)); meta$Sample <- sample_names(class_viz)
rtk_alpha <- left_join(rtk_alpha, meta, by = "Sample")

rtk_alpha %>%
  tidyr::pivot_longer(c(Observed, Shannon), names_to = "measure", values_to = "value") %>%
  ggplot(aes(Body_cav, value)) +
  geom_boxplot() + facet_wrap(~ measure, scales = "free_y") + theme_bw()
```

On this dataset the two routes are effectively interchangeable (Pearson *r* ≈ 1.00 for Observed, Shannon, InvSimpson and Chao1; group means match to \<0.1%), and a Kruskal–Wallis of richness across `Body_cav` lands on the same p-value (\~1e-9) either way. The payoff for `rtk` is the averaging over draws and the speed.

## Style and combine figures

> **In your free time — panels with cowplot**
>
> Colour by a second variable, theme each panel, then stack them into one figure and save at publication resolution.

``` r
a1 <- plot_richness(physeq, x = "Tissue", measures = c("Observed","Shannon"), color = "Population") +
  geom_boxplot() + scale_color_manual(values = c("darkgreen","coral","brown4","lightblue")) +
  theme_linedraw() + labs(title = "A.", y = " ")
a2 <- plot_richness(physeq, x = "Population", measures = c("Observed","Shannon"), color = "Tissue") +
  geom_boxplot() + theme_linedraw() + labs(title = "B.", y = " ")

plot_grid(a1, a2, rel_widths = c(2, 2))
ggsave("alpha_panel.png", width = 8, height = 5, units = "in", dpi = 300)  # set a real path, or dont
```

## Testing differences — parametric

The classic route: check the distribution, then two-way ANOVA and Tukey’s HSD.

``` r
datax <- p$data[p$data$variable == "InvSimpson", ]
shapiro.test(datax$value)              # p > 0.05 -> approximately normal

res_aov <- aov(value ~ Tissue + Population, data = datax)
summary(res_aov)
TukeyHSD(res_aov)
```

If Shapiro fails, transform the values (log, sqrt) or switch to the rank-based route below.

## Testing differences — non-parametric

> **Kruskal–Wallis + Dunn**
>
> Diversity values are often skewed and groups small, so a rank-based test is usually safer. Kruskal–Wallis across groups, then Dunn’s test for the pairwise comparisons.

``` r
library(dunn.test)
alpha_sd <- pssd2veg(physeq)
alpha_sd$shannon <- vegan::diversity(psotu2veg(physeq), index = "shannon")

kruskal.test(shannon ~ Tissue, data = alpha_sd)
dunn.test::dunn.test(alpha_sd$shannon, alpha_sd$Tissue, method = "holm")
```

It is pretty clear that we have a real difference in diversity across groups.

You can run Dunn directly on the `estimate_richness` output too:

``` r
rare_classvis<-rarefy_even_depth(subset_samples(cleanclassvis, Body_cav=="Bell" |Body_cav =="GVC") , rngseed=1, sample.size=min(sample_sums(cleanclassvis)), replace=F)

rich <- estimate_richness(rare_classvis)
kruskal.test(rich$Observed, sample_data(rare_classvis)$Body_cav)
```

> **Adjust your p-values**
>
> Holm and BH both correct for multiple comparisons — pick one and report which. Every pairwise test you add inflates the false-positive rate.
