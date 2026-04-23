---
layout: index
published: true
---

This hands-on course walks through an integrated **metagenomics and metabolomics analysis** of the [DIME dataset](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1) — a human dietary intervention study examining how bioactive compounds influence the gut microbiome and its metabolic output.

Rather than surveying every available tool, the course follows a single end-to-end research workflow: **MATAFILER4** for MAG recovery, a custom metabolomics pipeline for LC-MS data, **SpiecEasi** for cross-domain network inference, and **R** (vegan, DESeq2, DDT) for statistical analysis. Each module combines conceptual explanation with working code you can run on the DIME data or adapt to your own datasets.

<br>

#### What you will learn

- Recover and classify metagenome-assembled genomes (MAGs) from shotgun sequencing data
- Process and annotate untargeted LC-MS metabolomics data, including polyphenol metabolite features
- Perform diversity analysis, PERMANOVA, and differential abundance testing
- Build cross-domain co-occurrence networks linking microbial taxa to metabolites using SpiecEasi
- Understand why compositional data requires special treatment (CLR transformation, Aitchison distance)
- Reproduce and extend findings from the DIME study on dietary bioactives and gut health

<br>

#### Prerequisites

- Basic command-line familiarity (running commands, navigating directories)
- Foundational biology (genes, metabolism, microbes)
- Some experience with R or a similar statistical environment
- No prior metagenomics or metabolomics experience required

<br>

#### Course modules

| Module | Description |
|--------|-------------|
| **[Introduction](/Introduction/introduction/)** | The DIME dataset, why combine metagenomics and metabolomics, and who this course is for |
| **[Theory](/Theory/theory/)** | MAGs, untargeted metabolomics, compositional data, SpiecEasi network inference |
| **[Setup](/Setup/setup/)** | Installing MATAFILER4, R packages (SpiecEasi, DDT, vegan), and organising your data |
| **[Workflow](/Workflow/workflow/)** | Data processing pipeline: MAG recovery with MATAFILER4 and metabolomics preprocessing |
| **[Analysis](/Analysis/analysis/)** | Integrated R analysis: diversity, differential abundance, SpiecEasi networks, DDT |
| **[Conclusion](/Conclusion/conclusion/)** | Synthesising results, caveats, and applying the workflow to your own data |

<br>

<a class="btn btn-primary" href="/Introduction/introduction/"><i class="fa fa-play"></i> Start the course</a>
