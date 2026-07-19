# Reanalysing the DIME Dataset: Metagenomics & Metabolomics Analysis

An open online course that walks through an **integrated metagenomics and metabolomics analysis**, using the [DIME dataset](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1) — a human dietary-intervention study of how bioactive compounds influence the gut microbiome and its metabolic output — as a case study.

The material is aimed more broadly than just gut-microbiome researchers: the same principles apply to anyone working with shotgun metagenomics and metabolomics. It deliberately balances hands-on code with the underlying theory, so that when an analysis breaks on your own data you understand *why*. Rather than surveying every tool, it follows a single end-to-end research workflow:

- **MATAFILER4** for MAG (metagenome-assembled genome) recovery
- a custom LC-MS metabolomics pipeline built on **XCMS**, **CAMERA**, and **ProteoWizard**
- **SpiecEasi** for cross-domain network inference
- **R** (vegan, DESeq2, ALDEx2, phyloseq) and **DDT** for statistical analysis

> The most resource-intensive step, MAG recovery with MATAFILER4, needs an HPC cluster. If you don't have one, the pipeline output (MAGs) is provided so you can continue from that point using the [preprocessed data on Figshare](https://doi.org/10.6084/m9.figshare.29860841).

## What you will learn

- Recover and annotate MAGs from short-read shotgun sequencing data
- Process and annotate targeted and untargeted LC-MS metabolomics data, including polyphenol metabolites and short-chain fatty acids (SCFAs)
- Perform alpha/beta diversity analysis, PERMANOVA, and differential abundance testing
- Build cross-domain co-occurrence networks linking microbial taxa to metabolites with SpiecEasi
- Reproduce and extend findings from the DIME study on dietary bioactives and gut health

## Course modules

| Module | Description |
|--------|-------------|
| **Introduction** | The DIME dataset, why combine metagenomics and metabolomics, and who this course is for |
| **Theory** | MAGs, metabolomics, compositional data, and SpiecEasi network inference |
| **Setup** | Installing MATAFILER4, R packages (SpiecEasi, vegan, …) and organising your data |
| **Workflow** | Data processing: MAG recovery with MATAFILER4 and metabolomics preprocessing |
| **Numerical Ecology** | Diversity, normalisation, and compositionally-aware differential abundance |
| **Analysis** | Integrated R analysis: diversity, differential abundance, SpiecEasi networks |
| **Conclusion** | Synthesising results, caveats, and applying the workflow to your own data |

## Prerequisites

- Foundational biology (genes, metabolism, microbes)
- Basic command-line familiarity (running commands, navigating directories)
- Some experience with R or a similar statistical environment

Helpful but not required: familiarity with HPC/SLURM, computing environments (conda/mamba, Docker/Singularity), or prior 16S/metagenomics experience.

## Running locally

This is a [Jekyll](https://jekyllrb.com/) site. The quickest way to preview it with live reload is Docker:

```bash
docker compose up
```

Then open <http://localhost:4000>. Gems are cached in `.bundler/` and the site rebuilds automatically as you edit content. Course content lives under [`modules/`](/modules) — see [`modules/README.md`](modules/README.md) for the folder/section conventions.

## About this template

The site is built on [Course-in-a-Box](https://course-in-a-box.p2pu.org), a free tool by [Peer 2 Peer University](https://www.p2pu.org) for building and publishing online courses. To adapt the template itself, start with the [`_layouts`](/_layouts), [`_includes`](/_includes) and [`css`](/css) directories.

## License

- Template code — MIT License (© 2020 Peer 2 Peer University)
- Course content ("Modules") — [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
