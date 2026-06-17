---
layout: index
published: true
---

This course walks through an integrated **metagenomics and metabolomics analysis** of the [DIME dataset](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1) - a human dietary intervention study examining how bioactive compounds influence the gut microbiome and its metabolic output. 

However, we've intended for this course to have a broader scope than just those interested in the gut microbiome - the same principles should be useful for anyone working with shotgun metagenomics and metabolomics in general. We've also tried to balance some practical learning in a hands-on setting with theory. Practical experience with running code is essential, but just as important is understanding what the commands are doing, so that when it goes wrong in your own data (it will) you can figure out what to do. 

Rather than surveying every available tool, the course follows a single end-to-end research workflow: **MATAFILER4** for MAG recovery, a custom metabolomics pipeline for LC-MS data following **XCMS**, **CAMERA**, and **ProteoWizard**, **SpiecEasi** for cross-domain network inference, and **R** (vegan, DESeq2) and **DDT** for statistical analysis. Each module combines conceptual explanation with working code you can run on the DIME data. The most resource-intensive part of the course is running MATAFILER4. This is a MAG recovery pipeline which is only suitable for running on a high-performance computing cluster. In case this is not available, we provide the output of this pipeline (MAGs), so that you can continue the rest of the analysis from that point onwards. 

<br>

#### What you will learn

- Recover and annotate metagenome-assembled genomes (MAGs) from short read shotgun sequencing data
- Process and annotate targeted and untargeted LC-MS metabolomics data, including polyphenol metabolite features, and short-chain fatty acids (SCFAs)
- Perform alpha and beta diversity analysis, PERMANOVA, and differential abundance testing
- Build cross-domain co-occurrence networks linking microbial taxa to metabolites using SpiecEasi
- Reproduce and extend findings from the DIME study on dietary bioactives and gut health

<br>

#### Prerequisites

- Foundational biology (genes, metabolism, microbes). 
- Basic command-line familiarity (running commands, navigating directories). 
- Some experience with R or a similar statistical environment. 


#### Helpful

- Some familiarity with HPC (high-performance computing), including SLURM or similar. 
- A basic understanding of computing environments (e.g. conda or mamba environments, docker or singularity, etc.) - what they are and why they are a good idea. 
- A little prior metagenomics or metabolomics experience will go a long way, but is not required. If you have done a 16S analysis in the past, that would be enough experience. 

<br>

#### Course modules

| Module | Description |
|--------|-------------|
| **[Introduction](/metag-and-metabolo-course/modules/introduction/introduction/)** | The DIME dataset, why combine metagenomics and metabolomics, and who this course is for |
| **[Theory](/metag-and-metabolo-course/modules/theory/theory/)** | MAGs, metabolomics, compositional data, SpiecEasi network inference |
| **[Setup](/metag-and-metabolo-course/modules/setup/setup/)** | Installing MATAFILER4, R packages (SpiecEasi, DDT, vegan), and organising your data |
| **[Workflow](/metag-and-metabolo-course/modules/workflow/workflow/)** | Data processing pipeline: MAG recovery with MATAFILER4 and metabolomics preprocessing |
| **[Analysis](/metag-and-metabolo-course/modules/analysis/analysis/)** | Integrated R analysis: diversity, differential abundance, SpiecEasi networks, DDT |
| **[Conclusion](/metag-and-metabolo-course/modules/conclusion/conclusion/)** | Synthesising results, caveats, and applying the workflow to your own data |

<br>

<a class="btn btn-primary" href="/metag-and-metabolo-course/modules/introduction/introduction/"><i class="fa fa-play"></i> Start the course</a>
