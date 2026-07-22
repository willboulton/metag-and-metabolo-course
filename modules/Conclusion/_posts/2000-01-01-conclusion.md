---
title: "Conclusions and Next Steps"
description: "Synthesising the results, key caveats, and how to apply the workflow to your own metagenomics and metabolomics data."
learningResourceType: "reading"
---

## Congratulations! 

Over the course of this tutorial, you've learned:

1. **The theory** behind metagenomics and metabolomics, their integration, and why compositional data requires special treatment
2. **The practical workflow** for analysing real-world multi-omics data from raw reads to network-analysis
3. **Computational tools** — MATAFILER4 for MAG recovery, an outline of a custom metabolomics pipeline, SpiecEasi for network inference, NetCoMi and Vegan for statistical analysis, and numerical ecology
4. **Biological interpretation** of a complex dataset in the context of dietary bioactives and gut microbial metabolism

You've worked through a complete analysis pipeline; we've used several data products (e.g. abundance tables based on MGS, targeted metabolomics concentrations) using the DIME dataset as a real-world case study. However, you have enough theoretical knowledge to be able to run a similar analysis pipeline yourself. 

---

## Key Insights from the DIME Analysis

### What Did We Learn?

Based on the workflow you just completed, we were able to: 

1. Follow how large-scale metagenomics and metabolomics pipelines are run. 

2. Get an overview of the community and metabolite (SCFA) profiles across the participants and interventions. 

3. Summarise the changes in the community through an ordination.

4. Isolate the variance in the microbial community due to a per-participant effect, and separate this from the dietary intervention's effect. 

5. Correlate metabolites and species, and understand how generating a sparse covariance matrix reduces the false positive rate / multiple testing problem which otherwise hampers this kind of analysis. 

We found that: 

1. There were shifts in gut microbiome composition due to the high and low bioactive diets. 

2. However, these changes were not as pronounced as the changes to the metabolite profiles. 

3. We did not see much change in SCFA profiles, however, we were able to link bioactive intakes to certain SCFA producing species. 


---

## Important Lessons 

### 1. **Correlation != Causation**

Network analysis revealed correlations between organisms and metabolites, but this doesn't prove one causes the other. Consider:

- **Indirect relationships**: Organism A might produce metabolite X, but metabolite X is only detected because organism B processes it further
- **Common environmental drivers**: Two things might co-vary because they're both responding to the same ecological condition
- **Confounding factors**: The sampling method, sample handling, data quality (read depth, number of replicates), processing (including batch effects from some bioinformatics pipelines), and compositional factors may all generate confounding influences between factors

Use network results as hypotheses to test experimentally, not as definitive proof. 

### 2. **Data Quality Is An Important Determining Factor**

A key takeaway: your analysis is only as good as your data.

- **Sequencing depth matters**: Rare organisms may not be detected in shallow samples
- **MAG quality affects conclusions**: Fragmented or contaminated MAGs may lead to incorrect functional predictions
- **Metabolite annotation is partial**: Unidentified features could be important players you're missing
- **Normalisation choices matter**: Different normalisation approaches can affect conclusions

**What to do:** Always QC your data, understand your normalisation choices, and document assumptions.

### 3. **Integration Adds Complexity and Power**

Combining metagenomics and metabolomics:

**Advantages:**
- Tests hypotheses about functional relationships
- Reveals mechanisms behind community shifts
- Can identify keystone organisms (highly connected in networks)
- Generates experimental predictions

**Challenges:**
- More data types = more potential for error/artifact
- Correlations don't prove mechanism
- Requires expertise in both domains

### 4. **Choices Made During Analysis**

The choices you made in this analysis (parameters, thresholds, visualisations) shaped your conclusions. Different choices might yield slightly different results.

**Critical parameters (chosen directly or implicitly):**
- Minimum MAG completeness/contamination thresholds
- Metabolite filtering criteria (abundance, missing data)
- Normalisation methods
- Statistical significance cutoffs, multiple testing corrections
- Network correlation strength filters
- Visualisation parameters

Ideally, you could conduct sensitivity analyses: how do conclusions change if you alter key thresholds?


---

## Next Steps for Your Own Research

### Option 1: Dive Deeper into DIME

- Extract the MAGs and conduct comparative genomics with reference genomes. We didn't discuss any comparative genomics, however, that would certainly be an interesting project you could continue yourself with the data available
- Explore the metadata more fully - there are sleep-quality data and a much larger number of host dietary metadata available, which are worth investigating
- Look at personalisation - can you make personalised predictions of microbiome composition changes and health based on more complex modelling than the GLMs we used?

### Option 2: Apply to Your Own Data

- Identify a dataset combining metagenomics and metabolomics
- Adapt this workflow to your specific question; and try to apply the knowledge from the MATAFILER pipeline (or just use the pipeline itself)
- Compare results to existing publications on your system

Beyond this course, consider learning:

- **Advanced bioinformatics**: More sophisticated assembly, binning, and functional prediction methods, potentially incorporating metatranscriptomes, and long reads. 
- **Advanced statistics**: Multivariate analysis, causal inference, machine learning.
- **Computational skills**: R, Python, shell scripting, high-performance computing for larger datasets.

---

## Recommended Resources for Continued Learning

### Key Papers

**Metagenomics foundations:**
- [Quince et al. (2017) Review](https://doi.org/10.1038/nbt.3935)
- [The CAMI2 Assembly and Binning Competitions](https://doi.org/10.1038/s41592-022-01431-4)
- [MIMAG Standards](https://doi.org/10.1038/nbt.3893)

**Metabolomics foundations:**
- [Bauermeister et al. (2021) Review](https://doi.org/10.1038/s41579-021-00621-9)

**Multi-omics integration:**
- [SpiecEasi methodology paper](https://doi.org/10.1371/journal.pcbi.1004226)
- [The DIME study preprint](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1)

### Online Communities

- [Biostars](https://www.biostars.org) — active Q&A forum for bioinformatics
- [GitHub issues] — for tool-specific questions on repositories

### Workshops and Training

- The GIBA training course run by the Quadram Institute. 

---

## Final Thoughts

This course taught you a specific workflow, but the broader skill is understanding how to analyse complex biological data - particularly where there are choices in the analysis that affect interpretation of your results:

1. **Understand your data**: This is true for every single step of a workflow. Getting from sample preparation to abundance tables requires dozens of processing steps; it is a good idea to have some understanding of each. 
2. **Interpret carefully**: Correlation suggests mechanisms, but cannot prove an association or rule out the existance of a confounding variable. A lot of metagenomic analyses are quite exploratory in nature. A bold claim might require other forms of orthogonal evidence for your findings besides correlations. This is especially true when using network-based methods; remember that ultimately the raw data for these are often either a correlation or covariance matrix. 
3. **Communicate clearly**: Explain choices, limitations, and confidence in conclusions. Often, simple visualisations and statistics (box-plots, t-tests, linear models) are preferable to complex ones because they are more easily interpretable. A complicated analysis may only be necessary if the simple statistics are unclear or underpowered. 

---

## Feedback and Improvements

If you:
- Found sections unclear
- Identified errors or outdated information
- Have suggestions for improving topics
- Successfully applied this to your data 

We would be glad to hear about it. Please provide feedback to [this email](@mailto:will.boulton@quadram.ac.uk). 

---

## Thank You

Thank you for working through this course. We hope you now have both the practical skills and conceptual understanding to conduct your own metagenomics and metabolomics analyses.

**CONGRATULATIONS!**

---

## Appendix: Quick Reference

### Key Tools and Their Roles

| Tool | Purpose | Input | Output |
|------|---------|-------|--------|
| MATAFILER4 | End-to-end metagenomic pipeline | Raw FASTQ reads | MAGs, abundance tables, taxonomy, functional annotation |
| Custom metabolomics pipeline | MS data processing and annotation | Raw instrument files | Normalised feature matrix + annotations |
| SpiecEasi | Sparse network inference | MAG/metabolite abundance tables | Conditional independence network |
| vegan | Community ecology statistics | Distance matrices | PERMANOVA, ordination |
| phyloseq | Microbiome data structures | Feature + taxonomy + metadata tables | Analysis-ready objects |
| igraph | Network analysis and visualisation | Adjacency matrices | Network plots, centrality metrics |
| ggplot2 | Publication-quality graphics | Data frames | Figures |

### Key Parameters to Report

When publishing your analysis, report:
- MAG completeness and contamination cutoffs - the MIMAG thresholds are always easy to defend
- Extraction method - the Quigen extraction kits are pretty good and widely used. 
- Sequencing depth and preprocessing - (reads per sample),  - QC method used. 
- Assembly parameters (k-mer sizes, assembler) - MEGAHIT with k=,,,,, is a default and reasonable method. You should report a few assembly quality parameters such as the N50, the percentage of reads mapped to the assembly, for instance. 
- Metabolite filtering criteria (abundance threshold, missingness tolerance)
- Normalisation methods - for targeted metabolomics, this includes the use of standards. For sequencing data, converting MAG abundance tables to CPM (counts per million) is fine. 
- Statistical significance cutoffs (α value, p-value correction method) - a p-value of 0.05 is standard. When comparing multiple variables in a pairwise mannar, Benjimani-Hochberg q-value (again, of 0.05) is standard. 
- Network inference method and parameters (SpiecEasi: method, lambda range)
- Uploaded read sequence accession numbers - you'll get these from ENA, or NCBI SRA, when you upload your raw sequencing reads. These should generally not include human-identifiable reads - make sure you filter those out using something like hostile. 

