---
title: "Conclusions: Synthesising Results and Moving Forward"
---

## Congratulations! 

Over the course of this tutorial, you've learned:

1. **The theory** behind metagenomics and metabolomics, their integration, and why compositional data requires special treatment
2. **The practical workflow** for analysing real-world multi-omics data from raw reads to network-level insights
3. **Computational tools** — MATAFILER4 for MAG recovery, the custom metabolomics pipeline, SpiecEasi for network inference, vegan and DESeq2 for statistical analysis, and DDT for automated analysis workflows
4. **Biological interpretation** of complex, multi-layered datasets in the context of dietary bioactives and gut microbial metabolism

You've worked through a complete analysis pipeline from raw sequencing data to mechanistic hypotheses about microbial metabolism, using the DIME dataset as a real-world case study.

---

## Key Insights from the DIME Analysis

### What Did We Learn?

Based on the workflow you just completed, we discovered:

**[To be customized based on actual analysis results]**

*Example structure:*

1. **Community Composition Changes**
   - Dietary bioactives significantly altered gut microbiome composition
   - Specific taxa expanded under treatment (e.g., butyrate-producing Faecalibacterium)
   - Others contracted (e.g., potentially pathogenic species)

2. **Metabolic Shifts**
   - Polyphenol-derived metabolites appeared in treated samples (e.g., phenolic acids)
   - Short-chain fatty acid (SCFA) production increased
   - Host metabolites (e.g., secondary bile acids) changed with microbial composition

3. **Organism-Metabolite Links**
   - Network analysis connected known polyphenol-degrading taxa to polyphenol metabolite abundance
   - SCFA-producing MAGs correlated with SCFA levels
   - Specific MAGs showed functional annotations for genes encoding these pathways

4. **Mechanistic Understanding**
   - The results support a model where dietary bioactives:
     - Selectively promote certain bacterial taxa
     - These taxa encode polyphenol-degradation and SCFA-synthesis pathways
     - This results in both community remodeling and metabolic changes
     - Downstream effects on host health markers (to be validated)

---

## How These Results Compare to the Original Study

The DIME preprint reported similar findings:
- [Key finding 1 from published study]
- [Key finding 2 from published study]
- [Key finding 3 from published study]

Your analysis successfully **reproduced the main findings** using integrated metagenomics and metabolomics analysis. This demonstrates the power of combining these complementary data types.

**Where your analysis might differ:**
- Different MAG dereplication thresholds might yield slightly different genome sets
- Metabolite annotation depends on reference databases (which continue to improve)
- Statistical thresholds for filtering or significance are somewhat subjective choices

These differences highlight important considerations: **biology is real, but data analysis involves choices that affect conclusions.**

---

## Important Lessons Learned

### 1. **Correlation ≠ Causation**

Network analysis revealed correlations between organisms and metabolites, but this doesn't prove one causes the other. Consider:

- **Indirect relationships**: Organism A might produce metabolite X, but metabolite X is only detected because organism B processes it further
- **Common environmental drivers**: Two things might co-vary because they're both responding to the same ecological condition
- **Confounding factors**: The sampling method, sample handling, or other variables might explain correlations

**What to do:** Use network results as hypotheses to test experimentally, not as definitive proof.

### 2. **Data Quality Determines Results**

A key takeaway: your analysis is only as good as your data.

- **Sequencing depth matters**: Rare organisms may not be detected in shallow samples
- **MAG quality affects conclusions**: Fragmented or contaminated MAGs may lead to incorrect functional predictions
- **Metabolite annotation is partial**: Unidentified features could be important players you're missing
- **Normalization choices matter**: Different normalization approaches can affect conclusions

**What to do:** Always QC your data, understand your normalization choices, and document assumptions.

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

**What to do:** Embrace complexity, but validate with independent evidence.

### 4. **Your Analysis is One Possible Interpretation**

The choices you made in this analysis (parameters, thresholds, visualisations) shaped your conclusions. Different choices might yield slightly different results.

**Critical parameters you chose (directly or indirectly):**
- Minimum MAG completeness/contamination thresholds
- Metabolite filtering criteria (abundance, missing data)
- Normalisation methods
- Statistical significance cutoffs
- Network correlation strength filters
- Visualisation parameters

**What to do:** Be transparent about your choices. Better yet, conduct sensitivity analyses: how do conclusions change if you alter key thresholds?

---

## Applying This Workflow to Your Own Data

You can now apply this pipeline to any metagenomics + metabolomics dataset:

1. **Organise your data** in the same format as the DIME example
   - MAG abundance tables + taxonomy
   - Metabolite features + metadata
   - Sample metadata with experimental variables

2. **Run the metagenomics analysis:**
   ```bash
   matafiler4 --config your_config.yaml --output your_output/
   ```

3. **Preprocess metabolomics data** to abundance tables + annotations

4. **Integrate in R** using the same scripts adapted for your data:
   ```r
   # Modify input file paths and sample metadata
   source("scripts/01_prepare_data.R")
   source("scripts/03_network_analysis.R")
   # ... etc
   ```

5. **Generate figures and tables** for publication/presentation

6. **Validate conclusions** experimentally (e.g., isolate predicted organisms, measure metabolite production in culture)

### Adapting the Workflow

Consider these modifications for different questions:

- **Temporal data?** Use time-series methods instead of cross-sectional analysis
- **Multiple environments?** Stratify analysis by environment; compare network structures
- **Host health markers?** Add the host data as additional variables in network analysis
- **Functional data?** Use functional abundance tables (metagenomics) instead of taxonomic abundance
- **Targeted metabolomics?** Leverage higher identification rates; use absolute quantification instead of relative

---

## Next Steps for Your Own Research

### Option 1: Dive Deeper into DIME

- Extract the MAGs and conduct comparative genomics with reference genomes. We didn't discuss any comparative genomics, however, that would certainly be an interesting project you could continue yourself with the data available
- Explore the metadata more fully - there are sleep-quality data and a much larger number of host dietary metadata available, which are worth investigating
- Look at personalisation - can you make personalised predictions of microbiome composition changes and health based on more complex modelling than the GLMs we used?

### Option 2: Apply to Your Own Data

- Identify a dataset combining metagenomics and metabolomics
- Adapt this workflow to your specific question
- Compare results to existing publications on your system

### Option 3: Deepen Your Toolbox

Beyond this course, consider learning:

- **Advanced bioinformatics**: More sophisticated assembly, binning, and functional prediction methods, potentially incorporating metatranscriptomes, and long reads. 
- **Advanced statistics**: Multivariate analysis, causal inference, machine learning 
- **Computational skills**: Python, shell scripting, high-performance computing for larger datasets

---

## Recommended Resources for Continued Learning

### Key Papers

**Metagenomics foundations:**
- [Key metagenomics review to add]
- [Assembly/binning methodology to add]
- [MAG quality considerations to add]

**Metabolomics foundations:**
- [Key metabolomics review to add]
- [Metabolite annotation approaches to add]

**Multi-omics integration:**
- [SpiecEasi methodology paper]
- [Multi-omics integration review to add]
- [The DIME study preprint](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1)

### Online Communities

- [Biostars](https://www.biostars.org) — active Q&A forum for bioinformatics
- [GitHub issues] — for tool-specific questions on repositories

### Workshops and Training

- [Bioinformatics courses]
- [Annual microbiome conferences]
- [Methods courses on metabolomics]

---

## Final Thoughts

### The Big Picture

This course taught you a specific workflow, but the broader skill is understanding how to analyse complex biological data - particularly where there are choices in the analysis that affect interpretation of your results:

1. **Understand your data**: This is true for every single step of a workflow. Getting from sample preparation to abundance tables requires dozens of processing steps; it is a good idea to have some understanding of each. 
2. **Interpret carefully**: Correlation suggests mechanisms, but cannot prove an association or rule out the existance of a confounding variable. A lot of metagenomic analyses are quite exploratory in nature. A bold claim might require other forms of orthogonal evidence for your findings besides correlations. 
3. **Communicate clearly**: Explain choices, limitations, and confidence in conclusions. Often, simple visualisations and statistics (box-plots, t-tests, linear models) are preferable to complex ones because they are more easily interpretable. A complicated analysis may only be necessary if the simple statistics are unclear or underpowered. 

---

## Feedback and Improvements

This course is not set in stone. If you:
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
| DDT | Automated statistical analysis | phyloseq object | Diversity plots, differential abundance tables |
| vegan | Community ecology statistics | Distance matrices | PERMANOVA, ordination |
| DESeq2 / limma | Differential abundance | Count/intensity matrix | Log-fold changes, adjusted p-values |
| phyloseq | Microbiome data structures | Feature + taxonomy + metadata tables | Analysis-ready objects |
| igraph | Network analysis and visualisation | Adjacency matrices | Network plots, centrality metrics |
| ggplot2 | Publication-quality graphics | Data frames | Figures |

### Key Parameters to Report

When publishing your analysis, report:
- MAG completeness and contamination cutoffs - the MIMAG thresholds are always easy to defend
- Extraction method - the Quigen extraction kits are pretty good and widely used. 
- Sequencing depth and preprocessing - (reads per sample),  - 
- Assembly parameters (k-mer sizes, assembler) - MEGAHIT with k=,,,,, is a default and reasonable method. You should report a few assembly quality parameters such as the N50, the percentage of reads mapped to the assembly, 
- Metabolite filtering criteria (abundance threshold, missingness tolerance)
- Normalisation methods - for targeted metabolomics, this includes the use of standards. For sequencing data, converting MAG abundance tables to CPM (counts per million) is fine. 
- Statistical significance cutoffs (α value, p-value correction method) - a p-value of 0.05 is standard. When comparing multiple variables in a pairwise mannar, Benjimani-Hochberg q-value (again, of 0.05) is standard. 
- Network inference method and parameters (SpiecEasi: method, lambda range)
- Uploaded read sequence accession numbers - you'll get these from ENA, or NCBI SRA, when you upload your raw sequencing reads. These should generally not include human-identifiable reads - make sure you filter those out using something like hostile. 

### Checklist for Your Own Analysis

- [ ] Data quality assessment completed
- [ ] Preprocessing documented and code saved
- [ ] Analysis parameters recorded
- [ ] Sensitivity analysis conducted
- [ ] Results validated using independent approach
- [ ] Figures and tables generated
- [ ] Methods section written with sufficient detail for reproduction
- [ ] Limitations discussed
- [ ] Biological conclusions grounded in mechanism
