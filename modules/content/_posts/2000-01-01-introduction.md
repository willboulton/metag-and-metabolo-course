---
title: "Introduction: Welcome to Metagenomics and Metabolomics"
---

## What is This Course About?

This course introduces the principles and practice of **integrated metagenomics and metabolomics analysis**, using a real-world dataset as a case study. You'll learn how to:

- Process and analyze metagenomic data to characterize microbial communities
- Generate metagenome-assembled genomes (MAGs) that represent draft genomes from environmental samples
- Integrate metabolomic data to understand the metabolic outputs of these communities
- Co-analyze microbial and chemical data using network analysis and statistical approaches
- Interpret the results in biological context

Rather than covering every possible tool and method, this course focuses on a practical workflow that we use in our research group, demonstrated through the **DIME dataset** — a study on how dietary bioactives influence gut microbiome diversity and metabolism. The DIME study is available as a preprint and is shortly to be published in peer-reviewed form.

---

## Why Study Metagenomics and Metabolomics Together?

### The Microbiome Revolution
Over the past decade, metagenomics has transformed our understanding of microbial communities. Instead of studying only cultivable organisms, we can now sequence DNA directly from environmental samples (like the human gut) to capture the full diversity of microorganisms present. This enables us to:

- Identify previously unknown microbial species
- Understand the structure and composition of communities
- Link microbial genes to ecological functions

### The Metabolome Connection
However, DNA sequences don't directly tell us what these microbes are *doing*. This is where metabolomics comes in. Metabolomics measures the small molecules (metabolites) present in a sample — these are the products of microbial metabolism and can reflect:

- The functional output of microbial communities
- How microbes are responding to their environment
- The biochemical signals influencing host health

### Integration is Powerful
By combining both approaches, we can ask more sophisticated questions:

- **Which organisms are producing which metabolites?**
- **How do changes in community composition relate to changes in metabolism?**
- **Which metabolic pathways are active in the community?**

The DIME study provides an excellent example: when dietary bioactives are consumed, they shift gut microbial composition AND the community's metabolic profile. By analyzing both data types together, we can understand the mechanisms behind these changes.

---

## Who is This Course For?

This course is designed for:

- **PhD students** in biology, microbiology, bioinformatics, or related fields who want hands-on experience with multi-omics analysis
- **Master's students** seeking to understand advanced analytical workflows and their biological interpretation
- **Researchers with domain expertise** (e.g., microbiology, ecology, or chemistry) who are new to metagenomics and metabolomics but have some computational foundation
- **Early-career scientists** looking to develop skills in integrative data analysis

### Prerequisites

We assume you have:
- Basic familiarity with the command line / terminal
- Understanding of foundational biology (genes, proteins, metabolism, microbes)
- Some experience with data analysis or statistical software (R, Python, or similar)
- Comfort with reading scientific literature

If you're missing some of these, don't worry — we'll provide references and explanations, but this course is not a beginner's introduction to programming or biology.

---

## The DIME Dataset: Our Case Study

The **DIME study** ("Dietary bioactives increase gut microbiome diversity and alter host and microbial metabolite profiles") is a human dietary intervention study where participants consumed foods enriched with dietary bioactives (polyphenols, fiber, etc.) and researchers measured:

- **Metagenomic data**: Sequencing of fecal samples to profile the gut microbiome
- **Metabolomic data**: Mass spectrometry analysis of fecal metabolites
- **Host measurements**: Clinical and biochemical markers

This dataset is ideal for learning because it:
- Contains both metagenomics and metabolomics data from the same samples
- Comes from a controlled experiment with clear biological signals
- Uses real, published research that you can read and understand
- Addresses a timely question in microbiome research

---

## Course Structure

Here's what you'll learn:

1. **Theory** — Key concepts in metagenomics, metabolomics, and their integration, including compositional data challenges
2. **Setup** — Installing and configuring the tools you'll use (MATAFILER4, R packages including SpiecEasi and DDT)
3. **Workflow** — Running the complete analysis pipeline on the DIME data, from raw reads through to integrated results
4. **Analysis & Interpretation** — Statistical analysis (vegan, DDT), network visualization (SpiecEasi), and biological interpretation
5. **Conclusions** — Synthesizing results, caveats, and applying the workflow to your own data

---

## How to Use This Course

- **Follow the modules in order** — each builds on previous concepts
- **Hands-on practice** — you'll run actual commands and analyze real data
- **Refer back to theory** — when you encounter new concepts in the workflow, return to the theory section
- **Troubleshooting** — check the Troubleshooting section in Setup if you encounter issues

By the end of this course, you'll have:
- A complete analysis workflow you can apply to your own data
- Understanding of the biological insights gained from integrated analysis
- Confidence in interpreting metagenomics and metabolomics results

Let's get started!
