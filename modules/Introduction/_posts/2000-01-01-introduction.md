---
title: "Introduction: Welcome to Metagenomics and Metabolomics"
---

## What is This Course About?

This course introduces the principles and practice of **integrated metagenomics and metabolomics analysis**, using a real-world dataset as a case study. You'll learn how to:

- Process and analyse metagenomic data to characterise microbial communities
- Generate metagenome-assembled genomes (MAGs) that represent draft genomes from environmental samples
- Integrate metabolomic data to understand the metabolic outputs of these communities
- Interpret the results in biological context

Rather than covering every possible tool and method, this course focuses on a practical workflow that we use in our research group, demonstrated through the **DIME dataset** — a study on how dietary bioactives influence gut microbiome diversity and metabolism. The DIME study is available as a [preprint](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1.full.pdf) and is shortly to be published in peer-reviewed form.

We're also going to focus the practical on the network analysis rather than the two metagenomics and metabolomics pipelines. This is because running either of these pipelines is very resource-heavy computationally, and even just installing all the software could become a multi-day exercise depending on the speed of your internet connection. It would also need to run on an institute HPC; this might require talking to your system administrator about how node scratch space, module software, and SLURM queues are set up. We therefore focus on some of the downstream analyses that make use of 'data products' like MAG abundance and metabolite concentration tables. 

---

## Why Study Metagenomics and Metabolomics Together?

### The Microbiome 
Over the past decade, metagenomics has transformed our understanding of microbial communities. Instead of studying only cultivable organisms, we can now sequence DNA directly from environmental samples (like the human gut) to capture the full diversity of microorganisms present. This enables us to:

- Identify previously unknown microbial species
- Understand the structure and composition of communities
- Link microbial genes to ecological functions

### The Metabolome 
However, DNA sequences don't directly tell us what these microbes are *doing*. This is where metabolomics comes in. Metabolomics measures the small molecules (metabolites) present in a sample - these are the products of microbial metabolism and can reflect:

- The functional output of microbial communities
- How microbes are responding to their environment
- The biochemical signals influencing host health

### Integration of Complementary Data
By combining both approaches, we can ask more sophisticated questions:

- **Which organisms are producing which metabolites?**
- **How do changes in community composition relate to changes in metabolism?**
- **Which metabolic pathways are active in the community?**

The DIME study provides an excellent example: when dietary bioactives are consumed, they shift gut microbial composition AND the community's metabolic profile. By analysing both data types together, we can understand the mechanisms behind these changes.

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

We will provide references and explanations, but this course is not a beginner's introduction to programming or biology.

---

## The DIME Dataset: Our Case Study

The [**DIME study**](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1.full.pdf) ("Dietary bioactives increase gut microbiome diversity and alter host and microbial metabolite profiles") is a human dietary intervention study where participants consumed foods enriched with dietary bioactives (polyphenols, carotenoids, etc.), while keeping dietary fibre constant. Researchers measured:

- **Metagenomic data**: Whole genome shotgun sequencing of fecal samples to profile the gut microbiome
- **Metabolomic data**: Mass spectrometry analysis of urinary and faecal metabolites
- **Host measurements**: Clinical and biochemical markers, e.g. blood glucose monitoring, blood pressure measurements, pulse rate measurements

This dataset is ideal for learning because it:
- Contains both metagenomics and metabolomics data from the same samples
- Comes from a controlled experiment with clear biological signals
- Addresses a timely question in microbiome research

It is worth [reading the preprint](https://www.medrxiv.org/content/10.1101/2025.10.16.25338140v1.full.pdf) to familiarise yourself with the experimental setup. 

---

## Course Structure

Here's what you'll learn:

1. **Theory** — Key concepts in metagenomics, metabolomics, and their integration, including compositional data challenges
2. **Setup** — Installing and configuring the tools you'll use (MATAFILER4, R packages including SpiecEasi and DDT)
3. **Workflow** — Running the complete analysis pipeline on the DIME data, from raw reads through to integrated results
4. **Analysis & Interpretation** — Statistical analysis (vegan, DDT), network visualisation (SpiecEasi), and biological interpretation
5. **Conclusions** — Synthesising results, caveats, and applying the workflow to your own data

---

## How to Use This Course

- **Follow the modules in order** - each builds on previous concepts
- **Hands-on practice** - try to follow along with the 'Analysis' section, you should be able to run the code and reproduce results from the DIME study. 
- **Refer back to theory** - when you encounter new concepts in the workflow, return to the theory section

By the end of this course, you'll have:
- An overview of a metagenomics and metabolomics workflow, that you could adapt to your own data
- Understanding of some of the tools used to perform an integrated analysis
- Understanding of the biological insights gained from integrated analysis
- Confidence in interpreting metagenomics and metabolomics results

Good luck, and DON'T PANIC!