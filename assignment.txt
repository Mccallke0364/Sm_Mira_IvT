Ensure that you have snakemake dowloaded and the enviornment file (iwilladdthis)
Open up a new Terminal (^`)
type snakemake and press enter

# Sm_Mira_IvT

This is the parent repository for the final project for the fall 2024 semester of BIOL343. To embark on this project, fork the repository and edit it as you please. Grading specifications are described below.

## Tasks

You have been provided with 16 FASTQs filles in `/data/classes/2024/fall/biol343/course_files/20240923_LH00283_0144_B22NMVKLT3`. These reads originated from 8 different samples: 4 samples of miracidia that were hatched from eggs derived from the liver and 4 samples of miracidia that were hatched from eggs derived from the intestine. The goal of this experiment is to analyze whether miracidia from the intestine are transcriptomically different from those from the liver and, if so, identify differentially expressed genes. To answer the first question, you will perform a PCA of the normalized counts and analyze the sample clusters. To answer the second question, you will use DESeq2 to identify and plot DEGs.

>[!TIP]
> You may want to read [Winners vs. losers: Schistosoma mansoni intestinal and liver eggs exhibit striking differences in gene expression and immunogenicity](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1012268) before you get started. Our project is a follow up to that study.

All relevent metadata and QC data associated with the RNA extraction and library generation is provided in the `qc/` folder. You should review this information prior to starting the analysis, because it should inform your decision making on read trimming. Notably, the library was generated with the NEBNext Ultra II Directional mRNA kit, sequenced on a full 10B lane of a NovaSeq X with a 2x150bp format.

## Grade specifications

### C

To receive a C in this class, the following are required:

1. The full analysis in an IPython notebook that allows reproduction of the entire analysis. At least two plots (PCA of samples and volcano plot of DEGs) should be included. The analysis should include QC evaluation of:
      1. RNA extraction and library prep,
      2. reads before and after trimming,
      3. alignment,
      4. counting,
      5. PCA clustering,and
      6. differential expression.  
  > [!TIP]  
  > This means that you should include Markdown explaining and evaluating logs and QC metrics, as well as your thought process behind each command and argument.
2. An `environment.yaml` file that includes the specifications for generating a conda environment to use as the kernel for the IPython notebook.  
3. A single MultiQC report that includes QC from each command in the notebook.  
4. All of the associated data. Some of it may be too large to put on GitHub, but the count matrices should be included at the very least.  

  > [!IMPORTANT]  
  >The notebook should be organized, written, and formatted in a professional manner, as if you were going to provide it to a client. A rubric will be provided soon.

### B

To receive a B in this class, you are required to do everything required for a C, as well as the following:

1. Use git to version control the entire analysis. Commit messages should be ***informative***, and there should be at least one commit per remaining week of the semester, which means >4 commits.  

### A

To receive an A in this class, you are required to do everything required for a B and C, as well as the following:

1. A snakefile that allows complete reproduction of the analysis. This should be in the `pipeline/` directory.  
2. At least one YAML specifying conda environments to use with the `--use-conda` argument.  

## Other details

- Since the FASTQs are provided, you don't need an SRARunTable or to fetch any FASTQs.  
- The reads are paired-end, so your code from the previous homework/notebooks will need to be updated accordingly.

## Forking this repository

A fork is a repo that originates from someone else's repo but is now owned by you. It is different from a branch, which is owned by the repo owner. Forking this repository will allow you add this analyis to your personal portfolio of bioinfomatics analyses. To fork a repo, navigate to the repo on GitHub.com and click the Fork button on the top-right.
