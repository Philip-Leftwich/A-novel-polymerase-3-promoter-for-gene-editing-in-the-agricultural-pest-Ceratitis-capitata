# A novel polymerase 3 promoter for gene editing in the agricultural pest Ceratitis capitata

## Authors

Amber S. K. Hall† , Samuel Matthew Shackleton-Chavez† , Tracey Chapman , Philip T. Leftwich*  


[![DOI](https://zenodo.org/badge/DOI/[DOI-TO-BE-ASSIGNED].svg)](https://doi.org/[DOI-TO-BE-ASSIGNED])


## Overview

This repository contains bioinformatics analysis characterizing the 7SK RNA gene promoter in *Ceratitis capitata* (Mediterranean fruit fly) and related dipteran species. The research evaluates the 7SK polymerase III promoter as an alternative to the widely-used U6 promoter for CRISPR-Cas9 guide RNA expression in agricultural pest management applications.

## Repository Structure

```
project-root/
│
├── README.md                  # This file - project description and usage
├── CITATION.cff               # Machine-readable citation metadata  
├── LICENSE                    # Code (MIT) and data (CC BY 4.0) licenses
├── CLAUDE.md                  # Claude Code AI assistant guidance
├── fair-skill.md              # FAIR principles compliance guide
├── .Rprofile                  # renv activation hook
├── renv.lock                  # Exact R package versions
│
├── data/
│   ├── raw/                   # Unmodified input data (read-only)
│      ├── seqs_7sk.csv       # Primary sequence dataset
│      ├── alignments/        # FASTA sequence alignments
│      ├── plasmid_sequences/ # Cloning vectors (GenBank format)
│      └── sequencing_data/   # Raw FASTQ files with checksums

│
├── analysis/                  # Analysis scripts
│   ├── motif_string_match_script.R  # Main motif detection analysis
│   ├── crispresso_analysis_improved.sh  # Enhanced CRISPResso2 pipeline
│   └── crispresso_config.yaml       # Configuration file for CRISPResso2 analysis
│
├── R/                         # R functions and utilities (currently empty)
├── outputs/                   # Generated figures, tables, reports (currently empty)
├── docs/                      # Documentation and metadata
│   └── data-dictionary.md     # Variable definitions and data sources
```

## Data

### Sources and provenance


**Genomic sequences**: Retrieved from NCBI Genome database using the latest available genome assemblies as of April 2026. Reference sequence accessions are documented in the primary dataset (`data/raw/seqs_7sk.csv`).

**Experimental data**: Sequencing data generated at University of East Anglia sequencing facility using Illumina paired-end sequencing (150 bp reads).

**Plasmid sequences**: Designed and synthesized for this study, containing U6 and 7SK promoter constructs for comparative analysis.

### Access conditions

All data in this repository are openly available under Creative Commons Attribution 4.0 International (CC BY 4.0) license. No restrictions apply to access or reuse.

### Format notes


- **Primary data**: CSV format for interoperability
- **Sequence data**: Standard FASTA and GenBank formats
- **Raw sequencing**: FASTQ format with quality scores
- **Checksums**: MD5 hashes provided for data integrity verification

## Minimum Metadata Standard



## Software Environment


- **R version**: 4.4.1
- **Operating system used in development**: Windows 11 Pro 10.0.26200
- **Package management**: `renv` — restore with `renv::restore()`
- **Key dependencies**: tidyverse (2.0.0), stringdist (0.9.12)

### Restoring the environment

```r
# Install renv if not present
install.packages("renv")

# Restore all packages to the versions recorded in renv.lock
renv::restore()
```


## Reproducing the Analysis


The analysis consists of a single R script that performs motif detection and characterization:

```r

# Load required libraries and run analysis
source("motif_string_match_script.R")
```

### Analysis workflow

1. **Data loading**: Reads 7SK promoter sequences from CSV file
2. **Preprocessing**: Removes whitespace and standardizes sequence format
3. **Motif detection**: Searches for PSEA and TATA box motifs using:
   - Exact pattern matching with regular expressions
   - Fuzzy matching using Hamming distance for variant detection
4. **Position analysis**: Calculates motif positions relative to transcription start sites
5. **Output**: Returns annotated dataset with motif locations and match quality


### Key analysis parameters

- **PSEA motif**: `TAATTCCCAAGTGCTTATTTG` (21 bp insect-specific element)
- **TATA box pattern**: `(A|C|T)(A|T)TA(A|T)A` (regex for canonical TATA variants)
- **Search windows**: 
  - TATA: 30 bp window upstream of promoter end
  - PSEA: 70 bp window further upstream

## Experimental Validation Analysis

### CRISPResso2 Pipeline for Editing Efficiency Assessment

In addition to the computational motif analysis, this repository includes scripts for analyzing experimental CRISPR-Cas9 editing efficiency using CRISPResso2. This analysis compares guide RNA expression driven by 7SK versus U6 promoters targeting the *white* eye gene in *Ceratitis capitata*.

#### Running the CRISPResso2 Analysis

**Prerequisites:**
- HPC environment with conda/anaconda
- CRISPResso2 v2.3.3+
- fastp and cutadapt for read preprocessing
- Paired-end Illumina sequencing data (150 bp reads)

**Quick start:**
```bash
# Navigate to analysis directory
cd analysis

# Run the improved pipeline
bash crispresso_analysis_improved.sh

# Or customize settings by editing crispresso_config.yaml first
```

#### Analysis Workflow

1. **Environment setup**: Creates conda environment with required tools
2. **Quality control**: Uses fastp for adapter trimming and quality filtering
3. **Primer removal**: Removes PCR primers using cutadapt
4. **CRISPR analysis**: Runs CRISPResso2 with both wildtype and mutant reference sequences
5. **Results**: Generates editing efficiency reports, allele frequency tables, and visualization plots

#### Key Parameters

- **Target sequence**: White eye gene amplicon (293 bp)
- **Guide RNA**: `CGGGTGAAGGGTTTATCGGG` (20 nt)
- **PCR primers**: 
  - Forward: `ACGGCATATGACACAAAAACAG`
  - Reverse: `TTCTGCGATAGCTTTTTCAACA`
- **Reference alleles**:
  - Wildtype: Benakeion strain sequence
  - Mutant: Consensus 21 bp deletion in white gene

#### Expected Outputs

- `data/processed/`: Quality-trimmed FASTQ files
- `outputs/crispresso_results/`: CRISPResso2 analysis results
  - Editing efficiency summaries
  - Allele frequency tables
  - Mutation position plots
  - Quality control reports
- `outputs/logs/`: Analysis execution logs

#### Computational Requirements

- **Runtime**: ~5 minutes for both samples
- **Memory**: 8-16 GB RAM recommended
- **Storage**: ~2-5 GB for intermediate and output files
- **CPU**: 8+ cores recommended for parallel processing

## Results and Interpretation

The analysis identifies and characterizes key regulatory motifs in polymerase III promoters across dipteran species, enabling comparison of 7SK promoter efficiency with the standard U6 promoter for CRISPR applications. The experimental validation provides quantitative measures of guide RNA expression and editing efficiency for both promoter systems.

## Licence

<!-- REUSABLE: An explicit, machine-readable licence is required for R1.1. -->

**Code**: MIT License (see LICENSE file)  
**Data**: Creative Commons Attribution 4.0 International (CC BY 4.0)


## Contact



Philip Leftwich  
University of East Anglia  
ORCID iD: [https://orcid.org/0000-0001-9500-6592](https://orcid.org/0000-0001-9500-6592)  
Email: p.leftwich@uea.ac.uk

