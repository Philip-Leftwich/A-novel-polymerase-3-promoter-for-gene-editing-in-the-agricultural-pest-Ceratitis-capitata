# Data Dictionary

This document describes all datasets, variables, and file formats used in the 7SK RNA promoter analysis.

## Primary Dataset: `seqs_7sk.csv`

**Source**: Manual curation from genomic databases  
**Format**: CSV (Comma-separated values) and TSV (Tab-separated values)  
**Location**: `data/raw/seqs_7sk.csv`  
**Encoding**: UTF-8  

### Variables

| Variable | Type | Description | Units | Valid Range |
|----------|------|-------------|-------|-------------|
| `id` | Integer | Sequential record identifier | - | 1-9 |
| `ref_seq` | Character | Reference sequence accession number | - | NCBI accession format |
| `species` | Character | Species name (binomial nomenclature) | - | Standard taxonomic names |
| `seq_number` | Integer | Sequence number within species | - | ≥1 |
| `seq_id` | Character | Unique sequence identifier | - | Species code + "7SK" + number |
| `7sk_rna` | Character | Complete 7SK RNA gene sequence | nucleotides | DNA sequence (ATCG) |
| `7sk_promoter` | Character | Promoter region sequence upstream of 7SK gene | nucleotides | DNA sequence (ATCG) |

### Species Codes

| Code | Species | Common Name |
|------|---------|-------------|
| CCAP | *Ceratitis capitata* | Mediterranean fruit fly |
| AAEG | *Aedes aegypti* | Yellow fever mosquito |
| AALB | *Anopheles albimanus* | Malaria mosquito |
| AEAL | *Aedes albopictus* | Asian tiger mosquito |
| AARA | *Anopheles arabiensis* | Malaria mosquito |
| AGAM | *Anopheles gambiae* | African malaria mosquito |
| ASTE | *Anopheles stephensi* | Asian malaria mosquito |
| CQUI | *Culex quinquefasciatus* | Southern house mosquito |
| DM | *Drosophila melanogaster* | Fruit fly |

## Sequence Alignment Data

**Location**: `data/raw/alignments/`

### Files

- `7SK_gene_sequences.fasta` - Original unaligned sequences
- `7SK_gene_alignment_for_tree.fasta` - Multiple sequence alignment
- `7SK_seq_trimmed.fasta` - Quality-trimmed sequences
- `7SK_seq_trimmed_alignment.fasta` - Alignment of trimmed sequences
- `7SK gene consensus tree.fasta` - Consensus sequence for phylogeny

**Format**: FASTA  
**Headers**: Species identifier and sequence type

## Plasmid Sequences

**Location**: `data/raw/plasmid_sequences/`

### Files

- `uea_012-u6-white-eye-grna.gb` - U6 promoter construct (GenBank format)
- `uea_013-7sk-white-eye-grna.gb` - 7SK promoter construct (GenBank format)

**Purpose**: Cloning vectors for gene editing experiments

## Sequencing Data

**Location**: `data/raw/sequencing_data/`  
**Format**: FASTQ (gzip compressed)  
**Sequencing**: Paired-end Illumina sequencing  
**Public Access**: European Nucleotide Archive (ENA) project accession PRJEB111759  

### Files

- `7SK-gRNA_R1_001.fastq.gz` - Forward reads, 7SK construct
- `7SK-gRNA_R2_001.fastq.gz` - Reverse reads, 7SK construct  
- `U6-gRNA_R1_001.fastq.gz` - Forward reads, U6 construct
- `U6-gRNA_R2_001.fastq.gz` - Reverse reads, U6 construct
- `*.md5` - MD5 checksums for data integrity

**Quality**: Raw reads, not quality-filtered

## Metadata Standards Compliance

This dataset follows elements of:
- **MINSEQE** (Minimum Information about a high-throughput SEQuencing Experiment) for sequencing data
- **INSDC** (International Nucleotide Sequence Database Collaboration) standards for sequence identifiers

### Data Provenance

- **Genomic sequences**: Retrieved from NCBI Genome database
- **Reference assemblies**: Latest available genome assemblies as of April 2026
- **Sequencing data**: Generated at University of East Anglia sequencing facility; submitted to ENA under project PRJEB111759
- **Plasmid sequences**: Designed and synthesized for this study

### Missing Data Handling

- Empty sequence regions: Represented as gaps in alignments
- Failed sequencing reactions: Excluded from analysis
- Ambiguous nucleotides: Retained as 'N' characters where present

## File Format Specifications

### CSV/TSV
- Standard comma-separated values (CSV) or tab-separated values (TSV)
- UTF-8 encoding
- Header row with column names
- Suitable for cross-platform data exchange

### FASTA
- Standard nucleotide FASTA format
- Line width: 80 characters maximum
- Headers: `>Species_identifier|sequence_type`

### GenBank
- Standard GenBank flat file format
- Includes feature annotations for promoter regions, guide RNA sequences
- Vector backbone features annotated

### FASTQ
- Illumina 1.8+ quality encoding (Phred+33)
- Read length: 150 bp paired-end
- Quality scores: 0-40 range