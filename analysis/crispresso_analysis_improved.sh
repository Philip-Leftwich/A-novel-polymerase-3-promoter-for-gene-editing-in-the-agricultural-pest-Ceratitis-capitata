#!/bin/bash
################################################################################
# CRISPResso2 Analysis Pipeline for 7SK vs U6 Promoter Comparison
#
# Purpose: Analyze CRISPR-Cas9 editing efficiency in Ceratitis capitata using
#          7SK and U6 promoters to drive guide RNA expression targeting the
#          white eye gene
#
# Author: Philip Leftwich
# Date: 2026-04-20
# Version: 1.0.0
#
# Dependencies:
#   - CRISPResso2 v2.3.3+
#   - fastp
#   - cutadapt
#   - conda/anaconda
#
# Usage: bash crispresso_analysis_improved.sh [options]
################################################################################

set -euo pipefail  # Exit on any error, undefined variable, or pipe failure

################################################################################
# Configuration and Global Variables
################################################################################

# Script metadata
SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Environment settings
CONDA_ENV_NAME="crispresso_trim"
PYTHON_VERSION="3.10"
CRISPRESSO_VERSION="v2.3.3"

# Analysis parameters - these can be modified for different experiments
FORWARD_PRIMER="ACGGCATATGACACAAAAACAG"
REVERSE_PRIMER="TTCTGCGATAGCTTTTTCAACA"
GUIDE_SEQUENCE="CGGGTGAAGGGTTTATCGGG"

# Reference sequences
WILDTYPE_AMPLICON="ACGGCATATGACACAAAAACAGAAAGTACAACGTGTCGATCAAGTTATACAGGACCTCTCGCTGGGTAAATGTCAGAATACGTTGATTGGCGTGCCGGGTCGGGTGAAGGGTTTATCGGGTGGCGAGCGTAAGCGACTGGCATTCGCCTCGGAGGCGTTAACCGATCCGCCGCTATTAATTTGCGATGAACCTACCTCGGGTTTGGACTCGTTCATGGCACACAGCGTCGTACAGGTGTTGAAAAAGCTATCGCAGAA"
MUTANT_AMPLICON="ACGGCATATGACACAAAAACAGAAAGTACAACGTGTCGATCAAGTTATACAGGACCTCTCGCTGGGTAAATGTCAGAATACGTTGATTGGCGTGCCGGGTGGCGAGCGTAAGCGACTGGCATTCGCCTCGGAGGCGTTAACCGATCCGCCGCTATTAATTTGCGATGAACCTACCTCGGGTTTGGACTCGTTCATGGCACACAGCGTCGTACAGGTGTTGAAAAAGCTATCGCAGAA"

# Processing settings
N_PROCESSES=${N_PROCESSES:-8}  # Number of CPU cores to use
THREAD_COUNT=${THREAD_COUNT:-8}

# File paths - modify these for your environment
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw/sequencing_data"
PROCESSED_DATA_DIR="${PROJECT_ROOT}/data/processed"
OUTPUT_DIR="${PROJECT_ROOT}/outputs/crispresso_results"
LOG_DIR="${PROJECT_ROOT}/outputs/logs"

################################################################################
# Utility Functions
################################################################################

# Logging function
log() {
    local level="$1"
    shift
    local message="$*"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$level] $message" | tee -a "$LOG_DIR/crispresso_${TIMESTAMP}.log"
}

log_info() { log "INFO" "$@"; }
log_warn() { log "WARN" "$@"; }
log_error() { log "ERROR" "$@"; }

# Error handling
error_exit() {
    log_error "$1"
    exit 1
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Validate file exists and is readable
validate_file() {
    local file="$1"
    local description="$2"

    if [[ ! -f "$file" ]]; then
        error_exit "$description not found: $file"
    fi

    if [[ ! -r "$file" ]]; then
        error_exit "$description not readable: $file"
    fi

    log_info "Validated $description: $file"
}

# Create directory if it doesn't exist
ensure_dir() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        log_info "Created directory: $dir"
    fi
}

################################################################################
# Environment Setup Functions
################################################################################

# Setup conda environment
setup_environment() {
    log_info "Setting up CRISPResso2 environment..."

    # Load conda module (adjust for your HPC system)
    if command_exists module; then
        module load python/anaconda/2024.10/3.12.7 || {
            log_warn "Could not load anaconda module - assuming conda is already available"
        }
    fi

    # Check if conda environment exists
    if ! conda env list | grep -q "$CONDA_ENV_NAME"; then
        log_info "Creating new conda environment: $CONDA_ENV_NAME"
        conda create -n "$CONDA_ENV_NAME" \
            -c conda-forge -c bioconda \
            python="$PYTHON_VERSION" cutadapt fastp crispresso2 -y
    else
        log_info "Conda environment $CONDA_ENV_NAME already exists"
    fi

    # Activate environment
    log_info "Activating conda environment..."
    source activate "$CONDA_ENV_NAME" || conda activate "$CONDA_ENV_NAME"

    # Install/upgrade CRISPResso2 from GitHub
    log_info "Installing/upgrading CRISPResso2..."
    pip install "git+https://github.com/pinellolab/CRISPResso2@${CRISPRESSO_VERSION}"

    # Verify installation
    if ! command_exists CRISPResso; then
        error_exit "CRISPResso2 installation failed"
    fi

    log_info "CRISPResso2 version: $(CRISPResso --version 2>&1 | head -1)"
    log_info "Environment setup complete"
}

# Validate required tools
check_dependencies() {
    log_info "Checking dependencies..."

    local required_tools=("CRISPResso" "fastp" "cutadapt")
    for tool in "${required_tools[@]}"; do
        if ! command_exists "$tool"; then
            error_exit "Required tool not found: $tool"
        fi
        log_info "✓ $tool available"
    done
}

################################################################################
# Data Processing Functions
################################################################################

# Quality control and adapter trimming with fastp
run_fastp() {
    local input_r1="$1"
    local input_r2="$2"
    local output_r1="$3"
    local output_r2="$4"
    local sample_name="$5"

    log_info "Running fastp for $sample_name..."

    local fastp_report="${output_r1%.fastq.gz}_fastp_report"

    fastp \
        -i "$input_r1" \
        -I "$input_r2" \
        -o "$output_r1" \
        -O "$output_r2" \
        --detect_adapter_for_pe \
        --thread "$THREAD_COUNT" \
        --html "$fastp_report.html" \
        --json "$fastp_report.json" \
        --report_title "FastP Report - $sample_name" \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50

    if [[ $? -eq 0 ]]; then
        log_info "FastP completed successfully for $sample_name"
    else
        error_exit "FastP failed for $sample_name"
    fi
}

# Primer removal with cutadapt
remove_primers() {
    local input_r1="$1"
    local input_r2="$2"
    local output_r1="$3"
    local output_r2="$4"
    local sample_name="$5"

    log_info "Removing primers for $sample_name..."

    cutadapt \
        -g "$FORWARD_PRIMER" \
        -G "$REVERSE_PRIMER" \
        -o "$output_r1" \
        -p "$output_r2" \
        --minimum-length 50 \
        --cores "$THREAD_COUNT" \
        "$input_r1" "$input_r2" \
        > "${output_r1%.fastq.gz}_cutadapt.log"

    if [[ $? -eq 0 ]]; then
        log_info "Primer removal completed successfully for $sample_name"
    else
        error_exit "Primer removal failed for $sample_name"
    fi
}

# Run CRISPResso2 analysis
run_crispresso() {
    local fastq_r1="$1"
    local fastq_r2="$2"
    local amplicon_seqs="$3"
    local amplicon_names="$4"
    local output_dir="$5"
    local analysis_name="$6"

    log_info "Running CRISPResso2 analysis: $analysis_name"

    ensure_dir "$output_dir"

    CRISPResso \
        --fastq_r1 "$fastq_r1" \
        --fastq_r2 "$fastq_r2" \
        --amplicon_seq "$amplicon_seqs" \
        --amplicon_name "$amplicon_names" \
        --guide_seq "$GUIDE_SEQUENCE" \
        --output_folder "$output_dir" \
        --name "$analysis_name" \
        --n_processes "$N_PROCESSES" \
        --write_detailed_allele_table \
        --save_also_png \
        --min_frequency_alleles_around_cut_to_plot 0.2 \
        --expand_ambiguous_alignments

    if [[ $? -eq 0 ]]; then
        log_info "CRISPResso2 analysis completed: $analysis_name"
    else
        error_exit "CRISPResso2 analysis failed: $analysis_name"
    fi
}

################################################################################
# Main Analysis Pipeline
################################################################################

# Process a single sample (7SK or U6)
process_sample() {
    local sample_type="$1"  # "7SK" or "U6"
    local raw_r1="$2"
    local raw_r2="$3"

    log_info "Processing $sample_type sample..."

    # Create sample-specific directories
    local sample_processed_dir="$PROCESSED_DATA_DIR/$sample_type"
    ensure_dir "$sample_processed_dir"

    # Define output files for each step
    local fastp_r1="$sample_processed_dir/${sample_type}_fastp_R1.fastq.gz"
    local fastp_r2="$sample_processed_dir/${sample_type}_fastp_R2.fastq.gz"
    local trimmed_r1="$sample_processed_dir/${sample_type}_trimmed_R1.fastq.gz"
    local trimmed_r2="$sample_processed_dir/${sample_type}_trimmed_R2.fastq.gz"

    # Step 1: Quality control and adapter trimming
    run_fastp "$raw_r1" "$raw_r2" "$fastp_r1" "$fastp_r2" "$sample_type"

    # Step 2: Remove PCR primers
    remove_primers "$fastp_r1" "$fastp_r2" "$trimmed_r1" "$trimmed_r2" "$sample_type"

    # Step 3: CRISPResso2 analysis with wildtype reference only
    local wt_output_dir="$OUTPUT_DIR/${sample_type}_wildtype_${TIMESTAMP}"
    run_crispresso \
        "$trimmed_r1" "$trimmed_r2" \
        "$WILDTYPE_AMPLICON" \
        "WILDTYPE" \
        "$wt_output_dir" \
        "${sample_type}_wildtype_analysis"

    # Step 4: CRISPResso2 analysis with both wildtype and mutant references
    local dual_output_dir="$OUTPUT_DIR/${sample_type}_dual_reference_${TIMESTAMP}"
    run_crispresso \
        "$trimmed_r1" "$trimmed_r2" \
        "$WILDTYPE_AMPLICON,$MUTANT_AMPLICON" \
        "WILDTYPE,WHITE_MUTANT" \
        "$dual_output_dir" \
        "${sample_type}_dual_reference_analysis"

    log_info "$sample_type sample processing completed"

    # Return the paths for potential downstream use
    echo "$wt_output_dir:$dual_output_dir"
}

# Main function
main() {
    log_info "Starting CRISPResso2 analysis pipeline"
    log_info "Script: $SCRIPT_NAME"
    log_info "Working directory: $PWD"
    log_info "Project root: $PROJECT_ROOT"

    # Create necessary directories
    ensure_dir "$PROCESSED_DATA_DIR"
    ensure_dir "$OUTPUT_DIR"
    ensure_dir "$LOG_DIR"

    # Setup environment
    setup_environment
    check_dependencies

    # Define input files
    local sevensk_r1="$RAW_DATA_DIR/7SK-gRNA_R1_001.fastq.gz"
    local sevensk_r2="$RAW_DATA_DIR/7SK-gRNA_R2_001.fastq.gz"
    local u6_r1="$RAW_DATA_DIR/U6-gRNA_R1_001.fastq.gz"
    local u6_r2="$RAW_DATA_DIR/U6-gRNA_R2_001.fastq.gz"

    # Validate input files exist
    validate_file "$sevensk_r1" "7SK R1 FASTQ"
    validate_file "$sevensk_r2" "7SK R2 FASTQ"
    validate_file "$u6_r1" "U6 R1 FASTQ"
    validate_file "$u6_r2" "U6 R2 FASTQ"

    # Process 7SK sample
    log_info "Starting 7SK promoter analysis..."
    sevensk_results=$(process_sample "7SK" "$sevensk_r1" "$sevensk_r2")

    # Process U6 sample
    log_info "Starting U6 promoter analysis..."
    u6_results=$(process_sample "U6" "$u6_r1" "$u6_r2")

    log_info "Analysis pipeline completed successfully!"
    log_info "Results summary:"
    log_info "  7SK results: ${sevensk_results}"
    log_info "  U6 results: ${u6_results}"
    log_info "  Log file: $LOG_DIR/crispresso_${TIMESTAMP}.log"
}

################################################################################
# Script Execution
################################################################################

# Display help if requested
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    cat << EOF
CRISPResso2 Analysis Pipeline

Usage: $0 [options]

This script analyzes CRISPR-Cas9 editing efficiency comparing 7SK and U6
promoters driving guide RNA expression in Ceratitis capitata.

Environment Variables:
  N_PROCESSES    Number of CPU cores for CRISPResso2 (default: 8)
  THREAD_COUNT   Number of threads for fastp/cutadapt (default: 8)

Input Files (expected in data/raw/sequencing_data/):
  - 7SK-gRNA_R1_001.fastq.gz
  - 7SK-gRNA_R2_001.fastq.gz
  - U6-gRNA_R1_001.fastq.gz
  - U6-gRNA_R2_001.fastq.gz

Output Directories:
  - data/processed/          Intermediate processed files
  - outputs/crispresso_results/  Final CRISPResso2 results
  - outputs/logs/           Analysis logs

Options:
  -h, --help    Show this help message

EOF
    exit 0
fi

# Run main function if script is executed (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi