#!/bin/bash

# Ensure all necessary arguments are provided
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <data_file> <conditions_file> <output_directory> <project_name> <correlation_method>"
  echo "Arguments:"
  echo "  data_file           - Path to the data file (e.g., gene expression matrix)"
  echo "  conditions_file     - Path to the conditions file (shared across scripts)"
  echo "  output_directory    - Directory to save all output files"
  echo "  project_name        - Name of the project (used in file naming and titles)"
  echo "  correlation_method  - Correlation method for the heatmap (e.g., 'pearson' or 'spearman')"
  exit 1
fi

# Input arguments
DATA_FILE=$1
CONDITIONS_FILE=$2
OUTPUT_DIR=$3
PROJECT_NAME=$4
CORR_METHOD=$5

# Check if the input files exist
if [ ! -f "$DATA_FILE" ]; then
  echo "Error: Data file '$DATA_FILE' does not exist."
  exit 1
fi

if [ ! -f "$CONDITIONS_FILE" ]; then
  echo "Error: Conditions file '$CONDITIONS_FILE' does not exist."
  exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Log input parameters
echo "Data File: $DATA_FILE"
echo "Conditions File: $CONDITIONS_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Project Name: $PROJECT_NAME"
echo "Correlation Method: $CORR_METHOD"

module load R/4.2.0
# 1. Run PCA Script
echo "Running PCA analysis..."
Rscript /home/ttm3567/63_tylert/Analysis_Algorithms/PCA_Histo_CorrMatrix_pipeline/PCA_script_v4.R "$DATA_FILE" "$CONDITIONS_FILE" "$PROJECT_NAME" "$OUTPUT_DIR"

# 2. Run Histogram Script
echo "Running Histogram analysis..."
Rscript /home/ttm3567/63_tylert/Analysis_Algorithms/PCA_Histo_CorrMatrix_pipeline/Histogram_v3.R "$DATA_FILE" "$CONDITIONS_FILE" "$PROJECT_NAME" "$OUTPUT_DIR"

module unload R/4.2.0


# 3. Run Correlation Heatmap Script
echo "Running Correlation Matrix analysis..."
~/.conda/envs/corr_heatmap_env/bin/python /home/ttm3567/63_tylert/Analysis_Algorithms/PCA_Histo_CorrMatrix_pipeline/Corr_Heatmap_CmdLine_40over.py \
  -name "$PROJECT_NAME" \
  -df "$DATA_FILE" \
  -sl "$CONDITIONS_FILE" \
  -op "$OUTPUT_DIR" \
  -corr "$CORR_METHOD"

echo "All analyses completed. Results saved in '$OUTPUT_DIR'."
