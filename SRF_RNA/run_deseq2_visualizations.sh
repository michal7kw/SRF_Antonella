#!/bin/bash
#SBATCH --job-name=deseq2_visualizations
#SBATCH --output=logs/deseq2_visualizations.out
#SBATCH --error=logs/deseq2_visualizations.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA"
# Set working directory
cd $BASE_DIR

# Set default values
DDS_DIR="$BASE_DIR/results/deseq2"
COUNTS_DIR="$BASE_DIR/results/deseq2/count_matrix.csv"
METADATA="$BASE_DIR/metadata.csv"
OUTPUT_DIR="$BASE_DIR/results/deseq2_extra"
SPECIES="human"
RUN_PATHWAY=false
TOP_N_DEGS=50
VOLCANO_FC_CUTOFF=1
VOLCANO_PVAL_CUTOFF=0.05

# Function to display usage information
function display_usage {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -d, --dds_dir DIR         Directory containing DESeq2 results (default: $DDS_DIR)"
    echo "  -c, --counts FILE         Counts file (default: $COUNTS_DIR)"
    echo "  -m, --metadata FILE       Metadata file (default: $METADATA)"
    echo "  -o, --output_dir DIR      Output directory for visualization files (default: $OUTPUT_DIR)"
    echo "  -s, --species SPECIES     Species for gene annotation (human or mouse) (default: $SPECIES)"
    echo "  -p, --run_pathway         Run pathway analysis (default: $RUN_PATHWAY)"
    echo "  -n, --top_n_degs N        Number of top DEGs to include in heatmap (default: $TOP_N_DEGS)"
    echo "  -f, --fc_cutoff N         Log2 fold change cutoff for volcano plot (default: $VOLCANO_FC_CUTOFF)"
    echo "  -v, --pval_cutoff N       Adjusted p-value cutoff for volcano plot (default: $VOLCANO_PVAL_CUTOFF)"
    echo "  -h, --help                Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 -d path/to/deseq2 -c path/to/counts -m metadata.csv -o output_dir -s mouse -p"
    echo ""
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dds_dir)
            DDS_DIR="$2"
            shift 2
            ;;
        -c|--counts)
            COUNTS_DIR="$2"
            shift 2
            ;;
        -m|--metadata)
            METADATA="$2"
            shift 2
            ;;
        -o|--output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -s|--species)
            SPECIES="$2"
            shift 2
            ;;
        -p|--run_pathway)
            RUN_PATHWAY=true
            shift
            ;;
        -n|--top_n_degs)
            TOP_N_DEGS="$2"
            shift 2
            ;;
        -f|--fc_cutoff)
            VOLCANO_FC_CUTOFF="$2"
            shift 2
            ;;
        -v|--pval_cutoff)
            VOLCANO_PVAL_CUTOFF="$2"
            shift 2
            ;;
        -h|--help)
            display_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            display_usage
            exit 1
            ;;
    esac
done

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "Error: Rscript is not installed or not in your PATH"
    exit 1
fi

# Check if the R script exists
R_SCRIPT="$BASE_DIR/scripts/standalone_deseq2_visualizations.R"

if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: Cannot find standalone_deseq2_visualizations.R in ${R_SCRIPT}"
    echo "Please make sure the R script is in the same directory as this shell script."
    exit 1
fi

# Check if the directories exist
if [ ! -d "$DDS_DIR" ]; then
    echo "Warning: DESeq2 directory does not exist: ${DDS_DIR}"
    echo "Creating directory ${DDS_DIR}"
    mkdir -p "$DDS_DIR" || { echo "Error: Failed to create directory ${DDS_DIR}"; exit 1; }
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: ${OUTPUT_DIR}"
    mkdir -p "$OUTPUT_DIR" || { echo "Error: Failed to create directory ${OUTPUT_DIR}"; exit 1; }
fi

# Check if metadata file exists
if [ ! -f "$METADATA" ]; then
    echo "Warning: Metadata file does not exist: ${METADATA}"
fi

# Construct the command to run the R script
CMD=(
    "Rscript"
    "${R_SCRIPT}"
    "--dds_dir=${DDS_DIR}"
    "--counts=${COUNTS_DIR}"
    "--metadata=${METADATA}"
    "--output_dir=${OUTPUT_DIR}"
    "--species=${SPECIES}"
    "--top_n_degs=${TOP_N_DEGS}"
    "--volcano_fc_cutoff=${VOLCANO_FC_CUTOFF}"
    "--volcano_pval_cutoff=${VOLCANO_PVAL_CUTOFF}"
)

# Add run_pathway if true
if [ "$RUN_PATHWAY" = true ]; then
    CMD+=("--run_pathway=TRUE")
fi

echo "Running command: ${CMD[@]}"

# Run the R script
"${CMD[@]}"

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "DESeq2 visualization completed successfully!"
    echo "Results are available in: ${OUTPUT_DIR}"
else
    echo "Error: DESeq2 visualization failed with return code $?"
    exit 1
fi 