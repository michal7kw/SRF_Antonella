#!/bin/bash
set -e
set -o pipefail

# Get arguments
THREADS=$1
INDEX=$2
R1=$3
R2=$4
OUTPUT=$5
SAMPLE=$6
LOG=$7
TMP_DIR=$8

# Start timing
start_time=$(date +%s)

# Log basic information
echo "[$(date)] Starting alignment for sample $SAMPLE" > "$LOG"
echo "Host: $(hostname)" >> "$LOG"
echo "Working directory: $(pwd)" >> "$LOG"
echo "Available memory: $(free -h)" >> "$LOG"
echo "Available disk space: $(df -h .)" >> "$LOG"
echo "Number of threads: $THREADS" >> "$LOG"

# Log input file sizes
echo "Input file sizes:" >> "$LOG"
ls -lh "$R1" "$R2" >> "$LOG"

# Ensure directories exist
mkdir -p "$(dirname "$OUTPUT")"
mkdir -p "$(dirname "$LOG")"
mkdir -p "$TMP_DIR"

# Check bowtie2 version and availability
echo "Bowtie2 version:" >> "$LOG"
bowtie2 --version | head -n 1 >> "$LOG" 2>&1

# Check if index files exist
echo "Checking index files..." >> "$LOG"
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    if [ ! -f "${INDEX}.$ext" ]; then
        echo "[$(date)] Error: Bowtie2 index file ${INDEX}.$ext not found" | tee -a "$LOG" >&2
        exit 1
    fi
    ls -lh "${INDEX}.$ext" >> "$LOG"
done

# Check input files
echo "Checking input files..." >> "$LOG"
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "[$(date)] Error: Input fastq files not found" | tee -a "$LOG" >&2
    exit 1
fi

# Create a temporary output file
TMP_OUTPUT="$TMP_DIR/temp.sam"

# Log the full command
echo "[$(date)] Running bowtie2 command:" >> "$LOG"
echo "bowtie2 -p $THREADS \\
    -x $INDEX \\
    --local --very-sensitive-local \\
    --no-mixed --no-discordant \\
    -I 0 -X 1000 \\
    --rg-id '$SAMPLE' \\
    --rg 'SM:$SAMPLE' \\
    --rg 'PL:ILLUMINA' \\
    --rg 'LB:lib1' \\
    -1 $R1 \\
    -2 $R2 \\
    -S $TMP_OUTPUT" >> "$LOG"

# Set a timeout value (95% of allocated runtime in seconds)
TIMEOUT=$((1440 * 60 * 95 / 100))
echo "[$(date)] Starting alignment with $TIMEOUT second timeout" >> "$LOG"

# Run bowtie2 with timeout and proper error handling
# Use a separate process group to ensure all child processes are killed by timeout
timeout -k 30 $TIMEOUT bowtie2 -p $THREADS \
    -x $INDEX \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    -I 0 -X 1000 \
    --rg-id "$SAMPLE" \
    --rg "SM:$SAMPLE" \
    --rg "PL:ILLUMINA" \
    --rg "LB:lib1" \
    -1 "$R1" \
    -2 "$R2" \
    -S "$TMP_OUTPUT" 2>> "$LOG"

# Check the exit status
BOWTIE_STATUS=$?
if [ $BOWTIE_STATUS -ne 0 ]; then
    if [ $BOWTIE_STATUS -eq 124 ] || [ $BOWTIE_STATUS -eq 137 ]; then
        echo "[$(date)] Error: bowtie2 timed out after $TIMEOUT seconds" | tee -a "$LOG" >&2
    else
        echo "[$(date)] Error: bowtie2 failed with exit code $BOWTIE_STATUS" | tee -a "$LOG" >&2
    fi
    
    # Clean up temporary file
    rm -f "$TMP_OUTPUT"
    exit 1
fi

# Check if output was created and has size
if [ ! -s "$TMP_OUTPUT" ]; then
    echo "[$(date)] Error: Output file is empty or not created" | tee -a "$LOG" >&2
    exit 1
fi

# Move the temporary file to the final output location
mv "$TMP_OUTPUT" "$OUTPUT"

# Log output file size
echo "[$(date)] Output file size:" >> "$LOG"
ls -lh "$OUTPUT" >> "$LOG"

# Calculate and log runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "[$(date)] Alignment completed in $runtime seconds" >> "$LOG"

# Final memory usage
echo "[$(date)] Final memory state:" >> "$LOG"
free -h >> "$LOG"

echo "[$(date)] Alignment rule completed successfully" >> "$LOG"
