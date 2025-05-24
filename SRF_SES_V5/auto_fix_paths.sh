#!/bin/bash

# Auto-fix common path issues for Cut&Tag pipeline
echo "🔍 Auto-detecting and fixing common path issues..."

# Find potential genome index locations
echo "Searching for hg38 Bowtie2 indexes..."

POTENTIAL_PATHS=(
    "/home/michal/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as/GRCh38_noalt_as"
    "/home/michal/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as"
    "/home/michal/genomes/hg38/bowtie2/GRCh38_noalt_as"
    "/home/*/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as/GRCh38_noalt_as"
    "/home/*/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as"
    "/opt/genomes/hg38/GRCh38_noalt_as"
    "~/genomes/hg38/GRCh38_noalt_as"
    "/usr/local/share/genomes/hg38/GRCh38_noalt_as"
)

FOUND_INDEX=""
for path in "${POTENTIAL_PATHS[@]}"; do
    # Expand wildcards and tildes
    expanded_path=$(eval echo $path)
    if [[ -f "${expanded_path}.1.bt2" ]]; then
        FOUND_INDEX="$expanded_path"
        echo "✅ Found Bowtie2 index at: $FOUND_INDEX"
        break
    fi
done

if [[ -z "$FOUND_INDEX" ]]; then
    echo "❌ No Bowtie2 hg38 index found in common locations"
    echo ""
    echo "Options to fix this:"
    echo "1. Download pre-built index:"
    echo "   wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz | tar -xzC /home/michal/COMMON_FILES/bowtie2_indexes/hg38/"
    echo ""
    echo "2. Build from FASTA:"
    echo "   bowtie2-build hg38.fa /home/michal/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as"
    echo ""
    echo "3. Manually set path in script"
    exit 1
fi

# Check for hg38.chrom.sizes
echo "Checking for chromosome sizes file..."
CHROM_SIZES_PATHS=(
    "/home/michal/COMMON_FILES/hg38.chrom.sizes"
    "/home/*/COMMON_FILES/hg38.chrom.sizes"
    "~/genomes/hg38.chrom.sizes"
    "/opt/genomes/hg38.chrom.sizes"
)

FOUND_CHROM=""
for path in "${CHROM_SIZES_PATHS[@]}"; do
    expanded_path=$(eval echo $path)
    if [[ -f "$expanded_path" ]]; then
        FOUND_CHROM="$expanded_path"
        echo "✅ Found chromosome sizes at: $FOUND_CHROM"
        break
    fi
done

if [[ -z "$FOUND_CHROM" ]]; then
    echo "⚠️  Chromosome sizes file not found, will download automatically"
    FOUND_CHROM="/home/michal/COMMON_FILES/hg38.chrom.sizes"
fi

# Update the scripts with correct paths
echo ""
echo "🔧 Updating scripts with detected paths..."

# Update both scripts
for script in "create_big_wig_parallel_fixed.sh" "create_big_wig_parallel.sh"; do
    if [[ -f "$script" ]]; then
        echo "Updating $script..."
        sed -i "s|GENOME_INDEX=.*|GENOME_INDEX=\"$FOUND_INDEX\"|" "$script"
        sed -i "s|CHROM_SIZES=.*|CHROM_SIZES=\"$FOUND_CHROM\"|" "$script"
        echo "✅ Updated $script"
    fi
done

echo ""
echo "🎉 Path configuration complete!"
echo "   Genome index: $FOUND_INDEX"
echo "   Chromosome sizes: $FOUND_CHROM"
echo ""
echo "Ready to run:"
echo "   ./run_pipeline.sh" 