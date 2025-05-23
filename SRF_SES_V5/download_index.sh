# Create directory for the index
mkdir -p ~/bowtie2_indexes/hg38
cd ~/bowtie2_indexes/hg38

# Download pre-built index from NCBI
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip

# The index files will be named GRCh38_noalt_as.1.bt2, GRCh38_noalt_as.2.bt2, etc.
# Update your path:
GENOME_INDEX="$HOME/bowtie2_indexes/hg38/GRCh38_noalt_as"