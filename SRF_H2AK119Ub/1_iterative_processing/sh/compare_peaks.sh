#!/bin/bash

sample=$1
peaks1="analysis/peaks2/${sample}_broad_peaks_final.broadPeak"
peaks2="analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak"

# Create BED format for both peak sets
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak1_"NR,$7}' $peaks1 > peaks1.bed
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak2_"NR,$7}' $peaks2 > peaks2.bed

# Find unique peaks in each set
bedtools intersect -v -a peaks1.bed -b peaks2.bed > unique_to_first.bed
bedtools intersect -v -a peaks2.bed -b peaks1.bed > unique_to_second.bed

# Analyze unique peaks
echo "Peaks unique to first file:"
awk '{
    len=$3-$2;
    print "Length:", len, "Score:", $5
}' unique_to_first.bed | head -n 5

echo -e "\nPeaks unique to second file:"
awk '{
    len=$3-$2;
    print "Length:", len, "Score:", $5
}' unique_to_second.bed | head -n 5

# Get overlapping regions
bedtools intersect -wa -wb -a peaks1.bed -b peaks2.bed > overlapping.bed

echo -e "\nSample of overlapping peaks with scores:"
awk '{
    print "Peak1 score:", $5, "Peak2 score:", $10
}' overlapping.bed | head -n 5 