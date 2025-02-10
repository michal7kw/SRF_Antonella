# Description of the analysis

Analysis has be done for the H2AK119Ub.
The experiment was conducted on SNB19 glioma cell line (human) after the treatment with GFP (negative control) or with SES_YAF2, 
which is our construct of interest.
So, the samples are GFP 1, GFP 2 and GFP3; YAF 1, YAF 2 and YAF 3.

## Samples
- Control: GFP_1, GFP_2, GFP_3 
- Treatment: YAF_1, YAF_2, YAF_3

## Analysis Pipeline (1_iterative_processing)

1. **Differential Binding Analysis**
   - Performed using both DiffBind and DESeq2 approaches
   - Selected samples with similar peak counts (27k-37k peaks) for more reliable comparison
   - Used samples: GFP_1, GFP_3, YAF_2, YAF_3
   - Identified differentially bound regions between YAF and GFP conditions

2. **Peak Annotation and Enrichment**
   - Annotated peaks relative to genomic features
   - Generated visualizations including:
     - Peak annotation plots
     - TSS distance distributions
     - Genomic distribution of peaks

3. **Advanced Analysis**
   - Peak width distribution analysis
   - Signal intensity correlation between samples
   - Signal profile analysis around TSS regions
   - Peak clustering based on signal patterns
   - Summary statistics:
     - Total Peaks: 41,371
     - Median Peak Width: 1,442 bp
     - Number of Clusters: 3

4. **Gene Overlap Analysis**
   - Compared YAF-bound genes with SOX2 target genes
   - Generated gene lists:
     - All overlapping genes
     - YAF-specific genes
     - Regulatory region overlapping genes
   - Analyzed enrichment in regulatory regions

5. **GO Enrichment Analysis**
   - Performed on multiple gene sets:
     - All overlapping genes
     - YAF-only genes
     - Regulatory overlapping genes
     - Regulatory YAF-only genes


