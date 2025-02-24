V5 Peak Processing:
Two V5 ChIP-seq replicates were processed:
Rep1: Started with 2,246 peaks, after chromosome standardization: 2,057 peaks
Rep2: Started with 2,539 peaks, after chromosome standardization: 2,318 peaks
After merging replicates: 4,375 total peaks
After overlap merging: 1,089 final peaks
Width range: 501-3,774 bp

H2AK119Ub Peak Processing:
Six replicates were processed (3 GFP, 3 YAF):
GFP Rep1: 193,906 peaks
GFP Rep2: 243,427 peaks
GFP Rep3: 162,274 peaks
YAF Rep1: 290,169 peaks
YAF Rep2: 259,835 peaks
YAF Rep3: 221,327 peaks

Critical Finding - Blacklist Filtering:
97.5% (1,062) of V5 peaks were filtered out due to significant overlap with blacklist regions
Only 27 V5 peaks remained after blacklist filtering (1,089 - 1,062)

Technical Concerns:
The extremely high proportion of V5 peaks overlapping with blacklisted regions (97.5%) suggests potential quality issues with the V5 ChIP-seq data
The remaining 27 peaks may not be sufficient for meaningful downstream analysis

This could indicate:
Technical issues with the V5 ChIP-seq experiment
Possible contamination or artifacts
Non-specific binding issues

Recommendations:
Review the V5 ChIP-seq experimental protocol and quality metrics
Consider revalidating the V5 antibody specificity
Examine the raw sequencing data quality
Consider repeating the V5 ChIP-seq experiments if the quality issues cannot be resolved
The most concerning aspect is the extremely high proportion of V5 peaks overlapping with blacklisted regions. Blacklisted regions are genomic regions known to produce artifacts in ChIP-seq experiments, and having 97.5% of peaks in these regions is highly unusual and suggests significant technical issues with the V5 ChIP-seq data.

---

For V5 peaks: 97.5% (1,062 out of 1,089) peaks were filtered out
For H2AK119Ub peaks: 0% of peaks were filtered out in both subsequent filtering steps
This is a crucial difference - while the V5 peaks had severe issues with blacklisted regions, the H2AK119Ub peaks showed no overlap with blacklisted regions. This further supports that there might be specific technical issues with the V5 ChIP-seq experiment, since the H2AK119Ub data from the same analysis pipeline shows much better quality.

The blacklist filtering results show a stark contrast between V5 and H2AK119Ub peaks:

V5 peaks: 97.5% overlap with blacklisted regions (1,062 out of 1,089 peaks filtered)
H2AK119Ub peaks: 0% overlap with blacklisted regions (no peaks filtered)
This comparison strongly suggests that the quality issues are specific to the V5 ChIP-seq experiment and not a general problem with the analysis pipeline or blacklist filtering, since the H2AK119Ub data from the same pipeline shows no such issues.