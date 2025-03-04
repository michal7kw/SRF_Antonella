# SRA Data Download for SRF_SES_V5 Project

This directory contains scripts to download SRA data for the SRF_SES_V5 project, which includes data from glioblastoma cells (SNB19) with SES and SOX2 transduction.

## Dataset Information

The dataset includes 6 SRA runs:

| SRR Accession | BioSample    | GSM ID     | Cell Line | Treatment       |
|---------------|--------------|------------|-----------|-----------------|
| SRR18590288   | SAMN27279093 | GSM6008245 | SNB19     | SES-transduced  |
| SRR18590287   | SAMN27279092 | GSM6008246 | SNB19     | SES-transduced  |
| SRR18590286   | SAMN27279091 | GSM6008247 | SNB19     | SES-transduced  |
| SRR18590285   | SAMN27279090 | GSM6008248 | SNB19     | SOX2-transduced |
| SRR18590284   | SAMN27279089 | GSM6008249 | SNB19     | SOX2-transduced |
| SRR18590283   | SAMN27279088 | GSM6008250 | SNB19     | SOX2-transduced |

## Prerequisites

Before using the download script, you need to have the SRA Toolkit installed. You can download it from:
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

Make sure the SRA Toolkit binaries are in your PATH.

## Usage

The script is designed to run on a SLURM-based HPC system and uses job arrays for parallel downloading.

1. Submit the job to SLURM:
   ```
   sbatch download_sra.sh
   ```

2. The script will:
   - Launch 6 parallel jobs (one for each SRA file)
   - Download each SRA file using `prefetch`
   - Create separate log files for each download in the `logs` directory
   - Generate a metadata file (`sample_metadata.tsv`) with information about each sample

3. Monitor the job status:
   ```
   squeue -u $USER
   ```

## Converting SRA to FASTQ

By default, the script downloads SRA files but doesn't convert them to FASTQ format. To convert the downloaded SRA files to FASTQ, you can uncomment the relevant lines in the script or run:

```bash
fasterq-dump SRR18590288 --split-files --outdir ./
```

Replace `SRR18590288` with the appropriate accession number.

## SLURM Job Array Details

The script uses SLURM job arrays to download all files in parallel:

- Each array task (0-5) downloads one SRA file
- Each task has its own log files (download_sra_[JobID]_[ArrayTaskID].out/err)
- All tasks use the same resource allocation (1 CPU, 4GB memory)
- Maximum runtime is set to 120 hours

## Notes

- The download process may take a significant amount of time depending on your internet connection.
- The total size of all files is approximately 45.5 GB.
- Check the logs directory for any errors that might occur during download.
- You can adjust the SLURM parameters (memory, time, etc.) in the script header if needed. 