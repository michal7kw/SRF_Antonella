#!/usr/bin/env python3
"""
Automated Cut&Tag Data Download Script for GSE200062
Downloads both raw sequencing data and supplementary files
"""

import os
import subprocess
import requests
from urllib.parse import urljoin
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('cutntag_download.log'),
        logging.StreamHandler()
    ]
)

# Sample information
SAMPLES = {
    "GSM6008245": "SMNB19_SES_1",
    "GSM6008246": "SMNB19_SES_2", 
    "GSM6008247": "SMNB19_SES_3",
    "GSM6008248": "SMNB19_SOX2_1",
    "GSM6008249": "SMNB19_SOX2_2",
    "GSM6008250": "SMNB19_SOX2_3"
}

# GSM to SRR mapping
GSM_TO_SRR = {
    "GSM6008245": "SRR18590288",
    "GSM6008246": "SRR18590287",
    "GSM6008247": "SRR18590286",
    "GSM6008248": "SRR18590285",
    "GSM6008249": "SRR18590284",
    "GSM6008250": "SRR18590283"
}

def get_srr_from_gsm(gsm_id):
    """Get SRR accession from GSM ID using predefined mapping"""
    srr_id = GSM_TO_SRR.get(gsm_id)
    if srr_id:
        logging.info(f"Found {srr_id} for {gsm_id}")
        return srr_id
    else:
        logging.error(f"No SRR mapping found for {gsm_id}")
        return None

def download_sra_data(gsm_id, srr_id, output_dir="raw_data"):
    """Download raw sequencing data using SRA toolkit"""
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        logging.info(f"Downloading raw data for {gsm_id} ({srr_id})...")
        
        # Use prefetch to download SRA file
        cmd = ["prefetch", srr_id]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"Prefetch failed for {srr_id}: {result.stderr}")
            return False
        
        # Convert to FASTQ
        cmd = [
            "fastq-dump",
            "--split-files",
            "--gzip",
            "--outdir", output_dir,
            srr_id
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"Fastq-dump failed for {srr_id}: {result.stderr}")
            return False
        
        # Rename files for clarity
        old_r1 = os.path.join(output_dir, f"{srr_id}_1.fastq.gz")
        old_r2 = os.path.join(output_dir, f"{srr_id}_2.fastq.gz")
        new_r1 = os.path.join(output_dir, f"{gsm_id}_{srr_id}_R1.fastq.gz")
        new_r2 = os.path.join(output_dir, f"{gsm_id}_{srr_id}_R2.fastq.gz")
        
        if os.path.exists(old_r1):
            os.rename(old_r1, new_r1)
        if os.path.exists(old_r2):
            os.rename(old_r2, new_r2)
            
        logging.info(f"Successfully downloaded {gsm_id}")
        return True
        
    except Exception as e:
        logging.error(f"Error downloading {gsm_id}: {e}")
        return False

def download_supplementary_files(gsm_id, output_dir="supplementary_files"):
    """Download supplementary files (BED, BigWig) from GEO"""
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        logging.info(f"Downloading supplementary files for {gsm_id}...")
        
        # Construct FTP URL
        gsm_prefix = gsm_id[:7]
        base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_prefix}nnn/{gsm_id}/suppl/"
        
        # Get directory listing
        response = requests.get(base_url)
        if response.status_code != 200:
            logging.error(f"Cannot access supplementary files for {gsm_id}")
            return False
        
        # Parse HTML to find files (simple parsing)
        files = []
        for line in response.text.split('\n'):
            if '.bed.gz' in line or '.bw' in line or '.bigWig' in line:
                # Extract filename from href
                import re
                match = re.search(r'href="([^"]+\.(bed\.gz|bw|bigWig))"', line)
                if match:
                    files.append(match.group(1))
        
        # Download each file
        for filename in files:
            file_url = urljoin(base_url, filename)
            output_path = os.path.join(output_dir, filename)
            
            logging.info(f"Downloading {filename}...")
            response = requests.get(file_url, stream=True)
            
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                logging.info(f"Downloaded {filename}")
            else:
                logging.error(f"Failed to download {filename}")
        
        return True
        
    except Exception as e:
        logging.error(f"Error downloading supplementary files for {gsm_id}: {e}")
        return False

def download_sample(gsm_id, sample_name):
    """Download all data for a single sample"""
    logging.info(f"Processing {gsm_id} ({sample_name})")
    
    # Get SRR accession
    srr_id = get_srr_from_gsm(gsm_id)
    
    if not srr_id:
        logging.error(f"Could not find SRR for {gsm_id}")
        return False
    
    # Download raw data
    raw_success = download_sra_data(gsm_id, srr_id)
    
    # Download supplementary files
    supp_success = download_supplementary_files(gsm_id)
    
    return raw_success and supp_success

def main():
    """Main function to download all samples"""
    logging.info("Starting Cut&Tag data download for GSE200062")
    logging.info(f"Processing {len(SAMPLES)} samples")
    
    # Create output directories
    os.makedirs("raw_data", exist_ok=True)
    os.makedirs("supplementary_files", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    # Option 1: Sequential download (safer, slower)
    # for gsm_id, sample_name in SAMPLES.items():
    #     download_sample(gsm_id, sample_name)
    #     time.sleep(1)  # Be nice to NCBI servers
    
    # Option 2: Parallel download (faster, use with caution)
    max_workers = 6  # Don't overwhelm NCBI servers
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_gsm = {
            executor.submit(download_sample, gsm_id, sample_name): gsm_id 
            for gsm_id, sample_name in SAMPLES.items()
        }
        
        for future in as_completed(future_to_gsm):
            gsm_id = future_to_gsm[future]
            try:
                success = future.result()
                if success:
                    logging.info(f"Successfully completed {gsm_id}")
                else:
                    logging.error(f"Failed to complete {gsm_id}")
            except Exception as e:
                logging.error(f"Exception for {gsm_id}: {e}")
            
            # Small delay between completions
            time.sleep(1)
    
    logging.info("Download process completed!")
    
    # Summary
    logging.info("\n=== Download Summary ===")
    logging.info(f"Raw data files should be in: ./raw_data/")
    logging.info(f"Supplementary files should be in: ./supplementary_files/")
    logging.info(f"Check cutntag_download.log for details")

if __name__ == "__main__":
    main()