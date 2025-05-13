# %%a
import os
import csv

# Define paths relative to the workspace root
BASE_DIR = "D:/Github/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/"
metrics_dir = BASE_DIR, "5_peak_calling_v2"
summary_file_path = os.path.join(metrics_dir, "summary.txt")

# Sample names from 5_peak_calling_v2.sh
samples = ["GFP_1", "GFP_2", "GFP_3", "YAF_1", "YAF_2", "YAF_3"]

peak_counts_data = []

print(f"Looking for metric files in: {os.path.abspath(metrics_dir)}")
print(f"Will write summary to: {os.path.abspath(summary_file_path)}")

for sample_name in samples:
    metric_file_name = f"{sample_name}_metrics.csv"
    metric_file_path = os.path.join(metrics_dir, metric_file_name)

    if not os.path.exists(metric_file_path):
        print(f"Warning: Metric file not found for {sample_name} at {metric_file_path}")
        peak_counts_data.append((sample_name, "MetricFileNotFound"))
        continue

    try:
        with open(metric_file_path, 'r', newline='') as f_metric: # Added newline='' for csv reader
            reader = csv.reader(f_metric)
            found_peak_count = False
            for row in reader:
                if len(row) == 2 and row[0] == "Final_Peak_Count":
                    peak_counts_data.append((sample_name, row[1]))
                    found_peak_count = True
                    break
            if not found_peak_count:
                print(f"Warning: 'Final_Peak_Count' not found in {metric_file_path}")
                peak_counts_data.append((sample_name, "PeakCountKeyMissing"))
    except Exception as e:
        print(f"Error reading {metric_file_path}: {e}")
        peak_counts_data.append((sample_name, f"ReadError"))

# Ensure the output directory exists
try:
    os.makedirs(os.path.dirname(summary_file_path), exist_ok=True)
except OSError as e:
    print(f"Error creating directory {os.path.dirname(summary_file_path)}: {e}")
    # Potentially exit or handle if directory creation is critical and fails

try:
    with open(summary_file_path, 'w') as f_summary:
        f_summary.write("Sample\tPeak_Count\n") # Header
        for sample_name, count in peak_counts_data:
            f_summary.write(f"{sample_name}\t{count}\n")
    print(f"Summary successfully written to {summary_file_path}")
except Exception as e:
    print(f"Error writing summary file {summary_file_path}: {e}")