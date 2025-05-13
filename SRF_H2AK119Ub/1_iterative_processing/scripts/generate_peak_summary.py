import os
import csv

# Define paths relative to the project root (current working directory)
PROJECT_ROOT = os.getcwd()
METRICS_DIR_RELATIVE = "SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling_v2"
SUMMARY_FILE_NAME = "summary.txt"

# Absolute paths
metrics_dir_absolute = os.path.join(PROJECT_ROOT, METRICS_DIR_RELATIVE)
summary_file_path_absolute = os.path.join(metrics_dir_absolute, SUMMARY_FILE_NAME)

# Sample names from 5_peak_calling_v2.sh
samples = ["GFP_1", "GFP_2", "GFP_3", "YAF_1", "YAF_2", "YAF_3"]
peak_counts_data = []

print(f"Project root (CWD): {PROJECT_ROOT}")
print(f"Attempting to read metric files from directory: {metrics_dir_absolute}")
print(f"Summary file will be written to: {summary_file_path_absolute}")

for sample_name in samples:
    metric_file_name = f"{sample_name}_metrics.csv"
    metric_file_path_absolute = os.path.join(metrics_dir_absolute, metric_file_name)

    if not os.path.exists(metric_file_path_absolute):
        print(f"Warning: Metric file not found for {sample_name} at {metric_file_path_absolute}")
        peak_counts_data.append((sample_name, "MetricFileNotFound"))
        continue

    try:
        with open(metric_file_path_absolute, 'r', newline='') as f_metric:
            reader = csv.reader(f_metric)
            found_peak_count = False
            for row in reader:
                if len(row) == 2 and row[0] == "Final_Peak_Count":
                    peak_counts_data.append((sample_name, row[1]))
                    found_peak_count = True
                    break
            if not found_peak_count:
                print(f"Warning: 'Final_Peak_Count' key not found in {metric_file_path_absolute}")
                peak_counts_data.append((sample_name, "PeakCountKeyMissing"))
    except Exception as e:
        print(f"Error reading {metric_file_path_absolute}: {e}")
        peak_counts_data.append((sample_name, f"ReadError"))

# Ensure the output directory exists
try:
    os.makedirs(os.path.dirname(summary_file_path_absolute), exist_ok=True)
except OSError as e:
    print(f"Error creating directory {os.path.dirname(summary_file_path_absolute)}: {e}")
    # If directory creation fails, the script might not be able to write the summary.
    # Depending on requirements, one might want to exit or raise the error.

try:
    with open(summary_file_path_absolute, 'w') as f_summary:
        f_summary.write("Sample\tPeak_Count\n") # Header
        for sample_name, count in peak_counts_data:
            f_summary.write(f"{sample_name}\t{count}\n")
    print(f"Summary successfully written to {summary_file_path_absolute}")
except Exception as e:
    print(f"Error writing summary file {summary_file_path_absolute}: {e}")