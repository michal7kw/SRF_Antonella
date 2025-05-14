import os
import csv
import argparse

def count_lines_in_file(filepath):
    """Counts the number of lines in a given file."""
    try:
        with open(filepath, 'r') as f:
            return sum(1 for _ in f)
    except Exception as e:
        print(f"Error counting lines in {filepath}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Generate a summary of peak counts from metric files, with fallback to line counts from .broadPeak files.")
    parser.add_argument("data_directory",
                        help="Directory containing the '*_metrics.csv' and/or '*_broad_peaks_final.broadPeak' files, and where 'summary.txt' will be written.")
    
    args = parser.parse_args()

    data_dir_absolute = os.path.abspath(args.data_directory)
    SUMMARY_FILE_NAME = "summary.txt"
    summary_file_path_absolute = os.path.join(data_dir_absolute, SUMMARY_FILE_NAME)

    if not os.path.isdir(data_dir_absolute):
        print(f"Error: Provided data directory '{data_dir_absolute}' does not exist or is not a directory.")
        return

    samples = ["GFP_1", "GFP_2", "GFP_3", "YAF_1", "YAF_2", "YAF_3"]
    peak_counts_data = []

    print(f"Processing data from directory: {data_dir_absolute}")
    print(f"Summary file will be written to: {summary_file_path_absolute}")

    for sample_name in samples:
        metric_file_name = f"{sample_name}_metrics.csv"
        metric_file_path = os.path.join(data_dir_absolute, metric_file_name)
        
        peak_count_value = None
        status_message = ""

        if os.path.exists(metric_file_path):
            print(f"Found metrics file: {metric_file_path}")
            try:
                with open(metric_file_path, 'r', newline='') as f_metric:
                    reader = csv.reader(f_metric)
                    found_key = False
                    for row in reader:
                        if len(row) == 2 and row[0] == "Final_Peak_Count":
                            peak_count_value = row[1]
                            found_key = True
                            break
                    if not found_key:
                        status_message = "PeakCountKeyMissingInMetrics"
                        print(f"Warning: '{status_message}' in {metric_file_path}")
            except Exception as e:
                status_message = "MetricsFileReadError"
                print(f"Error reading {metric_file_path}: {e}")
        else:
            print(f"Warning: Metric file not found for {sample_name} at {metric_file_path}. Attempting to use .broadPeak file.")
            broadpeak_file_name = f"{sample_name}_broad_peaks_final.broadPeak"
            broadpeak_file_path = os.path.join(data_dir_absolute, broadpeak_file_name)

            if os.path.exists(broadpeak_file_path):
                print(f"Found broadPeak file: {broadpeak_file_path}")
                line_count = count_lines_in_file(broadpeak_file_path)
                if line_count is not None:
                    peak_count_value = str(line_count) # Store as string
                    print(f"Using line count ({line_count}) from {broadpeak_file_path} as peak count for {sample_name}.")
                else:
                    status_message = "BroadPeakFileReadError" # Error during line count
            else:
                status_message = "DataUnavailable (NoMetricsOrBroadPeak)"
                print(f"Warning: Neither metrics file nor broadPeak file found for {sample_name} in {data_dir_absolute}.")
        
        if peak_count_value is not None:
            peak_counts_data.append((sample_name, peak_count_value))
        else:
            peak_counts_data.append((sample_name, status_message))

    # Ensure the output directory exists (it should, as it's the input data_dir)
    try:
        os.makedirs(data_dir_absolute, exist_ok=True)
    except OSError as e:
        print(f"Error ensuring directory {data_dir_absolute} exists: {e}")
        # This might be redundant if data_dir_absolute must exist from the start.

    try:
        with open(summary_file_path_absolute, 'w') as f_summary:
            f_summary.write("Sample\tPeak_Count\n") # Header
            for sample_name_val, count_val in peak_counts_data:
                f_summary.write(f"{sample_name_val}\t{count_val}\n")
        print(f"Summary successfully written to {summary_file_path_absolute}")
    except Exception as e:
        print(f"Error writing summary file {summary_file_path_absolute}: {e}")

if __name__ == "__main__":
    main()