# Parallel Processing Configuration Options
# Edit the values below and source this file, then run the main script

echo "=== Parallel Processing Configuration Options ==="
echo ""

# Option 1: Conservative (Recommended for your setup)
# Good balance between parallelization and system stability
export MAX_JOBS=5
export THREADS_PER_JOB=4
echo "Option 1 (Current): $MAX_JOBS jobs × $THREADS_PER_JOB threads = $((MAX_JOBS * THREADS_PER_JOB)) cores"

# Option 2: Maximum Parallelization
# Run all 6 samples simultaneously with fewer threads each
# export MAX_JOBS=6
# export THREADS_PER_JOB=3
# echo "Option 2 (Max Parallel): $MAX_JOBS jobs × $THREADS_PER_JOB threads = $((MAX_JOBS * THREADS_PER_JOB)) cores"

# Option 3: Thread-Heavy
# Fewer parallel jobs but more threads per job (good for I/O heavy steps)
# export MAX_JOBS=3
# export THREADS_PER_JOB=6
# echo "Option 3 (Thread Heavy): $MAX_JOBS jobs × $THREADS_PER_JOB threads = $((MAX_JOBS * THREADS_PER_JOB)) cores"

# Option 4: Conservative
# Leave more cores free for system (good if you're doing other work)
# export MAX_JOBS=4
# export THREADS_PER_JOB=4
# echo "Option 4 (Conservative): $MAX_JOBS jobs × $THREADS_PER_JOB threads = $((MAX_JOBS * THREADS_PER_JOB)) cores"

echo ""
echo "Current configuration will process your 6 samples in approximately:"
echo "  - Time per sample: ~45-90 minutes (depending on file size)"
echo "  - Total time: ~60-120 minutes (vs 270-540 minutes sequential)"
echo "  - Speedup: ~4.5x faster than original script"
echo ""
echo "To use a different configuration:"
echo "  1. Edit this file and uncomment your preferred option"
echo "  2. Source it: source parallel_config.sh"
echo "  3. Run: ./create_big_wig_parallel.sh" 