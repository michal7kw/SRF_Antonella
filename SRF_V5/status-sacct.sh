#!/bin/bash

# Get job status from sacct
status=$(sacct -j "$1" --format=State --noheader | head -n 1 | awk '{print $1}')

if [[ $status == "COMPLETED" ]]; then
    echo "success"
elif [[ $status == "FAILED" || $status == "TIMEOUT" || $status == "CANCELLED+" || $status == "NODE_FAIL" || $status == "PREEMPTED" ]]; then
    echo "failed"
else
    echo "running"
fi
