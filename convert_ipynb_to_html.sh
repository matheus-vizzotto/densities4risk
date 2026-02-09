#!/bin/bash

# 1. Check if at least the notebook path was provided
if [ -z "$1" ]; then
    echo "Usage: ./export_nb.sh path/to/notebook.ipynb [custom_name]"
    exit 1
fi

NB_PATH=$1
# 2. Use the second argument if provided, otherwise default to "analysis"
BASE_NAME=${2:-"analysis"}
OUTPUT_DIR="reports/html"

# 3. Generate the timestamp (Format: 20260209)
TIMESTAMP=$(date +%Y%m%d)
FINAL_FILENAME="${BASE_NAME}_${TIMESTAMP}"

mkdir -p $OUTPUT_DIR

echo "Converting $NB_PATH to $OUTPUT_DIR/${FINAL_FILENAME}.html..."

# 4. Run the conversion
jupyter nbconvert --to html \
    --output-dir=$OUTPUT_DIR \
    --output="$FINAL_FILENAME" \
    --template lab \
    "$NB_PATH"

if [ $? -eq 0 ]; then
    echo "------------------------------------------------"
    echo "Success. File available at: $OUTPUT_DIR/${FINAL_FILENAME}.html"
else
    echo "Error: Conversion failed."
fi