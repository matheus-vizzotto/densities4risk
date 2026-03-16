#!/bin/bash

OUTPUT_DIR="environment"
TIMESTAMP=$(date +%Y%m%d)

mkdir -p $OUTPUT_DIR

FILENAME="environment_${TIMESTAMP}.yml"

echo "Exporting current conda environment..."

conda env export --from-history --no-builds > "$OUTPUT_DIR/$FILENAME"

if [ $? -eq 0 ]; then
    echo "----------------------------------------"
    echo "Environment saved to:"
    echo "$OUTPUT_DIR/$FILENAME"
else
    echo "Error: export failed."
fi
