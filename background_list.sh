#!/bin/bash

# Input file containing MiniAOD dataset patterns (one per line)
INPUT_FILE="datasets_bkg.txt"

# Output file
OUTFILE="nano_datasets.txt"
> "$OUTFILE"  # clear previous content

# Read each line in the input file
while IFS= read -r ds || [ -n "$ds" ]; do
    # Skip empty lines or lines starting with #
    [[ -z "$ds" || "$ds" =~ ^# ]] && continue

    echo "Processing: $ds"

    # Query DAS for the actual MiniAOD dataset(s)
    miniaod=$(dasgoclient -query="dataset dataset=$ds")

    for md in $miniaod; do
        echo "MiniAOD: $md"

        # Write MiniAOD line with colon
        echo "$md:" >> "$OUTFILE"

        # Query child NanoAOD dataset(s)
        nanos=$(dasgoclient -query="child dataset=$md")

        # Write each NanoAOD on its own line
        while IFS= read -r nano || [ -n "$nano" ]; do
            echo "$nano" >> "$OUTFILE"
        done <<< "$nanos"

        # Add empty line for separation (optional)
        echo "" >> "$OUTFILE"
    done

done < "$INPUT_FILE"

echo "Done! Results saved in $OUTFILE"

