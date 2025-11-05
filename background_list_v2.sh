#!/bin/bash

INPUT_FILE="datasets_bkg.txt"
OUTFILE="nano_datasets.html"
> "$OUTFILE"

echo "<html><body>" >> "$OUTFILE"

while IFS= read -r ds || [ -n "$ds" ]; do
    [[ -z "$ds" || "$ds" =~ ^# ]] && continue
    miniaod=$(dasgoclient -query="dataset dataset=$ds")

    for md in $miniaod; do
        echo "<h2><b>$md</b></h2>" >> "$OUTFILE"

        nanos=$(dasgoclient -query="child dataset=$md")
        while IFS= read -r nano || [ -n "$nano" ]; do
            [[ -n "$nano" ]] && echo "<p>$nano</p>" >> "$OUTFILE"
        done <<< "$nanos"

        echo "<hr>" >> "$OUTFILE"
    done
done < "$INPUT_FILE"

echo "</body></html>" >> "$OUTFILE"
echo "Done! HTML file saved in $OUTFILE"

