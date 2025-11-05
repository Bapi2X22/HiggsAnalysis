#!/bin/bash

# Input file: list of datasets (one per line)
DATASETS_FILE="Nano_bkg.txt"

total_size=0
total_events=0

printf "%-100s %15s %15s %15s %15s\n" "Dataset" "Events" "Size(GB)" "Files" "Size/Event(kB)"
echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------"

while IFS= read -r ds; do
    [[ -z "$ds" ]] && continue  # skip empty lines
    summary=$(dasgoclient -query="summary dataset=${ds}")
    nevents=$(echo "$summary" | grep -o '"nevents":[0-9]*' | cut -d: -f2)
    fsize=$(echo "$summary" | grep -o '"file_size":[0-9]*' | cut -d: -f2)
    nfiles=$(echo "$summary" | grep -o '"nfiles":[0-9]*' | cut -d: -f2)

    # Convert to human units
    size_gb=$(python3 -c "print(${fsize}/1e9)")
    size_per_evt=$(python3 -c "print((${fsize}/${nevents})/1024)" 2>/dev/null || echo 0)

    printf "%-100s %15d %15.2f %15d %15.2f\n" "$ds" "$nevents" "$size_gb" "$nfiles" "$size_per_evt"

    total_size=$((total_size + fsize))
    total_events=$((total_events + nevents))
done < "$DATASETS_FILE"

# Print totals
total_size_gb=$(python3 -c "print(${total_size}/1e9)")
avg_size_per_evt=$(python3 -c "print((${total_size}/${total_events})/1024)" 2>/dev/null || echo 0)

echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------"
printf "%-100s %15d %15.2f %15s %15.2f\n" "TOTAL" "$total_events" "$total_size_gb" "-" "$avg_size_per_evt"

