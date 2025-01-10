#!/bin/bash

# Directory containing JSON files
input_directory="tests/challenge_instances"

# Output file name (remains constant)
output_file="solution_output.json"

# Loop through each JSON file in the specified directory
for json_file in "$input_directory"/*.json; do
    # Check if the file exists
    if [[ -f "$json_file" ]]; then
        echo "Processing $json_file..."
        # Run the C++ program with the current JSON file as input
        ./opt_triangulation -i "$json_file" -o "$output_file"
    else
        echo "No JSON files found in $input_directory."
    fi
done

echo "All files processed."
