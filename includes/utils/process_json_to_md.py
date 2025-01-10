import os
import json

# Define the folder containing the JSON files
folder_path = "tests/challenge_instances"
output_md_file = "output_summary.md"

# Define the categories to check
categories = [
    "CONVEX_NO_CONSTRAINTS",
    "CONVEX_OPEN_CONSTRAINTS",
    "CONVEX_CLOSED_CONSTRAINTS",
    "NOT_CONVEX_PARALLEL_NO_CONSTRAINTS",
    "UNSPECIFIED_BOUNDARY"
]

# Create or overwrite the output markdown file
with open(output_md_file, "w") as md_file:
    md_file.write("# JSON Instance Summary\n\n")

    # Iterate through files in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".json"):
            json_file_path = os.path.join(folder_path, file_name)

            # Read the JSON file
            with open(json_file_path, "r") as json_file:
                try:
                    data = json.load(json_file)
                    instance_uid = data.get("instance_uid", "Unknown")
                    boundary_type = data.get("region_boundary", "UNSPECIFIED_BOUNDARY")

                    # Check if boundary type matches one of the categories
                    category_match = boundary_type if boundary_type in categories else "UNSPECIFIED_BOUNDARY"

                    # Write to the Markdown file
                    md_file.write(f"## {file_name}\n")
                    md_file.write(f"- **Instance UID**: {instance_uid}\n")
                    md_file.write(f"- **Boundary Category**: {category_match}\n\n")

                except json.JSONDecodeError:
                    md_file.write(f"## {file_name}\n")
                    md_file.write("- **Error**: Failed to parse JSON file.\n\n")

print(f"Summary written to {output_md_file}")
