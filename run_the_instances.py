import os
import subprocess

# Define the folder containing the JSON files
folder_path = "tests/challenge_instances"
cpp_executable = "./opt_triangulation"  # Adjust this path if necessary
output_file = "solution_output.json"  # Default output file name

# Get and sort JSON file names in alphabetical order
json_files = sorted(
    [f for f in os.listdir(folder_path) if f.endswith(".json") and f.startswith("ortho")]
)

# Store all running processes to monitor them if needed
processes = []

# Iterate through all JSON files in the folder
for file_name in json_files:
    if file_name.endswith(".json"):
        json_file_path = os.path.join(folder_path, file_name)

        # Start the C++ program asynchronously
        process = subprocess.Popen(
            [cpp_executable, "-i", json_file_path, "-o", output_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Add the process to the list
        processes.append((file_name, process))

# Optionally, wait for all processes to complete
for file_name, process in processes:
    stdout, stderr = process.communicate()
    if stderr:
        print(f"Error for {file_name}: {stderr}")
    else:
        print(f"Completed {file_name}: {stdout}")



