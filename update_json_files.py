import os
import json

#Path to the directory containing JSON files
directory = "tests/challenge_instances"

#Define configurations
configurations = {
    "local": {
        "method": "local",
        "parameters": {"L": 200, "alpha": 2.4, "beta": 0.2},
        "delaunay": True
    },
    "sa": {
        "method": "sa",
        "parameters": {"alpha": 2.4, "beta": 0.2, "L": 800, "batch_size": 6},
        "delaunay": True
    },
    "ant": {
        "method": "ant",
        "parameters": {"alpha": 2.4, "beta": 0.2, "xi": 1.0, "psi": 2.0, "lambda": 0.2, "kappa": 10, "L": 150},
        "delaunay": True
    },
    "auto": {
        "method": "auto",
        "parameters": {"alpha": 2.4, "beta": 0.2, "L": 800, "batch_size": 6},
        "delaunay": True
    }
}

#Select which configuration to apply (choose from: "local", "sa", "ant")
selected_configuration = "ant"  # Change this to "local" or "ant" as needed

#Validate selection
if selected_configuration not in configurations:
    raise ValueError(f"Invalid configuration selected: {selected_configuration}")

#Fields to add/update
new_fields = configurations[selected_configuration]

#Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".json"):
        filepath = os.path.join(directory, filename)
        
        #Read the JSON file
        with open(filepath, "r", encoding="utf-8") as file:
            data = json.load(file)
        
        #Update fields
        data.update(new_fields)
        
        #Write back to the same file
        with open(filepath, "w", encoding="utf-8") as file:
            json.dump(data, file, indent=4)
        
        print(f"Updated: {filename} with method: {selected_configuration}")
