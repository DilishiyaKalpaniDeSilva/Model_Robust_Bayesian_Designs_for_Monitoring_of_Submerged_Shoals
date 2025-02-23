#!/bin/bash

# Function to create directory tree structure
create_tree() {
  local dir="$1"
  local indent="$2"
  echo "$indent|-- $(basename "$dir")"
  for item in "$dir"/*; do
    local name
    name=$(basename "$item")
    if [ -d "$item" ]; then
      create_tree "$item" "$indent  |"
    else
      echo "$indent  |-- $name"
    fi
  done
}

# Use the current working directory as the root
root_directory="$(pwd)"

# Define the output file
output_file="directory_tree.txt"

# Call the create_tree function and redirect the output to the file
create_tree "$root_directory" "" > "$output_file"

# Display a message
echo "Directory tree structure with files has been saved to $output_file"
