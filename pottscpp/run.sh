#!/bin/bash

# Check if both C++ executable and Python script names are provided
if [ $# -ne 2 ]; then
  echo "Usage: ./run_scripts.sh <cpp_executable_name> <python_script_name>"
  exit 1
fi

# Extract arguments
CPP_EXEC="$1"
PYTHON_SCRIPT="../python/$2"

# Path to the bin folder for the C++ executable
BIN_DIR="./bin"

# Full path to the C++ executable
CPP_EXEC_PATH="$BIN_DIR/$CPP_EXEC"

# Check if the specified C++ executable exists
if [ ! -f "$CPP_EXEC_PATH" ]; then
  echo "Error: C++ executable '$CPP_EXEC_PATH' not found."
  exit 1
fi

# Check if the Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
  echo "Error: Python script '$PYTHON_SCRIPT' not found."
  exit 1
fi

# Run the C++ executable
echo "Running C++ executable: $CPP_EXEC_PATH"
"$CPP_EXEC_PATH"

# Check if the C++ execution succeeded
if [ $? -ne 0 ]; then
  echo "Error: C++ executable failed to execute."
  exit 1
fi

# Run the Python script
echo "Running Python script: $PYTHON_SCRIPT"
python3 "$PYTHON_SCRIPT"

# Check if the Python script executed successfully
if [ $? -ne 0 ]; then
  echo "Error: Python script failed to execute."
  exit 1
fi

echo "Both scripts executed successfully!"
