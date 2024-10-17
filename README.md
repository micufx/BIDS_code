# Mobile Basketball EEG Project

This project involves analyzing EEG data collected during basketball activities. The scripts in this repository are designed to process and analyze the data stored in the `data` directory using MATLAB.

## Directory Structure

- `data/raw`: Contains all the raw `.xdf` files and `.mat` event files.
- `scripts/`: Contains MATLAB scripts for data processing and analysis.

## Usage

Each script in the `scripts` directory accesses files from the `data` directory. Ensure that the `data` directory is properly populated with the necessary files before running any scripts.

## Requirements

- MATLAB
- fieldtrip toolbox (fieldtrip-20230215)
- EEGLAB toolbox (2023.0)
- xdfimport extension for eeglab
- dipfit extension for eeglab

## Running the Scripts

1. Open MATLAB.
2. Navigate to the `scripts` directory.
3. Run the desired script.

Example:
```matlab
cd('path/to/scripts');
run('example_script.m');
```

## License

This project is licensed under the MIT License.
