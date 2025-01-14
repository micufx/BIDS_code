# Motion capture & mobile EEG 

This project involves analyzing motion and EEG data collected during basketball free-throw shooting. The scripts in this repository are designed to process and analyze the data stored in the `data` directory using MATLAB.

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

**NOTE:**: The `fieldtrip` and `eeglab` directories have to be located on the top level directory to be added to the MATLAB path automatically.

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
