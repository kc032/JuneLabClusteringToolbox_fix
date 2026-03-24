# Clustering toolbox
## Design
The Clustering Toolbox available here was designed as a user interface that allows for clustering analysis of untargeted metabolomics data. This user interface contains 23 functionalities spanning from hierarchical clustering to clustering comparision to KEGG pathways analysis. Most functionalities within this UI contain user options allowing for customization of analyses. Furthermore, we provide a start-up window that allows the user to define the number of threads (e.g., processes they would like to use during analysis allowing for rapid clustering analysis). We recommend the user use at least one less thread than available on the computer, with the optimal being half of the available threads. 

## Ensemble Clustering combined with Clustering Optimization (ECCO)
Ensemble clustering combined with Clustering Optimization or ECCO is the main functionality within the clustering toolbox. This form of ensemble clustering is discussed in our pre-print [pre-print](https://www.biorxiv.org/content/10.1101/2022.11.03.515009v1.abstract) and can be adapted for your preferences. ECCO is built for ensemble clustering solutions of agglomerative hierarcical clustering algorithms. 

## Data Pre-processing notes
The Clustering Toolbox was built to mimic the pre-processing steps taken during untargeted metabolomics data analysis. We acknowledge that we do not provide all of the available pre-processing steps and recommend pre-processing prior to submission to UI and selecting 'None' and 'None' when prompted to select a transformation and scaling for your data. 

## Installation and set-up

### Option A — Conda (recommended, especially on Windows)

The file `environment.yml` defines the **ECCO_env** environment (Python 3.11–3.14, conda-forge). It matches what the GUI actually imports (no TensorFlow/Jupyter bloat).

1. Clone or download this repository.
2. From the project directory:
   - `conda env create -f environment.yml` or upload directly to Anaconda interface
   - `conda activate ECCO_env` or open Anaconda prompt with the ECCO_env selected
3. Navigate to the code location and start the UI: `python JuneLabClusteringGUI.py`

To share the environment on **Anaconda.org**: build the env locally, then `conda env export` (optionally from a clean env) and upload the YAML, or publish a package that documents `conda env create -f environment.yml`.

To refresh after `environment.yml` changes: `conda env update -f environment.yml --prune`.

If conda-forge does not yet offer your desired Python (e.g. 3.14) on your platform, edit `environment.yml` and set `python=3.12` or `python=3.11`, then run `conda env create` again.

### Option B — pip only

Use a recent Python (3.11+). Install dependencies with:

- `pip install -r requirements.txt`

`requirements.txt` uses **fpdf2** (PyPI name for the maintained `fpdf` API). Start the UI with `python JuneLabClusteringGUI.py`.

**NOTE:** If Excel files fail to load, ensure **openpyxl** (and **xlrd** for older `.xls` if needed) are installed; they are listed in `environment.yml` and `requirements.txt`.

## Example Files
The ExampleFiles directory contains **two example input files**. Please see further documentation for other input files. 

## Troubleshooting
Please submit an Issue in the Issues tab, and I will address as quickly as possible. 

