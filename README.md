# IridiumGeomagJerk
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6599767.svg)](https://doi.org/10.5281/zenodo.6599767) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository gives access to the data and scripts used in the analysis of rapid magnetic field changes as observed by the Iridium Satellite Constellation. 

## Please note
1. The manuscript where the data analysis and results are described is currently under review (see citation below). Please communicate with the authors and use discretion when using the materials available in the repository. 
2. The repository will be updated with more interactive features and scripts once the review process is complete. The updated version will enable the reproduction of all figures contained within the manuscript without needing 1) to download all the data, and 2) access to large computing resources.  

## Citing
When using any data, scripts or information from the repository, please cite the following peer-reviewed (in peer-review) papers:
1. Anderson, B. J., Angappan, R., Barik, A., Vines, S. K., Stanley, S., Bernasconi, P. N., et al. (2021). Iridium communications satellite constellation data for study of Earth's magnetic field. Geochemistry, Geophysics, Geosystems, 22, e2020GC009515. https://doi.org/10.1029/2020GC009515 (Source paper where the spherical harmonic coefficients for the geomagnetic field is provided along with the virtual geomagnetic observatory (VGO) data).
2. Angappan, R., Barik, A., Anderson, B. J., Vines, S. K., Stanley, S.(2021). Fast Global Wave Detection in Geomagnetic Jerks with Commercial Satellites. in review. (Source paper for the scripts and data provided in this repository).

## What is in the repository?
The repository contains three directories:
1. A directory with all the data products named, IridiumDataGeomagneticJerks. This directory contains three data products, namely, the spherical harmonic coefficients needed to reconstruct the geomagnetic field maps, the secular variation time series at each VGO, and a database of all jerks detected within the dataset. 
2. A directory with all the python scripts used to analyze the data. Details of the scripts are given below.
3. A directory with the interactive Jupyter notebook and data necessary to reproduce all the figures given in the paper. (will be published once peer review is complete). 

## Using the data and scripts

### The data
The data is presented in either .csv or .mat files. To access the info and headers in the .mat files you can use the following script:

```python
from scipy.io import loadmat
data = loadmat(filename)
locals().update(data)
```

### The scripts
The list below gives the purpose, outputs, and main options of each script:

1. `IridiumGeomagneticFieldDataMaps.py` - reads in the spherical harmonic coefficients and produces the global magnetic field maps. Also outputs the VGO magnetic field data. This script has an option to analyze the data with the m=0 signal removed or included.
2. `AccioSV.py` - reads in the output file from the script above and produces the secular variation (SV) time series at each VGO. A .mat fie containing the SV time series at each VGO is produced. 
3. `SmoothSV_AccioJerk.py` - reads in the output file from the script above and identifies geomagnetic jerks in the dataset. The file outputs the timing and amplitude of jerk at each VGO.
4. `JoinJerkData.py` - reads in the output file from the script above and joins all the different jerks from different VGOs into one dataset seperated by the concavity of the SV (i.e. peak or valley shaped jerk). This script has an option to analyze each different component of the magnetic field. 
5. `JoinPeakValley.py` - reads in the output file from SmoothSV_AccioJerk.py and joins all the different jerks from different VGOs into one dataset without seperating the peaks and valleys as seperate files. The peaks and valleys are still identifiable through the sign of the jerk amplitude. 
6. `OverallScript.sh` - runs through scripts 1 through 5 and produces all their outputs in one go. 

## Interactive Script
To be published once review process is complete.

## Questions
If you have any questions, please feel free to contact the corresponding author, Regupathi (Regu) Angappan at rangapp1@jhu.edu. Regu can also be found on Twitter, @ReguSphere.   
