# Dharat-tools

The word Dharat (ਧਰਤਿ) is a Gurmukhi (Punjabi) word, which means Earth/Land. 

This repository is for develpoing new tools to study natural hazards such as landslide and debris flow. In its current release a new method is developed to estimate the depth of the landslide slip surface and its volume estimation using Spline interpolaton method. The method can be provided with additonal 1D borehole or 3D surface data to futher improve the results. In future, further improvements and new tools will be continuously added with widerange range of applicability to understand complex Earth Science challanges. 

usage: Dharat-tools.py [-h] [--fill] --dem DEM --resolution RESOLUTION --boundary BOUNDARY [--save-interval SAVE_INTERVAL]
                       [--stop-delta-vol STOP_DELTA_VOL | --stop-max-depth STOP_MAX_DEPTH]
Steps to follow:

Step 1. Download the file with the name Dharat-tools.py and save it in your system drive where you want to run.

Step 2. Install essential python libraries (numpy, pandas, matplotlib, gdal, scipy).

Step 3. Create a new folder in current directory with name (input_data) containing DEM and landslide boundary shape file.

Step 4. Run Python Dharat-tools.py -h for help and instructions to use this code.

Any future changes will be updated in the code.

"THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK."

For more information contact: jaspreet_singh_1@sfu.ca

Reference: Singh, J and Sepúlveda, S., 2025. Assessing the landslide failure surface depth and volume: A new spline interpolation method. Engineering Geology. DOI: 10.1016/j.enggeo.2025.108319 
