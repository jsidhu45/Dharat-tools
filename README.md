# Dharat-tools

The word Dharat (ਧਰਤਿ) is a Gurmukhi (Punjabi) word, which means Earth/Land. 

This repository is for develpoing new tools to study Earth surface and its process. 
In its current release a new method is developed to estimate the depth of the landslide slip surface and its volume estimation using Spline interpolaton method. The method can be provided with additonal 1D borehole or 3D surface data to futher improve the results. In future, further improvements and new tools will be continuously added with widerange range of applicability to understand complex Earth Science challanges. 

Steps to follow:

Step 1. Download the file with the name Dharat-tools.py and save it in your system drive where you want to run.

Step 2. Install essential python libraries (numpy, pandas, matplotlib, gdal, scipy).

Step 3. Create a new folder in current directory with name (input_data) containing DEM and landslide boundary shape file.

Step 4. Run Python Dharat-tools.py and provide required information in command line.

Step 5. Wait for code to execute and final interpolated results will be saved in results directory.

**Initially run with smaller iterations and coarser DEM resolution to check computation time**

Currenty the code is in raw form, in future improvements will be updated to optimize and speed up the code. 

"THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK."

For more information contact: jaspreet_singh_1@sfu.ca

Singh, J and Sepúlveda, S., 2025. Assessing the landslide failure surface depth and volume: A new spline interpolation method. Engineering Geology. DOI: 10.1016/j.enggeo.2025.108319 
