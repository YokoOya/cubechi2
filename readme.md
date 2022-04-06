## What is This

This is a program to calculate the chi-squared values between 2 FITS cube files. 
'step2_cubechi2.cpp' is the main script.

The other python scripts can be used if the FITS files to compare are made by using 'FERIA' (https://github.com/YokoOya/FERIA). 
The python scripts work with python2. 


## Contents
- step1_makeInputfile.py
- step2_cubechi2.cpp
- step3-1_plot.py
- step3-2_texTable.py
- template.in



## step1_makeInputfile.py

This script outputs a text file to input to 'step2_cucbechi2'. 
'template.in' is a sample. 

1 observational FITS cube is compared with many FITS cubes of models made by FERIA. 
All the FITS cubes are required to have the same pixel numbers and pixel sizes. 



## step2_cubechi2.cpp

This script calculate the chi-squared values between 2 FITS cube files. 
'template.in' is the sample of the input file. 

'template.in' assumed: 
- The 1st, 2nd, and 3rd axes of the FITS files are Right Ascension (degree), Declination (degree), and the velocity (m/s). 
- The intensity is in Jy/beam. 
The chi-squared values can be calculated even with other manners for the axes and units. 
See below for the details of 'template.in'.



## step3-1_plot.py

This script makes plots from the output files of 'step2_cubechi2'. 


## step3-2_texTable.py

This script makes a tex file with tables from the output files of 'step2_cubechi2'. 
'aastex.cls' is required for typeset. 


## template.in

- Name of the file to be output
- Name of the 1st FITS file to be compared 
- Name of the 2nd FITS file to be compared
- The velocity axis of the 2nd FITS file is shifted by this value when it is compared with the 1st FITS file. (in m/s)
- RMS noise level in the 1st FITS file and threshold level for comparison
- Region to be compared. The minimum and maximum values of Right Ascension, Declination, and the velocity. 
- Fixed weight of the 2nd FITS file to the 1st FITS file. If this value is <= 0, the intensities in the 2nd FITS file are scaled to reduce the chi-squared values, preserving the relative intensity among pixels. 





