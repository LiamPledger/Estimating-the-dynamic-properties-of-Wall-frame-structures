# Estimating-the-dynamic-properties-of-Wall-frame-structures
This repository includes additional data, python scripts, and files used in the journal paper "Estimating the dynamic properties of wall-frame structures". The data is available for all to download and use.

Three datasets are available. The first is a dataset of the 110 regular wall-frame structures sized and assessed in the paper. The dataset includes the dimensions of the beams, columns and walls, as well as the no. of stories and bays. The second dataset is for the wall-frame structures with partial height walls. The same data is provided in this spreadsheet. The final dataset consists of the wall-frame structures assessed in the study with irregular storey heights. 

A MathCAD file is available for users. The MathCAD file provides a worked example of the procedure outlined in the paper. The values obtained from the approximate method can also be compared against the estimated period and mode-shape from the Wall-frame structural model developed in OpenSees. MathCAD software can be downloaded at the following link : https://www.mathcad.com/en. A PDF version of the MathCAD file is available as well. 

Two python scripts have been developed.
The first uses the approximate procedure outlined in the paper. This script takes the input structure dimensions from an excel spreadsheet and estimates the period and mode-shape of the structure. 
The second script is the structural model developed in OpenSeesPy which can be used by the user to obtain the period and mode-shape of the structure and compare with the approximate method. 


