# Estimating the Dynamic Properties of Wall-Frame Structures

This repository includes additional data, Python scripts, and files used in the journal paper **"Estimating the Dynamic Properties of Wall-Frame Structures"**. The data is available for all to download and use.

## Available Datasets

Three datasets are provided:

1. **Regular Wall-Frame Structures Dataset**:  
   A dataset of 110 regular wall-frame structures sized and assessed in the paper. It includes the dimensions of beams, columns, walls, the number of stories, and bays.
   
2. **Partial Height Wall-Frame Structures Dataset**:  
   This dataset includes wall-frame structures with partial height walls.
   
3. **Irregular Storey Heights Wall-Frame Structures Dataset**:  
   The final dataset consists of wall-frame structures with irregular storey heights.

These datasets are provided for the purposes of reproducibility and for users proposing alternative procedures to estimate the dynamic properties of wall-frame structures using a common dataset.

## MathCAD File

A MathCAD file is available, providing a worked example of the procedure outlined in the paper. The values obtained from the approximate method can be compared with the estimated period and mode shape from the Wall-frame structural model developed in OpenSees.

- Download MathCAD software: [MathCAD Official Site](https://www.mathcad.com/en)
- A PDF version of the MathCAD file is also available in this repository.

## Python Scripts

Two Python scripts are included:

1. **Approximate Procedure Script**:  
   This script uses the approximate procedure outlined in the paper (and also provided in the MathCAD file). It takes input structure dimensions from an Excel spreadsheet and estimates the period and mode shape of the structure.

2. **OpenSeesPy Structural Model Script**:  
   This script uses the structural model developed in OpenSeesPy to obtain the period and mode shape of the structure. It can be used to compare results with the approximate method.
