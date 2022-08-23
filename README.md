# LIMA: Lipid Mediator Analyst Application 

# Overview
This repository contains the code and required files to run the LIMA app locally using RStudio. The LIMA application is a user-friendly application to process and analyse lipid mediator profiling LC-MS data. The application built automatizes the calculation of lipid mediator concentrations and allow the analyses of the data by statistics, machine learning, and pathway analysis. This tool grants biologists and researchers an efficient manner to analyse lipid mediator data which can be used to make biological findings including identifying biomarkers and new therapeutic agonists for diseases

# Requirements
All the R scripts were created and used on **Mac OS Montery version 12.5**:

**R version**: 4.1.3 

**R Studio version**: 2022.07.1+554 

The scripts should be compatible with Windowds and Linux operating systems. 

For installing R and R Studio, follows the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/).

Packages: For the classifyre package, a zip file has been included in the LIMA App folder. The remainng packages' installation has been included in the LIMA_Final.R script.


# Contents

The repository contains two main folders. 

## [LIMA App](https://github.com/ZainabHaybe/Lipid-Mediator-Analyst/tree/main/LIMA%20App)

This folder contains the r script, classifyre package, and the coefficients folder required to run the LIMA App. 

**LIMA_App.R** script: is the r script that once opened on R Studio, the user can click "Run APP" to open and run the application. The script begins with the install packages and load packages for all the packages required for this script. Then it contains the user interface with the scripts for each page. Finally, is the server with the functions that does the processing for each function. 

## [coefficients](https://github.com/ZainabHaybe/Lipid-Mediator-Analyst/tree/main/LIMA%20App/coefficients)

Contains all the files with the additional information needed to process the data and run the application.

**MS2_CrossCoef.csv, MS3_CrossCoef.csv, MS4_CrossCoef.csv**: It is the cross reference coeffieinces which is the signal of the internal standards when the method was validated for the 3 mass spectometry instruments. 

**MS2_ConvCoef.csv, MS3_ConvCoef.csv, MS4_CrossCoef.csv**: It is the standard curve values for the mediators used to convert the area into amount in picograms for each of the 3 mass spectometry instruments. 

**Met_IntStd.csv**: The file containing the internal standards for each of the lipid mediators.

**Pathways2.csv**: The file containing the pathways for each internal standard.

**Rename_Comp.csv**: The correct naming information for some of the mediators during the data processing step.

**Rename_Met2.csv**: The file that standardizes the names of the mediators. 

**SPM_FA.csv**: The Fatty acid family for each of the mediator.





