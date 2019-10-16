# DengueZikaProteomics

########################

# Required Software for Workflow:
## R (version used: 3.6.0), installed from https://www.r-project.org/
## RStudio (version used: 1.2.1335), installed from https://www.rstudio.com/
## mapDIA (version used: 3.1.0), installed from https://sourceforge.net/projects/mapdia/files/
## WEKA (version used: 3.8), installed from https://www.cs.waikato.ac.nz/ml/weka/downloading.html

# Typical install time on a "normal" desktop computer: seconds
# Software available for Windows, Mac OS X and other platforms. No specific requirements.
# Tested on macOS High Sierra Version 10.13.3.

# Expected run time on a "normal" desktop computer for all steps: seconds

# For Step 1 (Data Post-Processing), a reduced Spectronaut csv output file is used as input. 
# Spectronaut is a commercial software, available at https://biognosys.com/shop/spectronaut

########################

#1 Post_Processing
## Software: RStudio Version 1.2.1335, R 3.6.0
## Script: Script_Data_Post_Processing.R
## Input files: Zika_Dengue_reduced.txt; ConditionSetup.txt; sample_meta_data.txt
## Settings: 0.6 non-missing values, lowess f=0.25
## Expected results: number of proteins: 303, number of peptides: 3329, number of fragments: 12887
## Output files: meta_data.txt; quant_data.txt; quantQC_data.txt

#2 mapDIA
## Software: mapDIA Version 3.1.0
## Input files: quant_data.txt; input.txt
## Output file: protein_level_data_mapDIA.txt (among others)

#3 PCA 
## Software: RStudio Version 1.2.1335, R 3.6.0
## Script: Script_PCA.R
## Input files: protein_level_data_mapDIA.txt; sample_meta_data.txt
## Output file: PrincipalComponent_Protein_Level_ClinicalVar.txt; meta_data_with_PC.txt

#4 Multiple_Linear_Regression
## Software: RStudio Version 1.2.1335, R 3.6.0
## Script: Multiple_Linear_Regression_Script.R
## Input files: protein_level_data_mapDIA.txt; meta_data_with_PC.txt
## Expected results: 26 significant proteins for pval_Zika <= 0.01; 71 significant proteins for pval_Zika <= 0.05
## Output file: pvalue_matrix.txt

#5 WEKA_Machine_Learning
## Software: WEKA 3.8
## Input file: Original_data_all_patients_sig_proteins_0_01_diagnosis.arff (26 significant proteins from step 4)
## class value: diagnosis (Zika/Dengue)
## classifiers used: logistic regression & J48, both with cross-validation (folds 10)
## expected results: for logistic regression, 78.6885 correctly classified for attributes (7 best proteins) P01024, P02671, P02766, P06396, P10720, P36955, Q14624; 75.4098 % for attributes (3 best proteins) P02671, P06396, P10720        
