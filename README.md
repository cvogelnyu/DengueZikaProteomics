# Dengue Zika Proteomics

# Required Software for Workflow:
R (version used: 3.6.0), installed from https://www.r-project.org/
RStudio (version used: 1.2.1335), installed from https://www.rstudio.com/
mapDIA (version used: 3.1.0), installed from https://sourceforge.net/projects/mapdia/files/
WEKA (version used: 3.8), installed from https://www.cs.waikato.ac.nz/ml/weka/downloading.html

Typical install time on a "normal" desktop computer: seconds
Software available for Windows, Mac OS X and other platforms. No specific requirements.
Tested on macOS High Sierra Version 10.13.3.

Expected run time on a "normal" desktop computer for all post processing and data analysis steps: seconds to minutes

For Step 1 (Data Post-Processing), a Spectronaut output file is used as input. Spectronaut is a commercial software package aimed at analyzing data-independent acquisition (DIA) proteomics experiments, available at https://biognosys.com/shop/spectronaut. Run time for DIA analysis is several hours on a "normal" desktop computer.


# 1 Post_Processing
## Software: 
RStudio Version 1.2.1335, R 3.6.0
## Script: 
Script_Data_Post_Processing.R
## Input files: 
Spectronaut_Output_Zika_BGS_Report_sample.csv: DIA proteomics data. This is a sample Spectronaut output file, since the original file is too large to upload on GitHub (4.43GB). It includes a random sample of 1000 rows from the original BGS report. 
ConditionSetup_meta.txt: Spectronaut run condition setup with diagnosis information for samples
merged_table_input.txt: reduced DIA proteomics data. From line 75 on, the complete file merged_table_input.txt is used instead of the sampled version.
## Expected results: 
number of proteins: 303
number of peptides: 3329
number of fragments: 12887
## Output files:
merged_table_input_sample.txt: sample version of the reduced DIA proteomics data 
quant_data.txt: normalized and quality controlled proteomics data on fragment level

# 2 mapDIA
## Software: 
mapDIA Version 3.1.0
## Input files: 
quant_data.txt: normalized and quality controlled proteomics data on fragment level, output from step 1
input.txt: script for mapDIA
mapdia: mapDIA unix executable file
## Output file: 
protein_level_data_mapDIA.txt: proteomics data on protein level

# 3 PCA 
## Software: 
RStudio Version 1.2.1335, R 3.6.0
## Script: 
Script_PCA.R: R script for PCA on protein level
## Input files: 
protein_level_data_mapDIA.txt: proteomics data on protein level, output from step 2
sample_meta_data.txt: metadata for serum samples
## Expected results:
values of first 5 principal components: 21.40  8.03  6.24  4.10  3.09
## Output file: 
PrincipalComponent_Protein_Level_ClinicalVar.txt: p-values, principal components correlated with metadata in linear regression 
meta_data_with_PC_protein_level.txt: metadata with principal components from PCA on protein level

# 4 Multiple_Linear_Regression
## Software: 
RStudio Version 1.2.1335, R 3.6.0
## Scripts: 
Linear_Regression_Diagnosis_Script.R: R script for linear regression (Zika/Dengue diagnosis)
Linear_Regression_Ambiguity_Script.R: R script for linear regression (ambiguous/unambiguous diagnosis)
## Input files: 
protein_level_data_mapDIA.txt: proteomics data on protein level
sample_meta_data.txt: metadata for serum samples
Test_Samples.csv: samples for test set
sample_meta_data_ambiguity.txt: metadata as above, but with added column on ambiguity based on additional serology results
## Expected results (Diagnosis): 
Analysis A: 13 significant proteins for pval_Zika < 0.01 (3 after BH correction)
Analysis B: 26 significant proteins for pval_Zika < 0.01 (11 after BH correction) 
## Expected results (Ambiguity): 
Analysis A: 3 significant proteins for pval_unambiguous < 0.01 (0 after BH correction); 
Analysis B: 6 significant proteins for pval_unambiguous < 0.01 (0 after BH correction) 
## Output files (Diagnosis): 
pvalue_matrix_AnalysisA_Diagnosis_Training.txt: p-values from linear regression (Analysis A, whole training set)
pvalue_matrix_AnalysisB_Diagnosis_Training.txt: p-values from multiple linear regression (Analysis B, subset of training set, samples with complete metadata)
## Output files (Ambiguity): 
pvalue_matrix_AnalysisA_Ambiguity_Training.txt: p-values from linear regression (Analysis A, whole training set)
pvalue_matrix_AnalysisB_Ambiguity_Training.txt: p-values from multiple linear regression (Analysis B, subset of training set, samples with complete metadata)


# 5 WEKA_Machine_Learning
## Software: 
WEKA 3.8
## Input file: 
Protein_data_mapDIA_norm_Trainingsset.arff: normalized map DIA output (protein data) for training set (82 samples)
Protein_data_mapDIA_norm_Testset.arff: normalized map DIA output (protein data) for test set (40 samples)
Logistic_Regression_3BestProteins.model: logistic regression model file
## class value: 
diagnosis (Zika/Dengue)
## attributes: 
P02671
P02775
P10720
## classifier: 
logistic regression, -R 6.0 -M -1 -num-decimal-places 4
## expected results: 
correctly classified instances training set 76.8293%, test set 77.5%
## output file:
Predictions_WEKA_Logistic_Regression_3BestProteins.csv: WEKA predictions for training and test set
Threshold_Zika_3proteins_test.arff: WEKA threshold curve data Zika 
Threshold_Dengue_3proteins_test.arff: WEKA threshold curve data Dengue 
WEKA_Run_Information.txt: WEKA run information data     
