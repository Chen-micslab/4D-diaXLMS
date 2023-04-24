# 4D-diaXL
## Introduction
4D-diaXL is a workflow which could allow the diaPASEF analysis on proteome-wide cross-linking study. 
## Guide to use 4D-diaXL
### 1. Analysis ddaPASEF data by pLink2 software
The DIA analysis in 4D-diaXL is library-based search, so you need to establish a experimental spectrum library on ddaPASEF mode.
#### a. Convert the .d file of ddaPASEF data to mgf file using DataAnalysis software (Bruker)
#### b. Process the mgf file with 'process_mgf_file.py'
Example:
```
python process_mgf_file.py --filedir './data/example.mgf'  --filename 'example'
```
It will generate a 'example_plink.mgf' which could be processed by pLink2 software.
#### c. Analyze the mgf file with pLink2. 
The 'example_plink.mgf' files were imported to pLink2(2.3.11) for database search.
#### d. Process the cross-linked results of pLink2 with 'process_plink_results.py'
Example:
```
python process_plink_results.py --inputdir './data/example.csv' --outputdir './data/example'
```
It will generate a 'example_crosslink_filter.csv' file.
### 2. Generate 4D crosslinking library
Here you could generate a 4D crosslinking library based on the 'example_plink.mgf' and 'example_crosslink_filter.csv' by running 'generate_4D_library.py'
Example:
```
python generate_4D_library.py --resultsdir './data/example_crosslink_filter.csv' --mgfdir './data/example_plink.mgf' --crosslinker 'DSS'
```
It will generate two library files 'example_crosslink_filter_DIANN_lib.csv' and 'example_crosslink_filter_normal_lib.csv', the 'example_crosslink_filter_DIANN_lib.csv' file could be directly used by DIA-NN software.
### 3. Merge multiple libraries 
If you have several fractions, there may be some overlap identification between multiple fractions, so it is necessary to remove the overlap. You can merge these 'normal_lib.csv' into one csv file, and then run the 'filter_library.py'
Example:
```
python filter_library.py ----filedir './data/merge.csv' 
```
It will generate two library files 'merge_filter_DIANN_lib.csv' and 'merge_filter_normal_lib.csv', the 'merge_filter_DIANN_lib.csv' file could be directly used by DIA-NN software.
### 4. Analysis diaPASEF data by DIA-NN software with 4D cross-linking library
DIA-NN (version 1.8.1) was used to process raw (.d) file using the 4D cross-linking spectral library. FASTA file was not selected. Other parameters was set as follows: ms1 accuracy, 15 ppm; ms2 accuracy, 15 ppm; precursor FDR, 1%; neural network class, single-pass model; quantification strategy, robust LC; cross-run normalisation, RT-dependent.
