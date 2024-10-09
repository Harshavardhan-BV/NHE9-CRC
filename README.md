# NHE9-CRC
Data analysis of colorectal cancer related to NHE9 and EMT

Some codes present in other repositories are given as submodules. 
- If you are using the code for the first time, you can clone the repository using the following command:
```bash
git clone --recurse-submodules https://github.com/Harshavardhan-BV/NHE9-CRC.git
```
- If you have already cloned the repository, you can update the submodules using the following command:
```bash
git submodule update --init --recursive
```

## EMT Scoring (76GS, KS, MLR) on Mircoarray data
- The script automatically downloads the data and performs EMT scoring
```bash
cd EMT_score_MicroArray
Rscript ./EMT_GEO.R <"GSEID">
```

## EMT Scoring (76GS, KS, MLR) on RNA-Seq data
```bash
cd EMT_Scoring_RNASeq
```
- The data have to be downloaded and put in Data folder
- Run the script
```bash
Rscript all_GSE_EMT.R
```

## ssGSEA
```bash
cd ssGSEA
```
- Copy the signature gmt file
```bash
cp ../Signatures/NHE9.gmt ./Signatures/
```
- Copy the expression matrix (TPM) to Data folder
```bash
cp /source/path ./Data/
```

## scRNA
```bash
cd scRNA
```
- The data have to be downloaded and put in Data folder
- Perform normalization and imputation on raw data
```bash
python scanpy_magic.py
```
- Perform AUCell scoring
```bash
Rscript AUCell.R
```
- Peform CMS Subtyping
```bash
Rscript CMSCaller.R
```