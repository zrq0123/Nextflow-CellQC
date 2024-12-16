# Nextflow-CellQC

## Overview
In this study, we deployed nextflow to design a workflow that integrates EmptyDrops and DoubletFinder for the detection of empty droplets and mixed cells.

## Installation
The versions of some packages used in the experiment are as follows:
nextflow version 24.10.0.5928
Seurat==5.1.0
ModEvA==3.2.0
DropletUtils==1.26.0
fields==16.3
parallel==4.4.0
KernSmooth==2.23.24
ROCR==1.0.11
sctransform==0.4.1


## Usage
For EmptyDrops
```bash
nextflow run doubleidentification.nf
--type empty
--file_input
--file_output 
--file_label 
--fdr 0.005


For DoubletFinder
```bash
nextflow run doubleidentification.nf
--type double
--file_input 
--file_output 
