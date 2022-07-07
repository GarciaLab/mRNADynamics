# mRNADynamics: Single-cell transcriptional dynamics in living multicellular organisms
LivemRNA is an image processing and transcriptional time-series analysis MATLAB software package. Its main purpose is image analysis of live mRNA data obtained using methods such as MS2 or PP7 in the early embryo of the fruit fly Drosophila melanogaster. This data usually consists of one or more channels of fluorescent puncta (MCP channel), one channel of a nuclear marker such as Histone-RFP (Nuclear channel), and one or more channels of fluorescently tagged transcription factors.

# 0. How to cite this code

Please, when you use this code cite “H. G. Garcia, M. Tikhonov, A. Lin, T. Gregor, Quantitative imaging of transcription in living Drosophila embryos links polymerase activity to patterning. Curr Biol 23, 2140-2145 (2013).”

# 1. Required Matlab toolboxes:

- Image Processing Toolbox
- Computer Vision Systems
- Curve Fitting
- Image Acquisition
- Image Processing
- Mapping
- Optimization (sometimes default but not always)
- Parallel Computing (required for many of our pipeline scripts)

# 2. Installation

The installation scripts will automatically create the folder structure described in the following figure:

![installation-structure](https://github.com/GarciaLab/mRNADynamics/blob/master/doc/installation-structure.jpg?raw=true)

1. Create a folder called “LivemRNA” as shown in the figure above. This will be our main repository of code and data.
2. Clone the  “mRNADynamics” repository from the “Master” branch into it.
3. Go into the “mRNADynamics\src” folder and run “InstallmRNADynamics”.
4. Close and re-start Matlab so that all the changes are implemented. You should see in MATLAB’s console output the message “mRNADynamics Startup script executed”.
5. Inside the “Data” folder one level up from “mRNADynamics” you will find two newfolders: “RawDynamicsData” and “DynamicsResults”. Also note that in “LivemRNA” the file “ComputerFolders.csv” has been created. If you want to change the location of any of the folders just edit this file.

# 1. Introduction and example usage

This text includes instructions for the acquisition and analysis of live mRNA data obtained using methods such as MS2 or PP7 in the early embryo of the fruit fly Drosophila melanogaster. Throughout the code we will assume that the data consists of one channel of fluorescent puncta (MCP channel) and one channel of a nuclear marker such as Histone-RFP (Nuclear channel). The following figure shows these channels as well as the flow of analysis all the way from raw images to data structures that can be used for making plots and analyzing data.

![scripts-example-usage](https://github.com/GarciaLab/mRNADynamics/blob/master/doc/scripts-summary.jpg?raw=true)

