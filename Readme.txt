This folder contains the core codes used in the making of 
"Stimulus selection enhances value-modulated somatosensory processing in the superior colliculus"
publication in PLOS BIOL 2025
Any enquiries should be directed to: YunWen Chu at ywchud@gmail.com

Prerequisites:

Kilosort2.5 and Deeplabcut is used for spike sorting and whisker tracking respectively, which will not be provided here. 
Please download their code base beforehand or use other equivalent codes.
This pipeline requires outputs from Kilosort2.5 for sorted spike infos and .csv files from DeeplabCut which contains coordinates of tracked whisker positions.

%----------------------------------------------------------------------
Data:

To demonstrate the entire processing pipeline from collected raw data,
an example of raw superior colliculus spike data along with the mouse whisking videos and tracked coordinates (version 2) is uploaded to Zenodo. This can be downloaded and should be placed in the data folder.
Preprocessed and simplified data of more mice (version 1) can also be accessed freely in the following Zenodo database:
https://zenodo.org/records/13346094?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjJiNzJjOTRjLTAxMWUtNDI1My1hOWM1LWFjZDdhNGVlZmQ4OCIsImRhdGEiOnt9LCJyYW5kb20iOiIxYTdhOTU2MDEzYmFlN2FjNzIxZGIzMWE5YmYyOTlhZiJ9.iShz1EPM-bfR678oTzT8Hu1npSKb5mq_Ird8BQKkKcaZMN-s7kS3aA-7PIwPPpLdIH_veMGVmgsbk7jHWZ4RWA

%----------------------------------------------------------------------
Instructions:

The general pipeline for (1) Preprocessing (2) Feature selection (3) Analysis (4) Visualization is described below:

(1) Preprocessing

Edit the config parameters and run the master script "~\StimulusSelection\Preprocessing\Preprocess_master_script.m" to extract essential digital/analog signals. 
This outputs 'Exp' which contains most of the experiment parameters and
'Sig' which contains the extracted analog/digital signals.

(2) Feature selection

Edit and run the script "~\StimulusSelection\Feature_Selection\FeatureSelection_script" to extract useful behavioral and neural features.
Behavioral features include mouse whisking kinematics (e.g. whisking angle/curvature), touches (single/multi-whisker), runspeed, lick and task performance.
Neural features include spike data and their post-stimulus time histogram (PSTH) alignment with different stimulus conditions.
Edit and run eventTriggered_NeuronPSTH.m for more options

(3) Analysis

This folder contains various functions used in spike analysis, such as current source density (CSD), SVM classification of stimulus prediction from spikes, feature occupancy control etc...

(4) Visualization 

This folder contains various functions used for data visualization such as line plotting with shaded error, raster, video playing etc....