# rsfMRI
This is my Matlab implementation of resting-state functional MRI time-series analysis to identify intrinsic brain network

The codes include motion correction, image registration, image segmentation, spatial normalization, image smoothing, time domain 
filtering, nuisance removal, spectral analysis, trend removal, correlation and causality analysis, principal component analysis 
to reduce dimentionality, general linear model and independent component analysis for brain network detection

### Top level
> fmri_preprocess_spm8_ycheng.m

Main function, runs through multiple sessions of one subject in a group folder. The steps and specifications for what is to be done is from the set up file <br />

<br />

> fmri_preprocess_specs.m

Set up file, initializing parameters involved in preprocessing steps to perform and opions for each <br />

<br />

### Some of the core components used in final script:
> fmri_extract_nuisance.m 

Function to extract nuisance from the time series data <br />
<br />

> fmri_regress_nuisance.m, fmri_regress_nuisance_1D

Function to regress out nuisance covariates <br />
<br />

> conn.m

Function to calculate the correlation between time-series data of voxels. <br />

<br />

> fmri_smooth_spm8.m

Function to perform smoothing on functional files in one directory <br />

<br />

> fmri_time_filt.m

Functions to filter time series, including low-pass, high-pass, band-pass, stop-band, and linear detrend filtering






