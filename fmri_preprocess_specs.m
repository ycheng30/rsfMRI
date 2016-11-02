% Preprocessing steps to perform and Options for each

%% Computer capability: Set to 1 if the computer has low RAM available (<2GB)
low_ram = 0;
prefix='';

%% Prepare input files
%PAR/REC to nifti convert
parrec2nii = 0;

%% Set origin - Setting this will make the program wait for user input!
set_origin = 0;

%% Slice time correct
slice_time_correct = 0;
slice_order=3;
    %1=ascending, 
    %2=descending; 
    %3=interleaved ascending starting at 1; 
    %4=interleaved ascending starting at 2;
    %5=interleaved descending starting at end; 
    %6=interleaved descending starting at end-1
    %7=custom file - use slice_order_file variable
slice_order_file = ''; % File for slice order - ignore for 
ref_slice=-1; % -1 will set reference slice as (spatially) middle slice

%% Motion correction
motion_correct = 1;

%% Coregistration
coregister = 1; %1=Anatomical to Functional, 2=Functional to Anatomical

%% Segment the MPRAGE
segment_anat = 1; 

%% Spatial Normalization
normalize = 1; %1=Use unified segmentation (only write), 2=Use EPI template when anatomical is not used (estimate and write)

%% Create brain mask
create_brain_mask_file = 1;
brain_mask_file='brain_mask.nii';

%% High pass time domain filtering
hp_filter = 1;
hp = 0.005; %Cutoff frequency 200s 

%% Smoothing for ICA
ica_smooth = 1;
ica_smooth_fwhm = -1; %-1 will set fwhm to twice the acquisition voxel size

%% Nuisance removal - typically performed only for seed based analysis
nuisance_remove = 0;
nuisance_file_postfix='nuisance';
nr_options.motion_params=1; %1;
nr_options.detrend_motion_params=0;
nr_options.filter_motion_params=1;
nr_options.filter_motion_params_cutoff=0.005;
nr_options.filter_motion_params_band='high';
nr_options.diff_motion_params=1; %Use motion params differential
nr_options.square_motion_params=0; %Use motion params square
nr_options.global_signal=1;
nr_options.wm_signal=2; %1=mean, 2=pca, 0=not included
nr_options.wm_pca_percent=0.96; % percent signal to use when using pca
nr_options.csf_signal=2; %1=mean, 2=pca, 0=not included
nr_options.csf_pca_percent=0.96; % percent signal to use when using pca

%% Smoothing
var_smooth = 0; %For FC analysis
smooth_fwhm = -1; % -1 will set fwhm to twice the acquisition voxel size

%% Band-pass time domain filtering
bp_filter = 1; %For FC analysis
bp = [0.01 0.1]; %Filter cut-offs in Hz
%filter_order=2; %Filter order