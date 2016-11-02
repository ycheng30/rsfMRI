function fmri_smooth_spm8(func_files,fwhm)
%Function to perform smoothing on functional files in one directory
%Usage
%   fmri_smooth_spm8(func_files,fwhm)
%       func_files cell array of strings
%       fwhm (3 element 1-D array) FWHM of gaussian smoothing kernel

% Ying Cheng, Modified Sep 1 2012 - changed to be spm8 compatible
% Ying Cheng, Oct 13, 2012

if(length(fwhm)==1), fwhm=repmat(fwhm,1,3); end;
matlabbatch{1}.spm.spatial{1}.smooth.fwhm = fwhm;
matlabbatch{1}.spm.spatial{1}.smooth.dtype = 0;
matlabbatch{1}.spm.spatial{1}.smooth.data=func_files;
save(fullfile(fileparts(func_files{1}),'smooth_job.mat'), 'matlabbatch');
spm_jobman('run', matlabbatch);


