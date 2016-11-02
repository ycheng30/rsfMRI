function fmri_extract_nuisance(func_files,tr,rp_file,wm_file,csf_file,brain_mask_file,options,nuisance_file)
%Function to extract nuisance time courses 
%CSF and WM voxel timecourses from the data are extracted and PCA is used to
%reduce dimensions. These PCA components are combined with the motion 
%parameters and global signal (as chosen) to produce a combined nuisance matrix 
%
%Usage
%   fmri_extract_nuisance(func_files,tr,rp_file,wm_file,csf_file,hp,brain_mask_file, ...
%                         options,nuisance_file)
%Input variables
%   func_files - cell array of functional files (from one run)
%   tr         - Repetition time (time taken for one dynamic(or volume) acquisition
%   rp_file    - realignment parameters (file with columns of realignment 
%                parameters (such as what is obtained from spm_realign or mcflirt)
%   wm_file    - White matter segmentation result in same space as functional 
%                (usually got from spm_segment or FSL's fast)
%   csf_file   - CSF segmentation result in same space as functional
%   hp         - is the filter cut-off in Hz (This is needed only in case the data was filtered
%                along time between realignment and motion extraction (we assume
%                butterworth 2nd order high-pass filter).
%   brain_mask_file - filename of the brain_mask image (usually got from
%                fmri_create_brain_mask)
%   options    - a structure with fields listed below which specify what
%                nuisances to use (skip those that are not needed - atleast
%                one option is necessary!)
%                options.motion_params=1; 0 or 1 - to specify if motion parameters should be considered as nuisance
%                options.filter_motion_params=0; 0 or 1 - to specify if
%                   motion_paraters need to be filtered - this is needed
%                   only in case the data is filtered along time between
%                   realignment and nuisance extraction/removal
%                options.filter_motion_params_cutoff=0.005; Cut-off
%                   frequency in Hertz for filtering (use the same cut-off
%                   used for the data)
%                options.filter_motion_params_band='high'; 'high' or 
%                   'low' or 'band' or 'stop' - use the same as for the data
%                options.detrend_motion_params=0; 0 or 1 - to speficy if
%                   motion parameters need to be detrended - this is needed
%                   only in case the data is detrended along time between
%                   realignment and nuisance extraction/removal
%                options.diff_motion_params=1; %Use motion params differential?
%                options.square_motion_params=0; %Use motion params square?
%                options.global_signal=1; %Use global signal?
%                options.wm_signal=2; %1=mean of WM regions, 2=PCA of WM voxel timecourses
%                options.wm_pca_percent=0.95; % percent signal to use when using pca on WM
%                options.csf_signal=2; %1=mean of CSF region, 2=PCA of CSF voxel timecourses
%                options.csf_pca_percent=0.95; % percent signal to use when using pca
%   nuisance_file - the output filename where the nuisance mat file will be written
%               The mat file has a structure variable nui. nui has two fields
%               nui.names - are the names of the nuisances
%               nui.tc - the matrix of timecourses of the nuisances

%Ying Cheng Aug 22, 2013 - Added a few more options - detrend mp,
%   ability to not specify fields in option, if it is not needed
%Ying Cheng Jul 28, 2013 - Cleaned code and commented v
%Ying Cheng May 20, 2013 - Toolbox (fmri_preproc_toolbox) compatibility
%Ying Cheng, Mar 7 2013 - added options (select which nuisances to extract)
%Ying Cheng, Mar 2, 2012 - cleaned up code, no major changes. Removed creation of brain mask (reads already created file).
%Ying Cheng, Dec 21, 2012 - check if CSF (or WM) mask is empty - ignore
%   empty masks and write nuisance names
%Ying Cheng, Dec 21, 2012 - Changed to 4D and separate CSF and WM PCA;
%   changed CSF threshold to 99% instead of 99.5%
%Ying Cheng, Aug 23, 2012 - Added line to create a brain mask; changed the way the script looks for tissue segmentation files
%Ying Cheng, Aug 19, 2012 - Changed for using in RS_LDDMM study
%Ying Cheng, July 10, 2012 - Back up for Change orient, COMPCOR MASKS (Erosion and Thresholding change.)
%Ying Cheng, May 26, 2012 - Principal Components explaining 95% variance
%Ying Cheng, March 05, 2012 - Added PCA and Detrending and Combined all the Nuisance Parameters

%% Read Functional data filenames in func_dir

disp(['Extracting Nuisances: ',fileparts(func_files{1})]);

P=strvcat(func_files);%#ok
so = get_orient(P(1,:));
nui.names=[];
nui.tc=[];

%% Motion parameters
if(isfield(options,'motion_params') && options.motion_params==1),
    nm=load(rp_file); %Read motion file
    if(isfield(options,'detrend_motion_params') && options.detrend_motion_params==1),
        nm=detrend(nm); %detrend motion parameters if specified so
    end;
    if(isfield(options,'filter_motion_params') && options.filter_motion_params==1)
        [b,a]=butter(2,2*tr*options.filter_motion_params_cutoff, ...
            options.filter_motion_params_band); %High pass filter if specified so
        nm=filtfilt(b,a,nm);
    end
    nm=detrend(nm,'constant'); %demeaning - done in all cases
    nui.names={'mp_x','mp_y','mp_z','mp_yaw','mp_pitch','mp_roll'}; %names for the nuisance vectors
    nui.tc=nm;
    if(isfield(options,'diff_motion_params') && options.diff_motion_params==1),
        nui.tc=[nui.tc, [0,0,0,0,0,0; diff(nui.tc)]];
        n=1;for i=(length(nui.names)+1):(length(nui.names)+6),
            nui.names{i}=['d',nui.names{n}];n=n+1; end;
    end;
    if(isfield(options,'square_motion_params') && options.square_motion_params==1),
        nui.tc=[nui.tc,nm.^2];
        n=1;for i=(length(nui.names)+1):(length(nui.names)+6),
            nui.names{i}=['s',nui.names{n}];n=n+1; end;
    end;
end; %End motion params



%% Read White matter segmentation result and make eroded mask
if(isfield(options,'wm_signal') && options.wm_signal~=0),%Make mask for either mean or PCA
    if(~strcmp(get_orient(wm_file),so)),% If reorienting is necessary
        warning('Mywarn:CheckOrientation','Orientations of the Functional and WM mask are not the same');
        change_orient(wm_file,so);
    end;
    Vm = spm_vol(wm_file);
    wm_mask = spm_read_vols(Vm);
    wm_mask = (wm_mask/max(wm_mask(:)))>0.99;
    wm_mask = double(imerode(logical(wm_mask),ones(3,3,3)));
    if(~any(wm_mask)), %IF WM mask is empty force ignore options.wm_signal
        options.wm_signal=0;
        warning('MyWarn:MaskEmpty','WM mask is empty');
    end;
    % Write out the eroded WM mask
    Vm.fname = fullfile(fileparts(Vm.fname),'wm_erode_mask.img');
    Vm.private.dat.fname = Vm.fname;
    spm_write_vol(Vm,wm_mask);
    clear Pm Vm
end;

%% Read CSF segmentation result and make eroded mask
if(isfield(options,'csf_signal') && options.csf_signal~=0), %Make mask for either mean or PCA
    if(~strcmp(get_orient(csf_file),so)),% If reorienting is necessary,
        warning('Mywarn:CheckOrientation','Orientations of the Functional and CSF mask are not the same');
        change_orient(csf_file,so);
    end;
    Vm = spm_vol(csf_file);
    csf_mask = spm_read_vols(Vm);
    csf_mask = (csf_mask/max(csf_mask(:)))>0.99;% Threshold the CSF mask at 99.5%
    csf_mask = csf_mask-bwperim(csf_mask,6); %Erode
    if(~any(csf_mask)), %If CSF mask is empty force ignore options.csf_signal
        options.csf_signal=0;
        warning('MyWarn:MaskEmpty','CSF mask is empty');
    end;
    % Write out eroded CSF mask
    Vm.fname = fullfile(fileparts(Vm.fname),'csf_erode_mask.img');
    Vm.private.dat.fname = Vm.fname;
    spm_write_vol(Vm,csf_mask);
    clear Vm;
end;


if((isfield(options,'wm_signal') && options.wm_signal ~=0) || ...
        (isfield(options','csf_signal') && options.csf_signal ~= 0) || ...
        (isfield(options,'global_sigal') && options.global_signal ~= 0))
    %% Prepare for PCA of WM + CSF
    V = spm_vol(P); %Read the functional images headers

    %% Read brain_mask image
    Vm=spm_vol(brain_mask_file);
    brain_mask=logical(spm_read_vols(Vm));
    clear Vm;

    %% Extract Nuisance Time courses from Data
    Y=spm_read_vols(V);
    if(options.wm_signal),
        tc_wm=zeros(size(Y,4),sum(wm_mask(:)));
        wm_mask=logical(wm_mask);
    end;
    if(options.csf_signal),
        tc_csf=zeros(size(Y,4),sum(csf_mask(:)));
        csf_mask=logical(csf_mask);
    end;
    if(options.global_signal),
        mY=zeros(size(Y,4),1);
    end;
    for i=1:size(Y,4),
        tY=Y(:,:,:,i);
        tY(isnan(tY))=0;
        if(options.global_signal~=0), mY(i) = mean(tY(brain_mask))'; end;
        if(options.wm_signal), tc_wm(i,:)=tY(wm_mask); end;
        if(options.csf_signal), tc_csf(i,:)=tY(csf_mask); end;
    end;
    if(options.global_signal), mY=detrend(mY,'constant'); end; %demean
    clear tY;

    %% Compute Mean/PCA of the WM and CSF separately
    %% WM signals
    if(options.wm_signal > 0)
        tc_wm = detrend(tc_wm,'constant'); %demean
        if(options.wm_signal==1) %Mean
            npca_wm=mean(tc_wm,2);
        else %PCA
            [u,s] = svd(tc_wm*tc_wm');
            eigen_values = diag(s);
            ss = norm(eigen_values)^2;
            ncomp=find(cumsum(eigen_values.^2)/ss > options.wm_pca_percent, 1, 'first' );
            npca_wm = u(:,1:ncomp);
            clear u s ncomp;
        end;
        nui.tc=[nui.tc,npca_wm];
        for i=1:size(npca_wm,2)
            nui.names{1+length(nui.names)}=['wm',num2str(i)];
        end;
    end;


    %% CSF signals
    if(options.csf_signal > 0)
        tc_csf = detrend(tc_csf,'constant'); %demean
        if(options.csf_signal==1), %Mean
            npca_csf=mean(tc_csf,2);
        else % PCA
            [u,s] = svd(tc_csf*tc_csf');
            eigen_values = diag(s);
            ss = norm(eigen_values)^2;
            ncomp=find(cumsum(eigen_values.^2)/ss >  options.wm_pca_percent, 1, 'first' );
            npca_csf = u(:,1:ncomp);
            clear u s ncomp
        end;
        nui.tc=[nui.tc,npca_csf];
        for i=1:size(npca_csf,2)
            nui.names{1+length(nui.names)}=['csf',num2str(i)];
        end;
    end;

    %% Global signal
    if(isfield(options,'global_signal') && options.global_signal==1),
        nui.tc=[nui.tc,mY];
        nui.names{1+length(nui.names)}='global';
    end;

end; %End wm || csf || gs ~= 0

%% Intensity Normalize the nuisance regressors (otherwise matrix inversion fails in Matlab
for i=1:size(nui.tc,2),
    nui.tc(:,i)=nui.tc(:,i)./max(squeeze(nui.tc(:,i)));
end;

%% Save the file
save(nuisance_file,'nui');