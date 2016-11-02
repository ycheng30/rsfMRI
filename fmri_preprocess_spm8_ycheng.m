function fmri_preprocess_spm8_ycheng(setup_file)

%Function to preprocess a single subject with multiple fmri runs.
%The steps and specifications for what is to be done is from the set up
%file (created by sin_sub_mul_sess_preprocess_setup)

% modified/updated by Ying Cheng 2013
% this function now works for multiple runs.

%% Initializing

load(setup_file);

disp(['Resting State Preprocessing on: ',fileparts(func_files{1}{1})]);


%% PAR/REC to NIFTI using r2agui -- did not check for multiple runs
if(parrec2nii)
    [p,f,e]=fileparts(anat_file);
    if(~strcmp(e,'.rec'))
        error('Please input rec file to convert PAR/REC');
    end;
    anat_file=rec_convert(anat_file,3);
    anat_file=anat_file{1};
    
    for ir=1:length(func_files),
        [p,f,e]=fileparts(func_files{ir}{1});
        if(strcmp(e,'.rec'))
            func_files{ir}=rec_convert(func_files{ir}{1},3);
        end;
    end;
end;

%% Set origin -- did not check for multiple runs
if(set_origin)
    %Anat file
    spm_image('init',anat_file);
    uiwait;
    
    
    
    %Func file
    global out_transform_file; %#ok I know this is bad programming but don't have the time to fix it!
    spm_image('init',func_files{1}{1});
    uiwait;
    
    %Reorient the rest of the func files (only from run1)
    if(exist(out_transform_file,'file'))
        % reorient_files=dir(fullfile(fileparts(func_files{1}{1}),'*_reorient.mat'));
        % for i=1:length(reorient_files),
        %     filedate(i,:)=datestr(reorient_files(i).date,'yyyymmddHHMMSS');
        % end;
        % [jnk,ii]=sortrows(filedate);
        % reorient_file=fullfile(fileparts(func_files{1}{1}),reorient_files(ii(end)).name);
        p=load(out_transform_file,'-ascii');
        
        % Ann Choe ----------------------------------------------------------
        %   loop through all runs, if there are multiple runs
        %
        
        for curr_run = 1:length(func_files)
            P=strvcat(func_files{curr_run});
            P(1,:)=[];
            reorient_nifti(P,p);
        end
        %---------------------------------------------------------------------
        
    end;
    clear global out_transform_file
    
    
end;

%% Open a sample image and try to get tr
vs = spm_vol(func_files{1}{1});
if(~exist('tr','var'))
    tr = vs.private.timing.tspace;
    if(tr > 50), %units are probably in msec
        tr = tr/1000;
    elseif(tr < 0.5)
        tr = tr.*1000;
    end;
    disp(['Verify TR is ',num2str(tr),' seconds']);
end;

%% Slice time correction
if(slice_time_correct)
    if mb<=1  % regular, not multiband
        load('rest_slice_time_job_template_spm8.mat'); % Load template
        matlabbatch{1}.spm.temporal.st.scans{1} = [];
        for ir=1:length(func_files)
            for i=1:length(func_files{ir}),
                [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
                matlabbatch{1}.spm.temporal.st.scans{ir}{i} = fullfile(func_dir,[prefix,func_file,ext]);
            end;
        end;
        matlabbatch{1}.spm.temporal.st.nslices=vs.dim(3); %Assumes that the 3rd dimension is slices
        matlabbatch{1}.spm.temporal.st.tr=tr;
        matlabbatch{1}.spm.temporal.st.ta=matlabbatch{1}.spm.temporal.st.tr - ...
            (matlabbatch{1}.spm.temporal.st.tr/matlabbatch{1}.spm.temporal.st.nslices);
        switch slice_order
            case 1
                matlabbatch{1}.spm.temporal.st.so = 1:matlabbatch{1}.spm.temporal.st.nslices;
            case 2
                matlabbatch{1}.spm.temporal.st.so = matlabbatch{1}.spm.temporal.st.nslices:-1:1;
            case 3
                matlabbatch{1}.spm.temporal.st.so = [1:2:matlabbatch{1}.spm.temporal.st.nslices, ...
                    2:2:matlabbatch{1}.spm.temporal.st.nslices];
            case 4
                matlabbatch{1}.spm.temporal.st.so = [2:2:matlabbatch{1}.spm.temporal.st.nslices, ...
                    1:2:matlabbatch{1}.spm.temporal.st.nslices];
            case 5
                matlabbatch{1}.spm.temporal.st.so = [matlabbatch{1}.spm.temporal.st.nslices:-2:1, ...
                    matlabbatch{1}.spm.temporal.st.nslices-1:-2:1];
            case 6
                matlabbatch{1}.spm.temporal.st.so = [matlabbatch{1}.spm.temporal.st.nslices-1:-2:1, ...
                    matlabbatch{1}.spm.temporal.st.nslices:-2:1];
            case 7
                matlabbatch{1}.spm.temporal.st.so = load(slice_order_file,'-ascii');
            case 8
                matlabbatch{1}.spm.temporal.st.so = [...
                    1:5:matlabbatch{1}.spm.temporal.st.nslices, ...
                    2:5:matlabbatch{1}.spm.temporal.st.nslices, ...
                    3:5:matlabbatch{1}.spm.temporal.st.nslices, ...
                    4:5:matlabbatch{1}.spm.temporal.st.nslices, ...
                    5:5:matlabbatch{1}.spm.temporal.st.nslices];
            otherwise
                error('Slice order has to be 1 - 7')
        end;
        matlabbatch{1}.spm.temporal.st.refslice = round(matlabbatch{1}.spm.temporal.st.nslices/2);  %??this depends on slice order too, need to correct
        save(fullfile(fileparts(func_files{1}{1}),'01_slice_time_correct_job.mat'),'matlabbatch');
        spm_jobman('run',matlabbatch);
        prefix=['a',prefix];
        clear matlabbatch;
    elseif mb>1  % multiband jhua
        P = [];
        for ir=1:length(func_files) % Each session
            P1 = [];
            for i=1:length(func_files{ir}), %Each dynamic
                [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
                P1{i} = fullfile(func_dir,[prefix,func_file,ext]);
%                 P{ir}{i} = fullfile(func_dir,[prefix,func_file,ext]);
            end;
            P{ir} = strvcat(P1);
        end;        
        
        % for historical reason, ignore multiband for ta calculation, see email
        V = spm_vol(P{1}(1,:));
        nslices = V(1).dim(3);
        ta = tr-tr/nslices;  %Acquisition Time (TA) {secs}
        timing(1) = ta / (nslices -1);  %time between last slices and next volume
        timing(2) = tr - ta;  %time between slices

        % slice timing, equidistant
        switch slice_order
            case 3  % Philips default
                timeperslice = 1000*tr/(nslices/mb); %ms
                oddslices = 1:2:nslices/mb;
                sliceorder(oddslices) = 0:timeperslice:(length(oddslices)-1)*timeperslice;  % ms
                evenslices = 2:2:nslices/mb;
                sliceorder(evenslices) = length(oddslices)*timeperslice:timeperslice:...
                    (length(oddslices)+length(evenslices)-1)*timeperslice;  % ms
                sliceorder = repmat(sliceorder, [1 mb]);
                refslice = sliceorder(evenslices(1));
            otherwise
                error('Slice order has to be 3 for multiband, other order not implemented')
        end;

        spm_slice_timing_mb(P, sliceorder, refslice, timing, 'a');
        prefix=['a',prefix];
    end
end;

%% Realignment
if(motion_correct)
    load('rest_realign_job_template_spm8.mat'); % Load template
    
    matlabbatch{1}.spm.spatial.realign.estimate.data{1} = [];
    for ir=1:length(func_files) % Each session
        for i=1:length(func_files{ir}), %Each dynamic
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            matlabbatch{1}.spm.spatial.realign.estimate.data{ir}{i} = ...
                fullfile(func_dir,[prefix,func_file,ext]);
        end;
        [func_dir,func_file]=fileparts(func_files{ir}{1});
        rp_file{ir}=fullfile(func_dir,['rp_',prefix,func_file,'.txt']); %If nuisance remove is selected
    end;
    save(fullfile(fileparts(func_files{1}{1}),'02_realign_job.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
    
end;

%% Coregistration
if(coregister)
    load('rest_coregister_job_template_spm8.mat'); % Load template
    [func_dir,func_file,ext] = fileparts(func_files{1}{3});
    % calculate the mean file, jhua
    img_mean = 0;
    n = 0;
    for ir=1:length(func_files) % Each session
        for i=floor(length(func_files{ir})/4):length(func_files{ir}), %Each dynamic, excluding starting scans
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            nii = load_nii(fullfile(func_dir,[prefix,func_file,ext]));
            img_mean = img_mean + nii.img;
            n = n + 1;
        end;
    end;
    img_mean = img_mean/n;
    nii.img = img_mean;
    mean_file = fullfile(func_dir,[prefix,func_file,'_mean',ext]);  % saved in the last run
    save_nii(nii,mean_file);
    clear nii;
    
%     matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = fullfile(func_dir,[prefix,func_file,ext]);
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = mean_file;   %jhua
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = anat_file;
    matlabbatch{1}.spm.spatial.coreg.estimate.other{1} = '';
    save(fullfile(fileparts(func_files{1}{1}),'03_coregister_job.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
end;

%% Segmentation
if(segment_anat),
    load('rest_segment_job_template_spm8.mat'); % Load template
    matlabbatch{1}.spm.spatial.preproc.data{1}=anat_file;
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{1}=fullfile(fileparts(which('spm')),'tpm','grey.nii');
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{2}=fullfile(fileparts(which('spm')),'tpm','white.nii');
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{3}=fullfile(fileparts(which('spm')),'tpm','csf.nii');
    matlabbatch{1}.spm.spatial.preproc.output.GM=[1,1,1];
    matlabbatch{1}.spm.spatial.preproc.output.WM=[1,1,1];
    matlabbatch{1}.spm.spatial.preproc.output.CSF=[1,1,1];
    save(fullfile(fileparts(func_files{1}{1}),'04_segment_job.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    
    [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
    wm_file=fullfile(anat_dir,['wc2',anat_file_nodir,ext]);
    csf_file=fullfile(anat_dir,['wc3',anat_file_nodir,ext]);
    clear matlabbatch;
end;

%% Spatial normalization
if(normalize),
    load('rest_normalize_job_template_spm8.mat'); % Load template
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [];
    
    n=1;
    for ir=1:length(func_files) %Each session
        for i=1:length(func_files{ir}), %Each dynamic
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n} = ...
                fullfile(func_dir,[prefix,func_file,ext]);
            n=n+1;
        end;
    end;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{n} = ...
        anat_file;
    [anat_dir,anat_file_nodir]=fileparts(anat_file);
    matlabbatch{1}.spm.spatial.normalise.write.subj.matname{1}= ...
        fullfile(anat_dir,[anat_file_nodir,'_seg_sn.mat']);
    save(fullfile(fileparts(func_files{1}{1}),'05_normalize_job.mat'),'matlabbatch');
    spm_jobman('run',matlabbatch);
    prefix=['w',prefix];
    clear matlabbatch;
end;



%% Create brain_mask - This is possible only if we had segmented and is needed only if we run one of the following!! - FIX THIS
if(create_brain_mask_file==1)
    fmri_create_brain_mask(anat_file,brain_mask_file);
end;


%% Remove NaN
for ir=1:length(func_files) % Each session
    for i=1:length(func_files{ir}), %Each dynamic
        [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
        nii = load_untouch_nii(fullfile(func_dir,[prefix,func_file,ext]));
        nii.img(isnan(nii.img)) = 0;
        save_untouch_nii(nii, fullfile(func_dir,[prefix,func_file,ext]));
    end;
end;


%% The following needs to be done session by session

for ir=1:length(func_files) 
    prefix_loop=prefix;
    
    
    %% HP filter or detrend
    if(hp_filter)
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop,func_file,ext]);
        end;
        fmri_time_filt(curr_func_files,tr,hp,2,brain_mask_file,low_ram);
        prefix_loop=['f',prefix_loop];
        clear curr_func_files func_dir func_file;
    end;
    
    %% Smooth the data - This is in case we want to run ICA
    if(ica_smooth)
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop,func_file,ext]);
        end;
        if(ica_smooth_fwhm==-1), ica_smooth_fwhm=round(abs(vs.mat(1,1).*2)); end;
        fmri_smooth_spm8(curr_func_files,ica_smooth_fwhm);
        clear curr_func_files func_dir func_file;
        prefix_loop2=['s',prefix_loop]; % prefix loop only for bp of smoothed data for ICA
        %don't change prefix here - no further processing necessary on these files
    end;
    
    %% Time domain filtering - for ICA as well
    if(bp_filter),
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop2,func_file,ext]);
        end;
        fmri_time_filt(curr_func_files,tr,bp,3,brain_mask_file,low_ram);
        
        %don't change prefix here - no further processing necessary on
        %these files
    end;
    
    
    %% Remove nuisance
    if(nuisance_remove)
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop,func_file,ext]);
        end;
        
        
        [func_dir,func_file,ext] = fileparts(func_files{ir}{1});
        nuisance_file = [func_dir,'/', prefix_loop, func_file, '-',nuisance_file_postfix,'.mat'];
        %         nuisance_file=[func_dir,nuisance_file_postfix,'.mat'];
        
        % Extract nuisances
        if(~motion_correct) %rp_file variable not defined
            tmp=dir(fullfile(fileparts(curr_func_files{1}),'rp*.txt'));
            rp_file{ir}=fullfile(fileparts(curr_func_files{1}),tmp(1).name);
        end;
        if(nr_options.wm_signal>0 && ~exist('wm_file','var')),
            [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
            wm_file=fullfile(anat_dir,['wc2',anat_file_nodir,ext]);
        else
            %wm_file='';
        end;
        if(nr_options.csf_signal>0 && ~exist('csf_file','var')),
            [anat_dir,anat_file_nodir,ext]=fileparts(anat_file);
            csf_file=fullfile(anat_dir,['wc3',anat_file_nodir,ext]);
        else
            %csf_file='';
        end;
        fmri_extract_nuisance(curr_func_files,tr,rp_file{ir},wm_file,csf_file,brain_mask_file,nr_options,nuisance_file);
        %extract_nuisance_20110310(func_dir{ir}, ['wa',func_flnm_mask{ir}],anat_dir,anat_flnm_mask,rp_file{ir},tr,hp,brain_mask_file,nuisance_file,nr_options);
        % Regress out nuisances
        if(low_ram~=1),
            fmri_regress_nuisance(curr_func_files,brain_mask_file,nuisance_file);
            %spm_regress_nuisance_dcn_20110310(func_dir{ir},['wa',func_flnm_mask{ir}],brain_mask_file,nuisance_file);
        else
            fmri_regress_nuisance_1D(curr_func_files,brain_mask_file,nuisance_file);
            %spm_regress_nuisance_1D_dcn_20110303(func_dir{ir},['wa',func_flnm_mask{ir}],brain_mask_file,nuisance_file);
        end;
        prefix_loop=['n',prefix_loop];
        clear curr_func_files func_dir func_file nuisance_file;
    end;
    
    
    
    %% Smooth the residuals
    if(var_smooth)
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop,func_file,ext]);
        end;
        if(smooth_fwhm==-1), smooth_fwhm=round(abs(vs.mat(1,1).*2)); end;
        fmri_smooth_spm8(curr_func_files,smooth_fwhm);
        prefix_loop=['s',prefix_loop];
        clear curr_func_files func_dir func_file;
    end;
    
    %% Time domain filtering
    if(bp_filter),
        for i=1:length(func_files{ir}),
            [func_dir,func_file,ext]=fileparts(func_files{ir}{i});
            curr_func_files{i}=fullfile(func_dir,[prefix_loop,func_file,ext]);
        end;
        fmri_time_filt(curr_func_files,tr,bp,3,brain_mask_file,low_ram);
        prefix_loop=['f',prefix_loop];
        clear curr_func_files func_dir func_file
    end;
    
end;
