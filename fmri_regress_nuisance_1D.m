function fmri_regress_nuisance_1D(func_files,brain_mask_file)
%Ying Cheng, Aug 19, 2012 for use with RS_LDDMM study
%Function to regress out nuisance covariates
%Nuisance covariates should be in the data directory (P's directory) in a
%file named nuisance.mat
%Output files are prepended with 'n' in the P filenames

%Ying Cheng, Feb 16, 2012

%% Read filenames of interest in func_dir
P=strvcat(func_files);
%vol_files = dir(fullfile(func_dir,[func_filenm_mask,'*.nii']));
%if(isempty(vol_files)),
%    vol_files = dir(fullfile(func_dir,[func_filenm_mask,'*.img']));
%end;
%for i_time = 1:length(vol_files),
%    P(i_time,:) = fullfile(func_dir,[vol_files(i_time).name,',1']);%#ok
%end;
%clear vol_files func_filenm_mask i_time

%% Read the brain mask
brain_mask = spm_read_vols(spm_vol(brain_mask_file));

%% Setup output headers
V = spm_vol(P);
Vo = V;
%hw = waitbar(0,'Copying header files');
for i_time = 1:size(V,1),
    Vo(i_time) = V(i_time);
    [pathname, filename, ext] = fileparts(V(i_time).fname);
    Vo(i_time).fname = fullfile(pathname,['n',filename,ext]);
    Vo(i_time).private.dat.fname = Vo(i_time).fname;
    if(strcmp(ext,'.img')),
        system(['cp ',strrep(V(i_time).fname,'.img','.hdr'),' ', ...
            strrep(Vo(i_time).fname,'.img','.hdr')]);
    end;
%    waitbar(i_time/size(V,1),hw);
end;
%close(hw);

%% Read nuisance vector
load(fullfile(fileparts(P(1,:)),'nuisance.mat'));
nt = size(nui.tc,1);%#ok
for i=1:size(nui.tc,2), nui.tc(:,i)=nui.tc(:,i)./max(abs(nui.tc(:,i))); end;
ntc = [nui.tc,ones(nt,1), (0:1/nt:(1-1/nt))'];% nuisances, constant, linear
precal_inv = ntc*((ntc'*ntc)\ntc');

%% Regress covariates
%hw=waitbar(0,'Regressing out nuisance covariates');
Y=zeros(V(1).dim(1),V(1).dim(2),size(V,1));
if(strcmp(ext,'.nii')),
    for i_time=1:size(V,1),
        spm_write_vol(Vo(i_time),zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)));
    end;
end;
for i_slice = 1:V(1).dim(3),
    %Load one slice (all time points)
    for i_time = 1:size(V,1),
        Y(:,:,i_time) = spm_slice_vol(V(i_time),spm_matrix([0,0,i_slice]), ...
            V(i_time).dim(1:2),0);%-mn_Y(i_time);
    end;
    %Get mean image (along time) to add back the image
    mY_time = mean(Y,3);
    %regress and get residuals
    for i_x = 1:size(Y,1),
        for i_y = 1:size(Y,2),
            if(brain_mask(i_x,i_y,i_slice) && any(Y(i_x,i_y,:))),
                tmpY=squeeze(Y(i_x,i_y,:)); 
                Y(i_x,i_y,:) = tmpY - precal_inv*tmpY;
            else
                Y(i_x,i_y,:) = zeros(1,1,size(V,1));
            end;
        end;
    end;
    %Write residuals after adding mean image back
    for i_time = 1:size(V,1),
        Y(:,:,i_time) = Y(:,:,i_time) + mY_time;
        spm_write_plane(Vo(i_time),Y(:,:,i_time),i_slice);
    end;
    %waitbar(i_slice/V(1).dim(3),hw);
end;
%close(hw);
