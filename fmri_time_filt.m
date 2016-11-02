function fmri_time_filt(func_files,tr,cutoffs,band,brain_mask_filename,low_ram)
%Function to filter time series in fMRI series of image
%Usage
%   fmri_time_filt(func_files,tr,cutoffs,band,brain_mask_filename,low_ram)
%       func_files is a array of cell strings with fMRI filenames 
%       tr - is the repetition time in seconds
%       cutoffs - is the cut-off frequencies for the filter in Hertz
%       band - type of filter 1=low-pass, 2=high-pass, 3=band-pass,
%           4=stop-band;5=linear-detrend (cutoffs is ignored for band=5)
%       brain_mask_filename - is the filename of the brain mask
%       low_ram - is set to 1 if the computer has low RAM (< 4GB)
%Output files are prepended with 'f' in the filenames of the input files

%July 10, 2012, Aug 22 - Added detrend capability and optimized 4D analysis
%July 10, 2012, Feb 5, 2012


%% Initialize
%Read filenames of interest in func_dir
if(iscell(func_files)), P=strvcat(func_files); else P=func_files; end;
disp(['Filtering ',fileparts(P(1,:))]);


%% Read data
V=spm_vol(P);
Vo=V;
switch band
case 1 % low pass
    [b,a]=butter(2,2*tr*cutoffs);
case 2 % high pass
    [b,a]=butter(2,2*tr*cutoffs,'high');
case 3 % band pass
    [b,a]=butter(2,[2*tr*cutoffs(1),2*tr*cutoffs(2)]);
case 4 % stop band
    [b,a]=butter(2,[2*tr*cutoffs(1),2*tr*cutoffs(2)],'stop');
case 5 % linear detrend
    b=[];a=[];
end;

brain_mask = logical(spm_read_vols(spm_vol(brain_mask_filename)));


%% Copy Header files
for i_time=1:size(V,1),
    [pathname,filename,ext]=fileparts(V(i_time).fname);
    Vo(i_time).fname=fullfile(pathname,['f',filename,ext]);
    Vo(i_time).private.dat.fname=Vo(i_time).fname;
    if(strcmp(ext,'.img')),
        copyfile(strrep(V(i_time).fname,'.img','.hdr'), ...
            strrep(Vo(i_time).fname,'.img','.hdr'));
    end;
end;

%% Filter time series and write files simultaneously
if(low_ram~=1) %When computer sufficient RAM (usually 4GB)
    Y=spm_read_vols(V);
    Y=reshape(Y,prod(V(1).dim),length(V));
    if(band<=4),
        Y(brain_mask,:)= filtfilt(b,a,Y(brain_mask,:)')' + ...
            repmat(mean(Y(brain_mask,:),2),1,length(V));
    else
        Y(brain_mask,:)=detrend(Y(brain_mask,:)')' + ...
            repmat(mean(Y(brain_mask,:),2),1,length(V));
    end;
    Y=reshape(Y,V(1).dim(1),V(1).dim(2),V(1).dim(3),length(V));
    for i_time=1:size(V,1),
        spm_write_vol(Vo(i_time),Y(:,:,:,i_time));
    end;
else %low RAM computers (for computers with < 4GB RAM - usually never used nowadays)
    if(strcmp(ext,'.nii')),
        for i_time=1:size(V,1),
            spm_write_vol(Vo(i_time),zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)));
        end;
    end;
    Y=zeros(V(1).dim(1),V(1).dim(2),size(V,1));
    for i_slice=1:V(1).dim(3),%Assume 1st Vol has same 3rd dim as the rest 
        %Load one slice (all time points)
        for i_time=1:size(V,1),
            Y(:,:,i_time)=spm_slice_vol(V(i_time),spm_matrix([0,0,i_slice]), ...
                V(i_time).dim(1:2),0);
        end;
        %Filtering
        for i_x=1:size(Y,1),
            for i_y=1:size(Y,2),
                if(brain_mask(i_x,i_y,i_slice))
                    Y(i_x,i_y,:)=filtfilt(b,a,Y(i_x,i_y,:).*100)./100 + ...
                        mean(Y(i_x,i_y,:));
                end;
            end;
        end;
        for i_time=1:size(V,1),
            spm_write_plane(Vo(i_time),Y(:,:,i_time),i_slice);
        end;
    end;
end;
