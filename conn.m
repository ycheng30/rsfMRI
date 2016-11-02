% To calculate correlations to the seeds
% Created by Ying Cheng
% multi runs, Ying Cheng 09/03/2014
% add tSNR map, Ying Cheng 09/03/2014

function conn(data_path,rundir,anatdir,data_filenm_mask,roi_filenm_mask,seedbased_output_dir,brain_mask_filenm, nruns)

if(nargin<8) || isempty(nruns),
    nruns = 1;
end;

% load brain mask
% vdata = spm_vol('H:\data\RS\Children\F_grpmask.img');
vdata = spm_vol(fullfile(data_path,anatdir,[brain_mask_filenm,'.nii']));
v0 = spm_read_vols(vdata);

% load images, convert to time series
tseries = [];
for ir = 1:nruns
    imgs = dir(fullfile(data_path,[rundir num2str(ir)],[data_filenm_mask,'*.nii']));
    imgs = {imgs.name}';

    tseries1 = zeros(sum(v0(:)>0), size(imgs,1));

    for iimg = 1:size(imgs,1)
        temp = spm_read_vols(spm_vol(fullfile(data_path,[rundir num2str(ir)],imgs{iimg})));
        temp = temp(:);
        tseries1(:, iimg) = temp(v0(:)>0);
    end
    
    tseries = cat(2, tseries, tseries1);
end

% temporal SNR (tSNR)
tsnr = mean(tseries,2) ./ std(tseries,[],2);

% load seed time course
roi1 = load(fullfile(data_path,seedbased_output_dir,[roi_filenm_mask{1},'_',data_filenm_mask,'_roi_tc.mat']));
roi = roi1.tc';
    
% calculate correlations, r to z transform
for ivox=1:size(tseries,1)
    cc=corrcoef(roi',tseries(ivox,:));
    Yo(ivox,1)=cc(1,2);
    Yo(ivox,1)=fisher_r2z(Yo(ivox,:));
end;        

% write the correlation map
result = v0;
result(v0~=0) = Yo;
v1 = vdata;
v1.fname = fullfile(data_path,seedbased_output_dir,['conn_',roi_filenm_mask{1},'_roi_tc.nii']);
v1.private.dat.fname=v1.fname;
v1.dt=[16 0];
spm_write_vol(v1, result);

% write the correlation map with NaN
result = v0;
result(v0~=0) = Yo;
result(v0==0) = NaN;
v1 = vdata;
v1.fname = fullfile(data_path,seedbased_output_dir,['conn_NaN_',roi_filenm_mask{1},'_roi_tc.nii']);
v1.private.dat.fname=v1.fname;
v1.dt=[16 0];
spm_write_vol(v1, result);

% write the tSNR map
result = v0;
result(v0~=0) = tsnr;
v1 = vdata;
v1.fname = fullfile(data_path,seedbased_output_dir,['tsnr_',roi_filenm_mask{1},'_roi_tc.nii']);
v1.private.dat.fname=v1.fname;
v1.dt=[16 0];
spm_write_vol(v1, result);


    

   


