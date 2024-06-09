%% LC visualisation pipeline for 7T multi-echo FLASH data 
% using two MT weighted images for the fusion: one MT-T1 and one MT-PD images
% V Malekian, FIL Physics, 09/06/2024

clc;clear;
%% Add path for NIFTI library and SPM12  
addpath('/export/home/vmalekian/Malab_lib/NIfTI_20140122/')
addpath(genpath('/export/home/vmalekian/Matlab_lib/spm12/'))

%% Specify all three MT-based data directories in the following order:
%%1- MT-T1 magnitude data directory 
%%2- MT-PD magnitude data directory 
dir_T1 =  '/export/home/vmalekian/LC_fMRI_EPI/misun_LC/code_test/mfc_3dflash_v1q_LC_4chos_FA24_0029/';
dir_PD =  '/export/home/vmalekian/LC_fMRI_EPI/misun_LC/code_test/mfc_3dflash_v1q_LC_4chos_FA8_0031/';

%% Load Wavelet Coefficients
load('/export/home/vmalekian/LC_fMRI_EPI/misun_LC/code_test/t1_slop_avg_norm.mat')  
load('/export/home/vmalekian/LC_fMRI_EPI/misun_LC/code_test/mt_slop_avg_norm.mat')  
load('/export/home/vmalekian/LC_fMRI_EPI/misun_LC/code_test/b_coarse.mat')  

%% Remove top and bottom slices 
unix(['fslroi ' (dir_T1) '*1.nii ' (dir_T1) 'e1_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_T1) '*2.nii ' (dir_T1) 'e2_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_T1) '*3.nii ' (dir_T1) 'e3_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_T1) '*4.nii ' (dir_T1) 'e4_roi1_mg 0 -1 0 -1 6 60'])

unix(['fslroi ' (dir_PD) '*1.nii ' (dir_PD) 'e1_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_PD) '*2.nii ' (dir_PD) 'e2_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_PD) '*3.nii ' (dir_PD) 'e3_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (dir_PD) '*4.nii ' (dir_PD) 'e4_roi1_mg 0 -1 0 -1 6 60'])


%% TE-based echo averaging
%%magn MT-T1
unix(['fslmaths ' dir_T1 'e1_roi1_mg -div 2.75 ' dir_T1 'e1_roi1_mg_w'])
unix(['fslmaths ' dir_T1 'e2_roi1_mg -div 5.13 ' dir_T1 'e2_roi1_mg_w'])
unix(['fslmaths ' dir_T1 'e3_roi1_mg -div 7.51 ' dir_T1 'e3_roi1_mg_w'])
unix(['fslmaths ' dir_T1 'e4_roi1_mg -div 9.89 ' dir_T1 'e4_roi1_mg_w'])
unix(['fslmaths ' dir_T1 'e1_roi1_mg_w -add ' dir_T1 'e2_roi1_mg_w -add '  dir_T1 'e3_roi1_mg_w -add ' dir_T1 'e4_roi1_mg_w ' dir_T1 'e_roi1_mg_w_sum'])
unix(['fslmaths ' dir_T1 'e_roi1_mg_w_sum -div 0.7928 ' dir_T1 'e_roi1_mg_w_mean'])

%%magn MT-PD
unix(['fslmaths ' dir_PD 'e1_roi1_mg -div 2.75 ' dir_PD 'e1_roi1_mg_w'])
unix(['fslmaths ' dir_PD 'e2_roi1_mg -div 5.13 ' dir_PD 'e2_roi1_mg_w'])
unix(['fslmaths ' dir_PD 'e3_roi1_mg -div 7.51 ' dir_PD 'e3_roi1_mg_w'])
unix(['fslmaths ' dir_PD 'e4_roi1_mg -div 9.89 ' dir_PD 'e4_roi1_mg_w'])
unix(['fslmaths ' dir_PD 'e1_roi1_mg_w -add ' dir_PD 'e2_roi1_mg_w -add '  dir_PD 'e3_roi1_mg_w -add ' dir_PD 'e4_roi1_mg_w ' dir_PD 'e_roi1_mg_w_sum'])
unix(['fslmaths ' dir_PD 'e_roi1_mg_w_sum -div 0.7928 ' dir_PD 'e_roi1_mg_w_mean'])


%% Bias correction
%%SPM12 TMP path
spm_tmp_path = '/export/home/vmalekian/Malab_lib/spm12/tpm/';
unix(['gunzip ' dir_T1 'e_roi1_mg_w_mean.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ dir_T1 'e_roi1_mg_w_mean.nii,1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_tmp_path 'TPM.nii,1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_tmp_path 'TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_tmp_path 'TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spm_tmp_path 'TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spm_tmp_path 'TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spm_tmp_path 'TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm('defaults','FMRI')
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

unix(['gunzip ' dir_PD 'e_roi1_mg_w_mean.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ dir_PD 'e_roi1_mg_w_mean.nii,1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_tmp_path 'TPM.nii,1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_tmp_path 'TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_tmp_path 'TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spm_tmp_path 'TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spm_tmp_path 'TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spm_tmp_path 'TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm('defaults','FMRI')
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% Brain masking
unix(['bet ' dir_T1 'me_roi1_mg_w_mean.nii.gz ' dir_T1 'me_roi1_mg_w_mean_brain -Z -f 0.5 -g 0  -m'])
unix(['fslcpgeom ' dir_T1 'me_roi1_mg_w_mean.nii.gz ' dir_T1 'me_roi1_mg_w_mean_brain.nii.gz'])

%% Registering MT-PD & T1 images
unix(['flirt -in ' dir_PD 'me_roi1_mg_w_mean -ref ' dir_T1 'me_roi1_mg_w_mean -out ' dir_PD 'me_roi1_mg_w_mean_temp -omat ' dir_PD 'me_roi1_mg_w_mean_temp.mat -interp spline -dof 6'])
unix(['convert_xfm -inverse -omat ' dir_PD 'me_roi1_mg_w_mean_temp_inv.mat ' dir_PD 'me_roi1_mg_w_mean_temp.mat'])
unix(['flirt -ref ' dir_PD 'me_roi1_mg_w_mean -in ' dir_T1 'me_roi1_mg_w_mean_brain_mask -out ' dir_PD 'me_roi1_mg_w_mean_brain_mask2mt -applyxfm -init ' dir_PD 'me_roi1_mg_w_mean_temp_inv.mat -interp nearestneighbour'])
unix(['fslmaths ' dir_PD 'me_roi1_mg_w_mean -mul ' dir_PD 'me_roi1_mg_w_mean_brain_mask2mt ' dir_PD 'me_roi1_mg_w_mean_brain'])
unix(['flirt -in ' dir_PD 'me_roi1_mg_w_mean_brain -ref ' dir_T1 'me_roi1_mg_w_mean_brain -out ' dir_PD 'me_roi1_mg_w_mean_brain2t1 -omat ' dir_PD 'me_roi1_mg_w_mean_brain2t1.mat -interp spline -dof 6'])
unix(['fslmaths ' dir_PD 'me_roi1_mg_w_mean_brain2t1 -thr 0 ' dir_PD 'me_roi1_mg_w_mean_brain2t1'])


%% T1 PD Multiplication
unix(['fslmaths ' dir_PD 'me_roi1_mg_w_mean_brain2t1 -mul ' dir_T1 'me_roi1_mg_w_mean_brain ' dir_PD 'me_roi1_mg_w_mean_masked_mul'])
unix(['fslmaths ' dir_PD 'me_roi1_mg_w_mean_masked_mul -sqrt ' dir_PD 'me_roi1_mg_w_mean_masked_sqrt'])


%% Load MT-PD and MT-T1 data 
string_mt='me_roi1_mg_w_mean_brain2t1.nii.gz';
a0mt=dir([dir_PD,string_mt]);
nf_mt=load_untouch_nii([dir_PD,a0mt(1).name]);

string_t1='me_roi1_mg_w_mean_brain.nii.gz';
a0t1=dir([dir_T1,string_t1]);
nf_t1=load_untouch_nii([dir_T1,a0t1(1).name]);

mt=nf_mt.img;
t1=nf_t1.img;

%% Apply wavelet coefficients
wl = 'db1';
Wt1 = dwt3(t1, wl);
Ft1 = (Wt1.dec);

Wmt = dwt3(mt, wl);
Fmt = (Wmt.dec);

mt_slop_avg_norm_mapped_comb = mt_slop_avg_norm_mapped_comb_total_mean;
t1_slop_avg_norm_mapped_comb = t1_slop_avg_norm_mapped_comb_total_mean;

b_t1_sub1 = b_coarse_mean(1);
b_mt_sub1 = b_coarse_mean(2);

Ctf{1,1,1}  = (2*b_t1_sub1/(b_mt_sub1+b_t1_sub1))*Ft1{1,1,1}+(2*b_mt_sub1/(b_mt_sub1+b_t1_sub1))*Fmt{1,1,1};
Ctf{1,1,2}  = t1_slop_avg_norm_mapped_comb(1)*Ft1{1,1,2}+mt_slop_avg_norm_mapped_comb(1)*Fmt{1,1,2};
Ctf{1,2,1}  = t1_slop_avg_norm_mapped_comb(2)*Ft1{1,2,1}+mt_slop_avg_norm_mapped_comb(2)*Fmt{1,2,1};
Ctf{1,2,2}  = t1_slop_avg_norm_mapped_comb(3)*Ft1{1,2,2}+mt_slop_avg_norm_mapped_comb(3)*Fmt{1,2,2};
Ctf{2,1,1}  = t1_slop_avg_norm_mapped_comb(4)*Ft1{2,1,1}+mt_slop_avg_norm_mapped_comb(4)*Fmt{2,1,1};
Ctf{2,1,2}  = t1_slop_avg_norm_mapped_comb(5)*Ft1{2,1,2}+mt_slop_avg_norm_mapped_comb(5)*Fmt{2,1,2};
Ctf{2,2,1}  = t1_slop_avg_norm_mapped_comb(6)*Ft1{2,2,1}+mt_slop_avg_norm_mapped_comb(6)*Fmt{2,2,1};
Ctf{2,2,2}  = t1_slop_avg_norm_mapped_comb(7)*Ft1{2,2,2}+mt_slop_avg_norm_mapped_comb(7)*Fmt{2,2,2};

Ctf{1,1,1}= 0.5*Ctf{1,1,1};Ctf{1,1,2}= 0.5*Ctf{1,1,2};
Ctf{1,2,1}= 0.5*Ctf{1,2,1};Ctf{1,2,2}= 0.5*Ctf{1,2,2};
Ctf{2,1,1}= 0.5*Ctf{2,1,1};Ctf{2,1,2}= 0.5*Ctf{2,1,2};
Ctf{2,2,1}= 0.5*Ctf{2,2,1};Ctf{2,2,2}= 0.5*Ctf{2,2,2};

%% Reconsruction
Wt1.dec = Ctf;
mtt1_rec_final = idwt3(Wt1);
string_magn='WaveletContarst.nii.gz';
save_nii(make_nii((mtt1_rec_final)),[dir_T1 string_magn])
unix(['fslcpgeom ' dir_T1 'me_roi1_mg_w_mean_brain.nii.gz ' dir_T1 'WaveletContarst.nii.gz'])
unix(['fslmaths ' dir_T1 'WaveletContarst -thr 0 ' dir_T1 'WaveletContarst'])

