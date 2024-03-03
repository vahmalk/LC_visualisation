%% LC visualisation pipeline for 7T multi-echo FLASH data 
% using three MT weighted images for the fusion: one MT-T1 and two MT-PD images
% V Malekian, FIL Physics, 28/02/2024

clc; clear;close all;
%% Add path for NIFTI library and SPM12  
addpath('/export/home/vmalekian/MI_Distortion/NIfTI_20140122/')
addpath(genpath('/export/home/vmalekian/Malab_lib/spm12/'))

%% Specify all three MT-based data directories in the following order:
%%1- MT-T1 magnitude data directory 
%%2- MT-PD magnitude data directory (first acquisition)
%%3- MT-PD magnitude data directory (second acquisition)
directory3='/export/home/vmalekian/LC/Phil_study/M700867/mfc_3dflash_v1q_LC_4chos_FA24_0010/';
directory5='/export/home/vmalekian/LC/Phil_study/M700867/mfc_3dflash_v1q_LC_4chos_FA8_0012/';
directory7='/export/home/vmalekian/LC/Phil_study/M700867/mfc_3dflash_v1q_LC_4chos_FA8_0014/';
%% Specify T1w MPM directory
directory_wb='/export/home/vmalekian/LC/Phil_study/M700867/t1w_mfc_3dflash_v1k_0017/';
directory_output = '/export/home/vmalekian/LC/Phil_study/M700867/';
%% Specify wavelet coefficient directory
directory_coef = '/export/home/vmalekian/LC/Phil_study/';
%%SPM12 TMP path
spm_tmp_path = '/export/home/vmalekian/Malab_lib/spm12/tpm/';

unix(['fslroi ' (directory3) '*1.nii ' (directory3) 'e1_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory3) '*2.nii ' (directory3) 'e2_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory3) '*3.nii ' (directory3) 'e3_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory3) '*4.nii ' (directory3) 'e4_roi1_mg 0 -1 0 -1 6 60'])

unix(['fslroi ' (directory5) '*1.nii ' (directory5) 'e1_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory5) '*2.nii ' (directory5) 'e2_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory5) '*3.nii ' (directory5) 'e3_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory5) '*4.nii ' (directory5) 'e4_roi1_mg 0 -1 0 -1 6 60'])

unix(['fslroi ' (directory7) '*1.nii ' (directory7) 'e1_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory7) '*2.nii ' (directory7) 'e2_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory7) '*3.nii ' (directory7) 'e3_roi1_mg 0 -1 0 -1 6 60'])
unix(['fslroi ' (directory7) '*4.nii ' (directory7) 'e4_roi1_mg 0 -1 0 -1 6 60'])

%% TE-based echo averaging

%%magn MT-T1
unix(['fslmaths ' directory3 'e1_roi1_mg -div 2.75 ' directory3 'e1_roi1_mg_w'])
unix(['fslmaths ' directory3 'e2_roi1_mg -div 5.13 ' directory3 'e2_roi1_mg_w'])
unix(['fslmaths ' directory3 'e3_roi1_mg -div 7.51 ' directory3 'e3_roi1_mg_w'])
unix(['fslmaths ' directory3 'e4_roi1_mg -div 9.89 ' directory3 'e4_roi1_mg_w'])
unix(['fslmaths ' directory3 'e1_roi1_mg_w -add ' directory3 'e2_roi1_mg_w -add '  directory3 'e3_roi1_mg_w -add ' directory3 'e4_roi1_mg_w ' directory3 'e_roi1_mg_w_sum'])
unix(['fslmaths ' directory3 'e_roi1_mg_w_sum -div 0.7928 ' directory3 'e_roi1_mg_w_mean'])

%%magn MT-PD
unix(['fslmaths ' directory5 'e1_roi1_mg -div 2.75 ' directory5 'e1_roi1_mg_w'])
unix(['fslmaths ' directory5 'e2_roi1_mg -div 5.13 ' directory5 'e2_roi1_mg_w'])
unix(['fslmaths ' directory5 'e3_roi1_mg -div 7.51 ' directory5 'e3_roi1_mg_w'])
unix(['fslmaths ' directory5 'e4_roi1_mg -div 9.89 ' directory5 'e4_roi1_mg_w'])
unix(['fslmaths ' directory5 'e1_roi1_mg_w -add ' directory5 'e2_roi1_mg_w -add '  directory5 'e3_roi1_mg_w -add ' directory5 'e4_roi1_mg_w ' directory5 'e_roi1_mg_w_sum'])
unix(['fslmaths ' directory5 'e_roi1_mg_w_sum -div 0.7928 ' directory5 'e_roi1_mg_w_mean'])

%%magn MT-PD
unix(['fslmaths ' directory7 'e1_roi1_mg -div 2.75 ' directory7 'e1_roi1_mg_w'])
unix(['fslmaths ' directory7 'e2_roi1_mg -div 5.13 ' directory7 'e2_roi1_mg_w'])
unix(['fslmaths ' directory7 'e3_roi1_mg -div 7.51 ' directory7 'e3_roi1_mg_w'])
unix(['fslmaths ' directory7 'e4_roi1_mg -div 9.89 ' directory7 'e4_roi1_mg_w'])
unix(['fslmaths ' directory7 'e1_roi1_mg_w -add ' directory7 'e2_roi1_mg_w -add '  directory7 'e3_roi1_mg_w -add ' directory7 'e4_roi1_mg_w ' directory7 'e_roi1_mg_w_sum'])
unix(['fslmaths ' directory7 'e_roi1_mg_w_sum -div 0.7928 ' directory7 'e_roi1_mg_w_mean'])

%% SPM12 Bias correction

unix(['gunzip ' directory3 'e_roi1_mg_w_mean.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ directory3 'e_roi1_mg_w_mean.nii,1']};
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

unix(['gunzip ' directory5 'e_roi1_mg_w_mean.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ directory5 'e_roi1_mg_w_mean.nii,1']};
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

unix(['gunzip ' directory7 'e_roi1_mg_w_mean.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ directory7 'e_roi1_mg_w_mean.nii,1']};
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

%% Averging MT-PD images and registering MT-PD & T1 images

%%registering second MT-PD image to the first one
unix(['flirt -in ' directory7 'me_roi1_mg_w_mean -ref ' directory5 'me_roi1_mg_w_mean -out ' directory7 'me_roi1_mg_w_mean_reg -omat ' directory7 'me_roi1_mg_w_mean_reg.mat -interp spline -dof 6'])
unix(['fslmaths ' directory7 'me_roi1_mg_w_mean_reg -thr 0 ' directory7 'me_roi1_mg_w_mean_reg'])
unix(['fslroi ' directory3 'me_roi1_mg_w_mean ' directory3 'me_roi1_mg_w_mean_regr 0 -1 0 -1 1 58'])
unix(['fslroi ' directory5 'me_roi1_mg_w_mean ' directory5 'me_roi1_mg_w_mean_regr 0 -1 0 -1 1 58'])
unix(['fslroi ' directory7 'me_roi1_mg_w_mean_reg ' directory7 'me_roi1_mg_w_mean_regr 0 -1 0 -1 1 58'])
unix(['fslmaths ' directory5 'me_roi1_mg_w_mean_regr -add ' directory7 'me_roi1_mg_w_mean_regr -div 2 ' directory5 'me_roi1_mg_w_mean_regr_mt'])

%%registering MT-T1 to the averaged MT-PD volume
unix(['flirt -in ' directory3 'me_roi1_mg_w_mean_regr -ref ' directory5 'me_roi1_mg_w_mean_regr_mt -out ' directory3 'me_roi1_mg_w_mean_regr2mt -omat ' directory3 'me_roi1_mg_w_mean_regr_t1_reg.mat -interp spline -dof 6'])
unix(['fslmaths ' directory3 'me_roi1_mg_w_mean_regr2mt -thr 0 ' directory3 'me_roi1_mg_w_mean_regr2mt'])

%% MT-T1 & PD Multiplication
unix(['fslmaths ' directory5 'me_roi1_mg_w_mean_regr_mt -mul ' directory3 'me_roi1_mg_w_mean_regr2mt ' directory5 'me_roi1_mg_w_mean_regr_mtt1_mul'])
unix(['fslmaths ' directory5 'me_roi1_mg_w_mean_regr_mtt1_mul -sqrt ' directory5 'me_roi1_mg_w_mean_regr_mtt1_sqrt'])

%% Load MT weighted data 
string_mt='me_roi1_mg_w_mean_regr_mt.nii.gz';
a0mt=dir([directory5,string_mt]);
nf_mt=load_untouch_nii([directory5,a0mt(1).name]);

string_t1='me_roi1_mg_w_mean_regr2mt.nii.gz';
a0t1=dir([directory3,string_t1]);
nf_t1=load_untouch_nii([directory3,a0t1(1).name]);


mt=nf_mt.img;
t1=nf_t1.img;

%% Wavelet decomposition 
wl = 'db1';
Wt1 = dwt3(t1, wl);
Ft1 = (Wt1.dec);
Wmt = dwt3(mt, wl);
Fmt = (Wmt.dec);

%% Load wavelet coefficient, obtained from one train dataset

load([directory_coef 'coef_mt_pd.mat'])
load([directory_coef 'coef_mt_t1.mat'])
load([directory_coef 'b_mt_pd.mat'])
load([directory_coef 'b_mt_t1.mat'])

coef_mt_d = coef_mt(2:8);
coef_t1_d = coef_t1(2:8);

outputMin = 1; 
outputMax = 2; 

mt_inputMin = min(coef_mt_d)
mt_inputMax = max(coef_mt_d); 
t1_inputMin = min(coef_t1_d); 
t1_inputMax = max(coef_t1_d); 
mt_slop_avg_norm_mapped = ((coef_mt_d - mt_inputMin) / (mt_inputMax - mt_inputMin)) * (outputMax - outputMin) + outputMin;
t1_slop_avg_norm_mapped = ((coef_t1_d - t1_inputMin) / (t1_inputMax - t1_inputMin)) * (outputMax - outputMin) + outputMin;

%% Scale MT-T1 seperately

Ctf{1,1,1}  = 1*Ft1{1,1,1};
Ctf{1,1,2}  = t1_slop_avg_norm_mapped(1)*Ft1{1,1,2};
Ctf{1,2,1}  = t1_slop_avg_norm_mapped(2)*Ft1{1,2,1};
Ctf{1,2,2}  = t1_slop_avg_norm_mapped(3)*Ft1{1,2,2};
Ctf{2,1,1}  = t1_slop_avg_norm_mapped(4)*Ft1{2,1,1};
Ctf{2,1,2}  = t1_slop_avg_norm_mapped(5)*Ft1{2,1,2};
Ctf{2,2,1}  = t1_slop_avg_norm_mapped(6)*Ft1{2,2,1};
Ctf{2,2,2}  = t1_slop_avg_norm_mapped(7)*Ft1{2,2,2};

Wt1.dec = Ctf;
t1_rec_final = idwt3(Wt1);
string_magn='t1_new2.nii.gz';
save_nii(make_nii((t1_rec_final)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 't1_new2.nii.gz'])
unix(['fslmaths ' directory3 't1_new2 -thr 0 ' directory3 't1_new2'])

hg = fspecial3('gaussian',[3 3 5],0.4);
volSmooth_t1 = imfilter(t1_rec_final,hg);
string_magn='t1_new2_sm.nii.gz';
save_nii(make_nii((volSmooth_t1)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 't1_new2_sm.nii.gz'])

%% Scale MT-PD seperately

Cmt{1,1,1}  = 1*Fmt{1,1,1};
Cmt{1,1,2}  = mt_slop_avg_norm_mapped(1)*Fmt{1,1,2};
Cmt{1,2,1}  = mt_slop_avg_norm_mapped(2)*Fmt{1,2,1};
Cmt{1,2,2}  = mt_slop_avg_norm_mapped(3)*Fmt{1,2,2};
Cmt{2,1,1}  = mt_slop_avg_norm_mapped(4)*Fmt{2,1,1};
Cmt{2,1,2}  = mt_slop_avg_norm_mapped(5)*Fmt{2,1,2};
Cmt{2,2,1}  = mt_slop_avg_norm_mapped(6)*Fmt{2,2,1};
Cmt{2,2,2}  = mt_slop_avg_norm_mapped(7)*Fmt{2,2,2};
Wmt.dec = Cmt;
mt_rec_final = idwt3(Wmt);
string_magn='mt_new2.nii.gz';
save_nii(make_nii((mt_rec_final)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 'mt_new2.nii.gz'])
unix(['fslmaths ' directory3 'mt_new2 -thr 0 ' directory3 'mt_new2'])

volSmooth_mt1 = imfilter(mt_rec_final,hg);
string_magn='mt_new2_sm.nii.gz';
save_nii(make_nii((volSmooth_mt1)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 'mt_new2_sm.nii.gz'])


%% MT-PD & MT-T1 fusion

inputMin = min([coef_t1_d coef_mt_d]); 
inputMax = max([coef_t1_d coef_mt_d]);

mt_slop_avg_norm_mapped_comb = ((coef_mt_d - inputMin) / (inputMax - inputMin)) * (outputMax - outputMin) + outputMin;
t1_slop_avg_norm_mapped_comb = ((coef_t1_d - inputMin) / (inputMax - inputMin)) * (outputMax - outputMin) + outputMin;

% x=2:8;
% figure,plot(x,mt_slop_avg_norm_mapped_comb,'r*'),hold on,plot(x,t1_slop_avg_norm_mapped_comb,'bo'),title('Sub-band Weights'),xlabel('SB'),ylabel('Weights'),legend('mt','t1')%,ylim([-20 100])

Ctf{1,1,1}  = 2*(b_t1_sub1/(b_mt_sub1+b_t1_sub1))*Ft1{1,1,1}+2*(b_mt_sub1/(b_mt_sub1+b_t1_sub1))*Fmt{1,1,1};
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


Wt1.dec = Ctf;
mtt1_rec_final = idwt3(Wt1);
string_magn='mtt1_new2.nii.gz';
save_nii(make_nii((mtt1_rec_final)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 'mtt1_new2.nii.gz'])
unix(['fslmaths ' directory3 'mtt1_new2 -thr 0 ' directory3 'mtt1_new2'])

hg = fspecial3('gaussian',[3 3 5],0.4);
volSmooth_t1 = imfilter(mtt1_rec_final,hg);
string_magn='mtt1_new2_sm.nii.gz';
save_nii(make_nii((volSmooth_t1)),[directory3 string_magn])
unix(['fslcpgeom ' directory5 'me_roi1_mg_w_mean_regr_mt.nii.gz ' directory3 'mtt1_new2_sm.nii.gz'])
unix(['fslmaths ' directory3 'mtt1_new2_sm -thr 0 ' directory_output 'LC_wavelet'])

%% Register to whole-brain MPM T1w data

unix(['fslmaths ' directory_wb 's*1.nii -div 2.75 ' directory_wb 'e1_roi1_mg_w'])
unix(['fslmaths ' directory_wb 's*2.nii -div 5.13 ' directory_wb 'e2_roi1_mg_w'])
unix(['fslmaths ' directory_wb 's*3.nii -div 7.51 ' directory_wb 'e3_roi1_mg_w'])
unix(['fslmaths ' directory_wb 's*4.nii -div 9.89 ' directory_wb 'e4_roi1_mg_w'])

unix(['fslmaths ' directory_wb 'e1_roi1_mg_w -add ' directory_wb 'e2_roi1_mg_w -add '  directory_wb 'e3_roi1_mg_w -add ' directory_wb 'e4_roi1_mg_w ' directory_wb 'e_roi1_mg_w_sum'])
unix(['fslmaths ' directory_wb 'e_roi1_mg_w_sum -div 0.7928 ' directory_wb 'e_roi1_mg_w_mean'])

unix(['fslreorient2std ' directory_wb 'e_roi1_mg_w_mean ' directory_wb 'e_roi1_mg_w_mean_reo'])

unix(['gunzip ' directory_wb 'e_roi1_mg_w_mean_reo.nii.gz'])
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[ directory_wb 'e_roi1_mg_w_mean_reo.nii,1']};
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

unix(['fslmaths ' directory_wb 'me_roi1_mg_w_mean_reo -thr 0 ' directory_wb 'WB'])
unix(['flirt -in ' directory3 'me_roi1_mg_w_mean_regr2mt -ref ' directory_wb 'me_roi1_mg_w_mean_reo -out ' directory_wb 'me_roi1_mg_w_mean_regr2mt2wb -omat ' directory_wb 'mraw_mean_reo_regis.mat -interp spline -dof 12'])
unix(['flirt -in ' directory3 'mtt1_new2_sm -ref ' directory_wb 'me_roi1_mg_w_mean_reo.nii  -out ' directory_wb 'mtt1_new2_sm2wb -applyxfm -init ' directory_wb 'mraw_mean_reo_regis.mat -interp spline']) 
unix(['fslmaths ' directory_wb 'mtt1_new2_sm2wb -thr 0 ' directory_wb 'LC_wavelet2WB'])



