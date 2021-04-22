%%

clear all
close all
clc

% Read data and standardize to a given template

addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath(genpath('C:\MatlabLibs\NIFTIread'));

addpath(pwd);
cDir = pwd;
dDir = 'U:\Data_Organized\UTurku\PRO3\';
% dDir = 'U:\Data_Organized\UTurku\SUPP\';
[studies,pID] = xlsread('U:\docs\parameters_PRO3');
% [studies,pID] = xlsread('U:\docs\parameters_SUPP');
pID = pID(2:end,2);

label = studies(:,3);

%%
% Read template for standardization
cd('X:\Annotations_Purysko\1');

T2 = load_untouch_nii('Pat1_T2_N4_Or.nii.gz'); template_T2 = T2.img;
ADC = load_untouch_nii('Pat1_T2_ADC_Aff.nii.gz'); template_ADC = ADC.img;
cgMask = load_untouch_nii('Pat1_T2_N4_CG_Or.nii.gz'); template_cgMask = repeatFuncOver3D(@im2bw,double(cgMask.img));
caMask = load_untouch_nii('Pat1_T2_N4_Or_Cancer.nii.gz'); template_caMask = repeatFuncOver3D(@im2bw,double(caMask.img));
mask = load_untouch_nii('Pat1_T2_N4_Prostate_Or.nii.gz'); template_mask = repeatFuncOver3D(@im2bw,double(mask.img));

templateVolMasked_T2 = double(template_T2).*template_mask;
templateVolMasked_ADC = double(template_ADC).*template_mask;
opts.temcancermasks = logical(template_caMask);
opts.docheck = false;
opts.dorescale = false;
% opts.numstdtopclip = 0;
%%


for i = 1:length(studies)
    disp(['reading study - ' cell2mat(pID(i))]);
    cd(dDir);
    
    T2 = load_untouch_nii([cell2mat(pID(i)) '_T2_Or.nii.gz']);T2 = (T2.img);
    ADC = load_untouch_nii([cell2mat(pID(i))  '_ADC_Rig.nii.gz']); ADC = double(ADC.img);
    mask = load_untouch_nii([cell2mat(pID(i)) '_Prostate_Or.nii.gz']);mask = repeatFuncOver3D(@im2bw,double(abs(mask.img)));
    caMask = load_untouch_nii([cell2mat(pID(i))  '_Lesion_Or.nii.gz']);caMask = repeatFuncOver3D(@im2bw,double(abs(caMask.img)));
    
    opts.incancermasks = logical(caMask);
    
    inputVolMask = double(T2).*mask;
    [~,stdMap,~] = int_stdn_landmarks(inputVolMask,templateVolMasked_T2,opts);
    close all;
    
    maskDil = dilateMaskVol(mask,30);
    inputVolMaskDil = double(T2).*maskDil;
    
    T2std = applystdnmap_rs(inputVolMaskDil,stdMap);
    
    ADC(ADC<0) = 0;
    inputVolMask = double(ADC).*mask;
    [~,stdMap,scalefactends] = int_stdn_landmarks(inputVolMask,templateVolMasked_ADC,opts);
    close all;
    
    maskDil = dilateMaskVol(mask,30);
    inputVolMaskDil = double(ADC).*maskDil;

    ADCstd = applystdnmap_rs(inputVolMaskDil,stdMap);
    
    indices = []; feats = [];
    for j = 1:size(mask,3)
        mask_ = mask(:,:,j);
        if max(mask_(:)>0);
%             T2feats_ =  computeTextureFeatures(T2std(:,:,j),mask_);
%             ADCfeats_ = computeTextureFeatures(ADC(:,:,j),mask_);

            T2feats_ =  computeCollageFeatures(T2std(:,:,j),mask_);
            ADCfeats_ = computeCollageFeatures(ADC(:,:,j),mask_);
            
            feats = [feats; T2feats_ ADCfeats_];
%             feats = T2feats_;
%             feats1 = [feats1; ADCfeats_];
            
            caMask_ = caMask(:,:,j);
%             cgMask_ = cgMask(:,:,j);
            
            temp = find(mask_==1);
            indices = [ indices; j*ones(size(temp)) temp caMask_(temp) label(i)*ones(size(temp))];  %  slice# pixel# cancerPresence#  BCRpresence#
        end
    end
    
%     load(['C:\Code\CCF_BCR\features\Pat' num2str(pID) '_featsPixel.mat']);
%     feats(:,64:end) = feats1;
    
    cd('C:\Code\CCF_BCR\features_UTurku\features_collage_PRO3');
%     cd('C:\Code\CCF_BCR\features_UTurku\features_SUPP');
    save(['Pat' num2str(i) '_featsPixel.mat'],'feats','indices');
    clear feats indices T2 ADC mask cgMask caMask inputVolMask T2std inputVolMaskDil
end