clear all
close all

addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath(genpath('C:\MatlabLibs\NIFTIread'));

dDir = 'X:\Annotations_Purysko\';
nP = 82;
T2 = load_untouch_nii([dDir num2str(nP) '\Pat' num2str(nP) '_T2_N4_Or.nii.gz']); T2 = T2.img;
% ADC = load_untouch_nii([dDir num2str(nP) '\Pat' num2str(nP) '_T2_ADC_Aff.nii.gz']); ADC = ADC.img;

load(['.\features_all\Pat' num2str(nP) '_featsPixel.mat']);

slices = unique(indices(:,1));

for i = 1:length(slices)
% i=6;
    tempInd = indices(indices(:,1)==slices(i),:);
    if max(tempInd(:,3))>0
        T2slice = T2(:,:,slices(i)); tempFeats = feats(indices(:,1)==slices(i),1:63);
%         ADCslice = ADC(:,:,slices(i)); tempFeats = feats(indices(:,1)==slices(i),64:126);
        mask = zeros(size(T2slice));
%         mask = zeros(size(ADCslice));
        caMask = mask;featMap = mask;
        
        mask(tempInd(:,2))=1;
        caMask(tempInd(tempInd(:,3)==1,2)) = 1;
        for nF = [20]
%             nF = 32;
            featMap(mask==1) = tempFeats(:,nF)/max(tempFeats(:,nF));
%             overlayProbMap(T2slice,mask,featMap,1,caMask);
%             overlayProbMap(ADCslice,mask,featMap,1,caMask);
%             overlayProbMap(ADCslice,caMask,featMap,1,caMask);
              overlayProbMap(T2slice,caMask,featMap,1,caMask);
            pause;
        end
    end
end
close all
%% This is for niki bloch's data

addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath(genpath('C:\MatlabLibs\NIFTIread'));

dDir = 'C:\Data\forDrBloch\T2w\';
dDirD = 'C:\Data\forDrBloch\DCE\';

% 
nP = 'BOU137';
T2 = load_untouch_nii([dDir nP '_T2_N4.nii.gz']); T2 = T2.img;
maskV = load_untouch_nii([dDir nP '_Prostate_Or_Rsmp.nii.gz']); maskV = maskV.img;
% % 
% T2 = load_untouch_nii([dDirD nP '\ANGLEDPROSTATEWITHSP_7_Aff.nii.gz']); T2 = T2.img;
% maskV = load_untouch_nii([dDir nP '\ANGLEDPROSTATEWITHSP_3_Cancer.nii.gz']); maskV = maskV.img;
% ADC = load_untouch_nii([dDir num2str(nP) '\Pat' num2str(nP) '_T2_ADC_Aff.nii.gz']); ADC = ADC.img;

load(['C:\Data\forDrBloch\code\features\T2w\' nP '_T2_featsPixel.mat']);
% load(['C:\Data\forDrBloch\code\features\DCE\' nP '_DCE3_featsPixel.mat']);

slices = unique(indices(:,1));

for i = 1:length(slices)
    tempInd = indices(indices(:,1)==slices(i),:);
    if (max(tempInd(:,3))==1)
        mask = im2bw(maskV(:,:,slices(i))); 
        
        T2slice = T2(:,:,slices(i)); tempFeats = feats(indices(:,1)==slices(i),1:62);
%         ADCslice = ADC(:,:,slices(i)); tempFeats = feats(indices(:,1)==slices(i),64:126);
%         mask = zeros(size(ADCslice));
        caMask =  zeros(size(T2slice));featMap =  zeros(size(T2slice));
        inds = find(mask==1);
        caMask(inds(tempInd(:,3)==1)) = 1; 
        mask = double(I1); tempFeats = tempFeats(mask(inds)==1,:);
        for nF = 1:62
            nF
%             nF = 4;
            featMap(mask==1) = tempFeats(:,nF)/max(tempFeats(:,nF));
            overlayProbMap(T2slice,mask,featMap,1,caMask);
%             overlayProbMap(ADCslice,mask,featMap,1,caMask);
            pause;
        end
    end
end

