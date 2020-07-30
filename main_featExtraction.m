

%%
clear all;
close all;

addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath(genpath('C:\MatlabLibs\NIFTIread'));


[A,B,C ] = xlsread('patientsRSNA.xlsx');
pats = A(:,1); label = A(:,5:7);

%%

templatePath = 'C:\Code\DeltaRadiomics\dataRSNA\484\TP2';

template_T2 = mha_read_volume([templatePath filesep 'T2.mha']);
% template_ADC = mha_read_volume([templatePath filesep 'ADC.mha']);
% template_mask = mha_read_volume([mDir filesep studies{1,i} filesep 'PM.mha']);

template_caMask = mha_read_volume([templatePath filesep 'T2-label.mha']); template_caMask = logical(template_caMask);
template_mask = imdilate(template_caMask,strel('disk',18));

templateVolMasked_T2 = double(template_T2).*template_mask;
% templateVolMasked_ADC = double(template_ADC).*template_mask;
opts.temcancermasks = logical(template_caMask);
opts.docheck = false;
opts.dorescale = false;

datapath = 'C:\Code\DeltaRadiomics\dataRSNA\';
featStats = zeros(size(A,1), 75*8);
%%


for i =  9:size(A,1)
    
    disp(['reading study - ' num2str(A(i,1))]);
    
    T2 = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'T2.mha']);
    ADC = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'ADC.mha']);
    %     mask = mha_read_volume([mDir filesep studies{1,i} filesep 'PM.mha']);
    caMask = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'T2-label.mha']); caMask = logical(caMask);
    caMask_ADC = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'ADC-label.mha']); caMask_ADC = logical(caMask_ADC);
    %     header = mha_read_header([mDir filesep studies{1,i} filesep 'T2.mha']);
    mask = imdilate(caMask,strel('disk',18));
    
    if ~isequal(size(T2),size(caMask))
        disp('T2W mask and image are of different size!');
        pause;
    end
    
    if ~isequal(size(ADC),size(caMask_ADC))
        disp('ADC mask and image are of different size!');
        pause;
    end
    
    opts.incancermasks = logical(caMask);
    %
    if A(i,8)==1
        T2 = lpfbiascorr3(T2);
    end
    
    inputVolMask = double(T2).*mask;
    [~,stdMap,~] = int_stdn_landmarks(inputVolMask,templateVolMasked_T2,opts);
    close all;
    
    maskDil = dilateMaskVol(mask,30);
    inputVolMaskDil = double(T2).*maskDil;
    
    T2std = applystdnmap_rs(T2,stdMap);
    
    T2w_feats = []; ADC_feats = [];
    
    for j = 1:size(T2std,3)
        caMask_ = caMask(:,:,j);
        if max(caMask_(:))>0
            
            T2w_feats = [T2w_feats; computeTextureFeatures(T2std(:,:,j),caMask_)];
            %
        end
    end
    
    
    
    for j = 1:size(ADC,3)
        caMask_ADC_ = caMask_ADC(:,:,j);
        if max(caMask_ADC_(:))>0
            
            ADC_feats = [ADC_feats; computeTextureFeatures(ADC(:,:,j),caMask_ADC_)];
            %
        end
    end
    
    
    featStats(i,:) = [computeROIstatistics(T2w_feats) computeROIstatistics(ADC_feats)] ;
end


