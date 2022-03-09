

%%
clear all;
close all;

addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath(genpath('C:\MatlabLibs\NIFTIread'));


[A,B,C ] = xlsread('D:\extraCapFat\new_UH_RP_datasheet_Feb6_20.xlsx');
pats = {B{2:end,1}}; label = A(:,1);

%%

templatePath = 'G:\.shortcut-targets-by-id\1CPvrtnpodCQ7_BNcO6dTfs1xhZNu3-kJ\Template';

template_T2 = niftiread([templatePath filesep 'T2W.nii']);
template_ADC = niftiread([templatePath filesep 'ADC.nii']);
template_mask = niftiread([templatePath filesep 'PM.nii']);

% template_caMask = mha_read_volume([templatePath filesep 'T2-label.mha']); template_caMask = logical(template_caMask);
% template_mask = imdilate(template_caMask,strel('disk',18));

templateVolMasked_T2 = double(template_T2).*double(template_mask);
% templateVolMasked_ADC = double(template_ADC).*template_mask;
% opts.temcancermasks = logical(template_caMask);
opts.docheck = false;
opts.dorescale = false;

datapath = 'G:\.shortcut-targets-by-id\1lRjm8H415GVMV9hWhB-yTApVJlWcKIsZ\UH_45PAT_PPF';

%%


for i =  1:length(pats)
    
    disp(['reading study - ' pats{i}]);
    
    T2 = mha_read_volume([datapath filesep pats{i}(13:end) '_T2.mha']);
    ADC = mha_read_volume([datapath filesep pats{i}(13:end) '_ADC.mha']);
    mask = mha_read_volume([datapath filesep pats{i}(13:end) '_T2_prostate_label.mha']);
    hdr = mha_read_header([datapath filesep pats{i}(13:end) '_T2.mha']);
    %     caMask = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'T2-label.mha']); caMask = logical(caMask);
    %     caMask_ADC = mha_read_volume([datapath filesep num2str(A(i,1)) filesep 'TP2' filesep 'ADC-label.mha']); caMask_ADC = logical(caMask_ADC);
    %     header = mha_read_header([mDir filesep studies{1,i} filesep 'T2.mha']);
    %     mask = imdilate(caMask,strel('disk',18));
    
    %     if ~isequal(size(T2),size(caMask))
    %         disp('T2W mask and image are of different size!');
    %         pause;
    %     end
    %
    %     if ~isequal(size(ADC),size(caMask_ADC))
    %         disp('ADC mask and image are of different size!');
    %         pause;
    %     end
    
    %     opts.incancermasks = logical(caMask);
    %
    %     if A(i,8)==1
%     T2 = lpfbiascorr3(T2);
    %     end
    
    inputVolMask = double(T2).*double(mask);
    [~,stdMap,~] = int_stdn_landmarks(inputVolMask,templateVolMasked_T2,opts);
    close all;
    
    %     maskDil = dilateMaskVol(mask,30);
    %     inputVolMaskDil = double(T2).*maskDil;
    
    T2std = applystdnmap_rs(T2,stdMap);
    
    mha_write_volume([datapath filesep pats{i}(13:end) '_T2_std.mha'],T2std,hdr.PixelDimensions,hdr.Offset);
end

