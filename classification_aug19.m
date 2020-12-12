addpath(genpath('C:\MatlabLibs\FEAST\FEAST-master'));
addpath('C:\MatlabLibs\libsvm-3.20\matlab');
addpath(genpath('C:\Code\General'));
addpath(genpath('C:\Matlab code from repo'));
addpath('C:\MatlabLibs\MATLAB-SMOTE-master');
addpath('C:\Code\risk_stratification_DL_Amogh');
addpath('D:\extraCapFat');

%%

clear all;
load('D:\extraCapFat\feature_matrix\SUPP\feat_stats_T2std_Ecmask.mat');
load('D:\extraCapFat\feature_matrix\SUPP\SUPP_label.mat');
X_train = feat_stats; Y_train = label;


load('D:\extraCapFat\feature_matrix\UH\UH_feature_stats_T2std_ecMask.mat');
load('D:\extraCapFat\feature_matrix\UH\UH_label.mat');

X_test = feat_stats; Y_test = label;


%%

% clear all;
% close all;
clc;

X_train = feat_stats_T2;
Y_train = label;

featsAvg = simplewhiten(X_train);
featsAvg_n = featsAvg;
pVals = zeros(size(Y_train));

for i = 1:size(featsAvg,2)
    featsAvg_n(:,i) = rescale(featsAvg(:,i));
%     pVals(i) = ranksum(featsAvgN(Y_train==0,i),featsAvgN(Y_train==1,i));
end

for i = 1:size(featsAvg_n,2)
%     featsAvgN(:,i) = rescale_range(featsAvg(:,i),-1,1);
    pVals(i) = ranksum(featsAvg_n(Y_train==0,i),featsAvg_n(Y_train==1,i));
end

selectF = find(pVals<=0.05);
% 
featsAvg_sf = featsAvg_n(:,selectF);
label = Y_train;
featsAvgN = featsAvg_sf;
% [augFeats,augLabel] = ADASYN(featsAvg_sf,Y_train);
% 
% featsAvgN = [featsAvg_sf;augFeats];
% label = [Y_train;augLabel];

%% 

for cID = 3:3
    for fsID = 1:1
        disp(['running ' num2str(cID) '-' num2str(fsID)]);
        output = crossVal_(featsAvgN,label,cID,fsID);
        save(['models\trainData_' num2str(cID) '_' num2str(fsID) '.mat'],'output');
    end 
end
 
%% prepare test data

featsAvgT = simplewhiten_rs(X_test,X_train);

featsAvgT_n = zeros(size(featsAvgT));
for i = 1:size(featsAvgT,2)
    featsAvgT_n(:,i) = rescale_range_rs(featsAvgT(:,i),featsAvg(:,i));
%     pVals(i) = ranksum(featsAvgN(Y_train==0,i),featsAvgN(Y_train==1,i));
end
% featsAvg = simplewhiten(X_test);
 featsAvgNT = featsAvgT_n(:,selectF);

% 
% for i = 1:size(featsAvgNt,2)
%     featsAvgN(:,i) = rescale_range(featsAvgNt(:,i),-1,1,allFeats(1:size(X_train,1),i));
% end

label = Y_test;

%%
addpath('C:\MatlabLibs\libsvm-3.20\matlab');
aucMatrix =[];
clc
aucPlot = struct;
for cID = 3:3
    for fsID = 1:1
        load(['D:\extraCapFat\models\trainData_' num2str(cID) '_' num2str(fsID) '.mat']);
        disp(num2str(output.aucAll(1:5,[3 6])));
        temp = input('enter the row aucAll  - ');
        x = output.aucAll(temp,1:2);
        B_best = output.B_all{x(1),x(2)};
        B_best = B_best(output.aucAll(temp,4));
        testLabels = label;
        
        switch cID
            case 1
                  posteriors = glmval(B_best.b,featsAvgNT(:,B_best.f),'logit');
%                 [predLabel,posteriors,~] = predict(B_best.b,featsAvgN(:,B_best.f));
            case 2
                [~,posteriors] = predict(B_best.b,featsAvgNT(:,B_best.f));
            case 3
                [predLabel,~,posteriors] = svmpredict(testLabels,featsAvgNT(:,B_best.f),B_best.b , '-b 1');
        end
        
        
        [X,Y,~,AUC,optPT,~,~] = perfcurve(testLabels,posteriors(:,1),1);disp(['auc  - ' num2str(AUC)]);
        aucPlot(fsID).X = X;
        aucPlot(fsID).Y = Y;
        pause;
        aucMatrix = [ aucMatrix; output.aucAll(temp,[3 6]) AUC];
    end
end

