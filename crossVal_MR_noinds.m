function varOutput = crossVal_MR_noinds(featsAvgN,label,classifierID,featSelectID)

addpath(genpath('C:\MatlabLibs\FEAST\FEAST-master'));
addpath('C:\MatlabLibs\libsvm-3.20\matlab');

% featsAvgN - normalized features matrix (n x f)
% label - label vector (n x 1) (1 and 0 labels)
% classifierID: 1 - linear regression, 2 - random forests 3 - svm
% featSelectID (feature selection method): 1 - mrmr, 2 - jmi, 3 - cmim

aucAll = [];  tblAll = [];
k = 3;runs = 150;B = struct;
maxF = 5; B_all = cell(maxF,runs); % select the max number of features
discreteSteps = 4; edges = [min(featsAvgN(:)):(max(featsAvgN(:))-min(featsAvgN(:)))/discreteSteps : max(featsAvgN(:))]; % to discretize for using the FEAST toolbox

varOutput = struct;

for f = 1:maxF
    topFeatFreq = [];
    for r = 1:runs
        temp1 = zeros(1,k);
        inds = crossvalind('Kfold',length(label),k);
        for i = 1 : k
            
            trainDataSelect = featsAvgN(inds~=i,:); trainLabelsSelect = label(inds~=i);
            trainDataSelect = discretize(trainDataSelect,edges);
            %             topFeats = mrmr_miq_d(trainDataSelect,trainLabelsSelect,f);
            
            switch featSelectID
                case 1
                    fsID = 'mrmr';
                case 2
                    fsID = 'jmi';
                case 3
                    fsID = 'cmim';
                case 4
                    fsID = 'none';
            end
            
            nf = f;
             maxFeats = maxF;
%             if classifierID == 2
%                 nf = f;
%                 maxFeats = maxF^2;
%             else
%                 nf = f;
%                 maxFeats = maxF;
%             end
            
            tempF = zeros(1,maxFeats);
            
            
            if strcmp(fsID,'none')
                topFeats = 1:nf;
            else
                topFeats = feast(fsID,nf,trainDataSelect,trainLabelsSelect);
            end
            
            tempF(1,(maxFeats+1-nf):maxFeats) = topFeats;
            topFeatFreq = [topFeatFreq;tempF];
            
            trainData = featsAvgN(inds~=i,topFeats);
            
            trainLabels = label(inds~=i);
            
            
            testData = featsAvgN(inds==i,topFeats);
            testLabels = label(inds==i);
            
            reorder = randperm(size(trainData,1));
            trainData = trainData(reorder,:);
            trainLabels = trainLabels(reorder,:);
            
            B(i).f = topFeats;
            switch classifierID
                case 1 % logistic regression classifier
                    
                    b = glmfit(trainData,trainLabels,'binomial');
                    B(i).b = b;
                    posteriors = glmval(b,testData,'logit');
                    [X,Y,~,temp1(i)] = perfcurve(testLabels,posteriors(:,1),1);disp(['auc of the fold - ' num2str(temp1(i))]);
                    %                     obj = fitcdiscr(trainData,trainLabels);B(i).b = obj;
                    %                     [~,posteriors] = predict(obj,testData);
%                 case 2 % random forests classifier
%                     
%                     B(i).b = TreeBagger(50,trainData,trainLabels,'OOBPrediction','on','Method','classification','NVarToSample','all');
%                     [~,posteriors] = predict(B(i).b,testData);
                
               case 2 % RF
                    B(i).b = fitctree(trainData,trainLabels);
                    [~,posteriors] = predict(B(i).b,testData);
                    [X,Y,~,temp1(i)] = perfcurve(testLabels,posteriors(:,2),1);disp(['auc of the fold - ' num2str(temp1(i))]);
                    
                case 3 % SVM classifier
                    B(i).b = svmtrain(trainLabels,trainData, '-t 2 -c 100 -b 1 -q'); %#ok<SVMTRAIN>
                    [~,~,posteriors] = svmpredict(testLabels,testData,B(i).b , '-b 1');
                    [X,Y,~,temp1(i)] = perfcurve(testLabels,posteriors(:,1),1);disp(['auc of the fold - ' num2str(temp1(i))]);
            end
            
            if length(unique(testLabels))<=1
  
                temp1(i) = 0;
            end
            
            
            %             posteriorsAll = glmval(B(i).b,featsAvgN(:,B(i).f),'logit');
            %
            %             testLabelsAll = label-1;
            %
            %             [~,~,~,AUC] = perfcurve(testLabelsAll,posteriorsAll(:,1),1);
            %
            %             if AUC>.7
            %                 pause;
            %             end
            
        end
        B_all{f,r} = B;
        C = find(temp1==max(temp1));
        aucAll = [aucAll; f r mean(temp1) C(1) max(temp1) std(temp1)];
    end
    A = topFeatFreq(topFeatFreq~=0);
    tbl = tabulate(A);
    
    if size(tbl,1)~=size(featsAvgN,2)
        temp2 = size(tbl,1)+1:size(featsAvgN,2);
        tbl = cat(1,tbl,[ temp2' zeros(length(temp2),2)]);
    end
    
    tblAll = cat(2,tblAll,tbl(:,2));
end

[~,y] = sort(aucAll(:,3),'descend');
aucAll = aucAll(y,:);

varOutput.aucAll = aucAll;
varOutput.tblAll = tblAll;
varOutput.B_all = B_all;
