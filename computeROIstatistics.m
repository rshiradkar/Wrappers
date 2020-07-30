function featStats = computeROIstatistics(feats)

% feats - n x d, n is the pixels in ROI, d is number of features
% featStats - computes mean, median, std, 25th, 75th percentile, skewness
% kurtosis; so it is 1 x d*7

% featStats = [mean(feats,1) median(feats,1) std(feats,0,1) prctile(feats,25,1)...
%     prctile(feats,75,1) skewness(feats,1,1) kurtosis(feats,1,1)];

featStats = [mean(feats,1) std(feats,0,1) skewness(feats,1,1) kurtosis(feats,1,1)];
    
