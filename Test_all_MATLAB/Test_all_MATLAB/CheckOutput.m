clear

% Loads the output from the HPFO algorithm
betaMeans = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hbeta.tsv','\t');
thetaMeans = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/htheta.tsv','\t');
sigmaMeans = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hsigma.tsv','\t');
rhoMeans = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hrho.tsv','\t');
% xiRate = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/thetarate_shape.tsv','\t');
% xiShape = dlmread('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/thetarate_rate.tsv','\t');
obsItem = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/obsItem.tsv','\t');
obsUser = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/obsUser.tsv','\t');

% Sorts matrices by item and user id
betaMeans = sortrows(betaMeans,2);
thetaMeans = sortrows(thetaMeans,2);
sigmaMeans = sortrows(sigmaMeans,2);
rhoMeans = sortrows(rhoMeans,2);
% xiRate = sortrows(xiRate,2);
% xiShape = sortrows(xiShape,2);

obsItem = sortrows(obsItem,1);
obsUser = sortrows(obsUser,1);

% Check that all ids match
for i = 1:size(obsItem,1)
    if betaMeans(i,2) ~= obsItem(i,1) 
        disp(strcat('betaMeans: ',i));
    end
    
    if rhoMeans(i,2) ~= obsItem(i,1) 
        disp(strcat('rhoMeans: ',i));
    end
end
for i = 1:size(obsUser,1)
    if thetaMeans(i,2) ~= obsUser(i,1) 
        disp(strcat('thetaMeans: ',i));
    end
    
    if sigmaMeans(i,2) ~= obsUser(i,1) 
        disp(strcat('sigmaMeans: ',i));
    end
end