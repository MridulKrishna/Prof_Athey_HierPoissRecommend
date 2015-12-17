clear
clc

dir = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/';
dirOut = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/';
dirCodes = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/';

%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact5';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-seed5-scfact5';
name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact10';

factor = 1;

% Loads the output from the HPFO algorithm
betaRate = dlmread(strcat(dirOut,name,'/hbeta_rate.tsv'),'\t');
betaShape = dlmread(strcat(dirOut,name,'/hbeta_shape.tsv'),'\t');
beta = dlmread(strcat(dirOut,name,'/hbeta.tsv'),'\t');
thetaRate = dlmread(strcat(dirOut,name,'/htheta_rate.tsv'),'\t');
thetaShape = dlmread(strcat(dirOut,name,'/htheta_shape.tsv'),'\t');
theta = dlmread(strcat(dirOut,name,'/htheta.tsv'),'\t');
sigmaRate = dlmread(strcat(dirOut,name,'/hsigma_rate.tsv'),'\t');
sigmaShape = dlmread(strcat(dirOut,name,'/hsigma_shape.tsv'),'\t');
sigma = dlmread(strcat(dirOut,name,'/hsigma.tsv'),'\t');
rhoRate = dlmread(strcat(dirOut,name,'/hrho_rate.tsv'),'\t');
rhoShape = dlmread(strcat(dirOut,name,'/hrho_shape.tsv'),'\t');
rho = dlmread(strcat(dirOut,name,'/hrho.tsv'),'\t');

% Loads user observables
obsUser = dlmread(strcat(dir,'obsUser.tsv'),'\t');
U = size(obsUser,1);

% Id map for users
u_c2n = containers.Map('KeyType','uint64','ValueType','uint64');
u_n2c = containers.Map('KeyType','uint64','ValueType','uint64');
for u = 1:U
        u_c2n(obsUser(u,1))=u;
        u_n2c(u)=obsUser(u,1);
    end

%sigmaOrd = sortrows(sigma,factor+2);

match = zeros(size(sigma,1),9+size(sigma,2));
matchComp = zeros(size(sigma,1),37);

for i = 1:size(sigma,1)
    hhid = sigma(i,2);
    hhid
    match(i,1) = hhid;
    matchComp(i,1) = hhid;
    matchComp(i,2:end) = obsUser(u_c2n(hhid),2:end);
    
    match(i,2) = matchComp(i,2);
    match(i,3) = matchComp(i,3);
    match(i,4) = matchComp(i,4);
    match(i,5) = 1*matchComp(i,5)+2*matchComp(i,6)+3*matchComp(i,7)+4*matchComp(i,8)+5*matchComp(i,9)+6*matchComp(i,10);
    match(i,6) = 1*matchComp(i,12)+2*matchComp(i,13)+3*matchComp(i,14)+4*matchComp(i,15)+5*matchComp(i,16);
    match(i,7) = 1*matchComp(i,18)+2*matchComp(i,19)+3*matchComp(i,20)+4*matchComp(i,21)+5*matchComp(i,22)+7*matchComp(i,23);
    match(i,8) = 1*matchComp(i,25)+2*matchComp(i,26)+3*matchComp(i,27)+4*matchComp(i,28)+5*matchComp(i,29);
    match(i,9) = 1*matchComp(i,31)+2*matchComp(i,32)+3*matchComp(i,33)+4*matchComp(i,34)+5*matchComp(i,35)+6*matchComp(i,36)+7*matchComp(i,37);
   
    match(i,10:end) = sigma(i,:);
end