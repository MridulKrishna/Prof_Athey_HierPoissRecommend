clear
clc

dir = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/';
dirOut = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/';
dirCodes = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/';

%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact5';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-seed5-scfact5';
name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact10';

factor = 9;

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

thetaOrd = sortrows(theta,factor+2);

top = zeros(20,10);
topComp = zeros(20,37);

for i = 1:20
    hhid = thetaOrd(size(theta,1)+1-i,2);
    hhid
    top(i,1) = hhid;
    topComp(i,1) = hhid;
    topComp(i,2:end) = obsUser(u_c2n(hhid),2:end);
    
    top(i,2) = topComp(i,2);
    top(i,3) = topComp(i,3);
    top(i,4) = topComp(i,4);
    top(i,5) = 1*topComp(i,5)+2*topComp(i,6)+3*topComp(i,7)+4*topComp(i,8)+5*topComp(i,9)+6*topComp(i,10);
    top(i,6) = 1*topComp(i,12)+2*topComp(i,13)+3*topComp(i,14)+4*topComp(i,15)+5*topComp(i,16);
    top(i,7) = 1*topComp(i,18)+2*topComp(i,19)+3*topComp(i,20)+4*topComp(i,21)+5*topComp(i,22)+7*topComp(i,23);
    top(i,8) = 1*topComp(i,25)+2*topComp(i,26)+3*topComp(i,27)+4*topComp(i,28)+5*topComp(i,29);
    top(i,9) = 1*topComp(i,31)+2*topComp(i,32)+3*topComp(i,33)+4*topComp(i,34)+5*topComp(i,35)+6*topComp(i,36)+7*topComp(i,37);
   
    top(i,10) = thetaOrd(size(theta,1)+1-i,factor+2);
end