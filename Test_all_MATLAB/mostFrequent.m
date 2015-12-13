clear
clc

dir = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/';
dirOut = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/';
dirCodes = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/';

name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact5';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-seed5-scfact5';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact10';

factor = 23;

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

% Loads the brand and flavor codes
fileID = fopen(strcat(dirCodes,'/flavors.tsv'));
flavors = textscan(fileID,'%d8 %s','Delimiter','\t');
fclose(fileID);

fileID = fopen(strcat(dirCodes,'/brands.tsv'));
brands = textscan(fileID,'%d8 %s','Delimiter','\t');
fclose(fileID);

betaOrd = sortrows(beta,factor+2);

topNum = zeros(20,6);
top = cell(20,7);

for i = 1:20
    str = num2str(betaOrd(size(beta,1)+1-i,2),'%u');
    top{i} = str;
    
    if length(str) == 12
        topNum(i,1) = str2double(str(1:1));
        topNum(i,2) = str2double(str(2:3));
        topNum(i,3) = str2double(str(4:5));
        topNum(i,4) = str2double(str(6:6));
        topNum(i,5) = str2double(str(7:7));
        topNum(i,6) = str2double(str(8:12));
    else
        topNum(i,1) = str2double(str(1:2));
        topNum(i,2) = str2double(str(3:4));
        topNum(i,3) = str2double(str(5:6));
        topNum(i,4) = str2double(str(7:7));
        topNum(i,5) = str2double(str(8:8));
        topNum(i,6) = str2double(str(9:13));        
    end
    
    top{i,1} = brands{2}{topNum(i,1)};
    if topNum(i,2) == 0
        top{i,2} = '';
    else
        top{i,2} = flavors{2}{topNum(i,2)};
    end
    if topNum(i,3) == 0
        top{i,3} = '';
    else
        top{i,3} = flavors{2}{topNum(i,3)};
    end
    
    if topNum(i,4) == 1
        top{i,4} = 'low fat';
    end
    if topNum(i,5) == 1
        top{i,5} = 'no fat';
    end
    top{i,6} = topNum(i,6);
    top{i,7} = betaOrd(size(beta,1)+1-i,factor+2);
end