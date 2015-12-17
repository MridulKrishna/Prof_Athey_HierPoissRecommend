clear
clc

dir = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/';
dirOut = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb';

% Loads train, test, and validation sets
train = dlmread(strcat(dir,'train.tsv'),'\t');
test = dlmread(strcat(dir,'test.tsv'),'\t');
validation = dlmread(strcat(dir,'validation.tsv'),'\t');

strain = size(train,1);
stest = size(test,1);
svalidation = size(validation,1);

nRuns = 13;
names = cell(nRuns,1);
names{1,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb';
names{2,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact5';
names{3,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact10';
names{4,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std';
names{5,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std-scfact5';
names{6,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std-scfact10';
names{7,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones';
names{8,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones-scfact5';
names{9,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones-scfact10';
names{10,1} = 'n2599-m275-k0-uc36-ic50-batch-hier-vb';
names{11,1} = 'n2599-m275-k111-uc0-ic0-batch-hier-vb';
names{12,1} = 'n2599-m275-k25-uc36-ic0-batch-hier-vb';
names{13,1} = 'n2599-m275-k25-uc0-ic50-batch-hier-vb';

Ks = 25*ones(nRuns,1);
Ks(10) = 0;
Ks(11) = 111;
UCs = 36*ones(nRuns,1);
UCs(11) = 0;
UCs(13) = 0;
ICs = 50*ones(nRuns,1);
ICs(11) = 0;
ICs(12) = 0;

testll = zeros(nRuns,1);
validationll = zeros(nRuns,1);

for n = 1:1
    name = names{n,1};

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
    xiRate = dlmread(strcat(dirOut,name,'/thetarate_rate.tsv'),'\t');
    xiShape = dlmread(strcat(dirOut,name,'/thetarate_shape.tsv'),'\t');
    xi = dlmread(strcat(dirOut,name,'/thetarate.tsv'),'\t');
    etaRate = dlmread(strcat(dirOut,name,'/betarate_rate.tsv'),'\t');
    etaShape = dlmread(strcat(dirOut,name,'/betarate_shape.tsv'),'\t');
    eta = dlmread(strcat(dirOut,name,'/betarate.tsv'),'\t');
    
    % Loads observed item and user characteristics
    obsItem = dlmread(strcat(dir,'obsItem.tsv'),'\t');
    obsUser = dlmread(strcat(dir,'obsUser.tsv'),'\t');
    
    % Stores the number of factors, users, and items
    U = size(obsUser,1);
    I = size(obsItem,1);
    K = size(beta,2)-2;
    L = size(sigma,2)-2;
    M = size(rho,2)-2;
    
    % Sorts matrices by item and user id
    beta = sortrows(beta,2);
    theta = sortrows(theta,2);
    sigma = sortrows(sigma,2);
    rho = sortrows(rho,2);
    
    eta = sortrows(eta,2);
    xi = sortrows(xi,2);
    etaRate = sortrows(etaRate,2);
    xiRate = sortrows(xiRate,2);
    etaShape = sortrows(etaShape,2);
    xiShape = sortrows(xiShape,2);

    obsItem = sortrows(obsItem,1);
    obsUser = sortrows(obsUser,1);
    
    % Check that all ids match
    if (UCs(n) > 0 && Ks(n) > 0 )
    for i = 1:size(obsItem,1)
        if beta(i,2) ~= obsItem(i,1) 
            disp(strcat('beta: ',num2str(i)));
        end

        if rho(i,2) ~= obsItem(i,1) 
            disp(strcat('rho: ',num2str(i)));
        end
        
        if eta(i,2) ~= obsItem(i,1) 
            disp(strcat('eta: ',num2str(i)));
        end
        
        if etaRate(i,2) ~= obsItem(i,1) 
            disp(strcat('etaRate: ',num2str(i)));
        end
        
        if etaShape(i,2) ~= obsItem(i,1) 
            disp(strcat('etaShape: ',num2str(i)));
        end
    end
    end
    if ( ICs(n) > 0 && Ks(n) > 0 )
    for i = 1:size(obsUser,1)
        if theta(i,2) ~= obsUser(i,1) 
            disp(strcat('theta: ',num2str(i)));
        end

        if sigma(i,2) ~= obsUser(i,1) 
            disp(strcat('sigma: ',num2str(i)));
        end
        
        if xi(i,2) ~= obsUser(i,1) 
            disp(strcat('xi: ',num2str(i)));
        end
        
        if xiShape(i,2) ~= obsUser(i,1) 
            disp(strcat('xiShape: ',num2str(i)));
        end
        
        if xiRate(i,2) ~= obsUser(i,1) 
            disp(strcat('xiRate: ',num2str(i)));
        end
    end
    end

    disp('Finished checking that they match')
    
    % Saves the keys
    if ICs(n) > 0
        i_c2n = containers.Map('KeyType','uint64','ValueType','uint64');
        i_n2c = containers.Map('KeyType','uint64','ValueType','uint64');
        for i = 1:I
            i_c2n(obsItem(i,1))=i;
            i_n2c(i)=obsItem(i,1);
        end
    end
    if UCs(n) > 0
        u_c2n = containers.Map('KeyType','uint64','ValueType','uint64');
        u_n2c = containers.Map('KeyType','uint64','ValueType','uint64');
        for u = 1:U
            u_c2n(obsUser(u,1))=u;
            u_n2c(u)=obsUser(u,1);
        end
    end
    
    % Computes likelihood of validation set
    llValidation = 0;
    
    for i = 1:svalidation
        validation(i,2);
        
        rate = 0;
        if Ks(n) > 0
            itheta = theta(u_c2n(validation(i,1)),3:end);
            ibeta = beta(i_c2n(validation(i,2)),3:end);
            rate = rate + itheta*ibeta';
        end
        if ICs(n) > 0
            isigma = sigma(u_c2n(validation(i,1)),3:end);
            x = obsItem(i_c2n(validation(i,2)),2:end);
            ietaRate = etaRate(i_c2n(validation(i,2)),3);
            ietaShape = etaShape(i_c2n(validation(i,2)),3);
            rate = rate + isigma*x'*ietaRate/(ietaShape-1);
%             rate = rate + isigma*x';
        end
        if UCs(n) > 0
            irho = rho(i_c2n(validation(i,2)),3:end);
            w = obsUser(u_c2n(validation(i,1)),2:end);
            ixiRate = xiRate(u_c2n(validation(i,1)),3);
            ixiShape = xiShape(u_c2n(validation(i,1)),3);
            rate = rate + irho*w'*ixiRate/(ixiShape-1);
%             rate = rate + irho*w';
        end
        
        rating = validation(i,3);

        if (rating < 150)
            indLL = rating*log(rate)-rate-log(factorial(rating));
        else
            indLL = rating*log(rate)-rate-logfact(rating);
        end
            
%         if i < 100
%             disp(indLL)
%             disp(rate)
%             disp(rating)
%         end
        
        llValidation = llValidation + indLL;
        
        if ( i == 3942 )
            disp('a')
            disp(i_c2n(validation(i,2)))
            disp(ietaRate)
            disp(ietaShape)
            disp(itheta.*ibeta)
            disp('isigma')
            m1 = [isigma.*x*ietaRate/(ietaShape-1);isigma;x;ietaRate/(ietaShape-1)*ones(1,size(x,2));...
                ietaRate*ones(1,size(x,2));ietaShape*ones(1,size(x,2))]'
            disp('irho')
            m2 = [irho.*w*ixiRate/(ixiShape-1);irho;w;ixiRate/(ixiShape-1)*ones(1,size(w,2))]'
            i
            disp(validation(i,1))
            disp(validation(i,2))
            disp(strcat('s: ',num2str(rate),' y: ',num2str(rating),' ll: ',num2str(indLL)));
            asf = 1+2;
        end
        
%         if i < 100
%             format long
%             disp(indLL)
%             format long
%             disp(rate)
%             format long
%             disp(rating)
%         end
    end
  
    mllValidation = llValidation/svalidation;
    
    validationll(n) = mllValidation;
    
    % Computes likelihood of test set
    lltest = 0;
    
%     mat = zeros(stest,5);
    
    for i = 1:stest
        test(i,2);
        
        rate = 0;
        if Ks(n) > 0
            itheta = theta(u_c2n(test(i,1)),3:end);
            ibeta = beta(i_c2n(test(i,2)),3:end);
            rate = rate + itheta*ibeta';
        end
        if ICs(n) > 0
            isigma = sigma(u_c2n(test(i,1)),3:end);
            x = obsItem(i_c2n(test(i,2)),2:end);
            ietaRate = etaRate(i_c2n(test(i,2)),3);
            ietaShape = etaShape(i_c2n(test(i,2)),3);
            rate = rate + isigma*x'*ietaRate/(ietaShape-1);
        end
        if UCs(n) > 0
            irho = rho(i_c2n(test(i,2)),3:end);
            w = obsUser(u_c2n(test(i,1)),2:end);
            ixiRate = xiRate(u_c2n(test(i,1)),3);
            ixiShape = xiShape(u_c2n(test(i,1)),3);
            rate = rate + irho*w'*ixiRate/(ixiShape-1);
        end

        rating = test(i,3);

        if (rating < 150)
            indLL = rating*log(rate)-rate-log(factorial(rating));
        else
            indLL = rating*log(rate)-rate-logfact(rating);
        end
  
        lltest = lltest + indLL;
        
%         disp(test(i,1))
%         disp(test(i,2))
%         disp(strcat('s: ',num2str(rate),' y: ',num2str(rating),' ll: ',num2str(indLL)));
        
%         if i < 100
%             format long
%             disp(indLL)
%             format long
%             disp(rate)
%             format long
%             disp(rating)
%         end
        
%         mat(i,1) = rating;
%         mat(i,2) = rate;
%         mat(i,3) = test(i,1);
%         mat(i,4) = test(i,2);
%         mat(i,5) = indLL;
    end
   
    mlltest = lltest/stest;
    
    testll(n) = mlltest;
    
end

% mat = sortrows(mat,[1,2]);
% mat2 = dlmread(strcat(dirOut,name,'/likelihoodsTest.tsv'),'\t');
% mat2 = sortrows(mat2,[1,2]);
% diff = mat - mat2;
% todo = [mat,mat2,diff];
% for i = 1:7998
%     if diff(i,1) ~= 0
%         disp(strcat('col1: ',num2str(i)))
%     end
%     if diff(i,2) > 1e-3
%         disp(strcat('col2: ',num2str(i)))
%     end
%     if diff(i,5) > 1e-1
%         disp(strcat('col5: ',num2str(i)))
%     end
% end
