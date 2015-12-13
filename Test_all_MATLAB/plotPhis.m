clear

dir = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/';
dirOut = '~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/output/';
%name = 'n2599-m275-k25-uc36-ic50-batch-hier-vb';

names = cell(10,1);
names{1,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb';
names{2,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact5';
names{3,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-scfact10';
names{4,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std';
names{5,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std-scfact5';
names{6,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-std-scfact10';
names{7,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones';
names{8,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones-scfact5';
names{9,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones-scfact10';
names{10,1} = 'n2599-m275-k25-uc36-ic50-batch-hier-vb-ones-scfact100';

for n = 1:10
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

    % Loads observed item and user characteristics
    obsItem = dlmread(strcat(dir,'obsItem.tsv'),'\t');
    obsUser = dlmread(strcat(dir,'obsUser.tsv'),'\t');

    % Sorts matrices by item and user id
    beta = sortrows(beta,2);
    theta = sortrows(theta,2);
    sigma = sortrows(sigma,2);
    rho = sortrows(rho,2);

    obsItem = sortrows(obsItem,1);
    obsUser = sortrows(obsUser,1);

    % Check that all ids match
    for i = 1:size(obsItem,1)
        if beta(i,2) ~= obsItem(i,1) 
            disp(strcat('beta: ',i));
        end

        if rho(i,2) ~= obsItem(i,1) 
            disp(strcat('rho: ',i));
        end
    end
    for i = 1:size(obsUser,1)
        if theta(i,2) ~= obsUser(i,1) 
            disp(strcat('theta: ',i));
        end

        if sigma(i,2) ~= obsUser(i,1) 
            disp(strcat('sigma: ',i));
        end
    end

    disp('Finished checking that they match')

    U = size(obsUser,1);
    I = size(obsItem,1);
    K = size(beta,2)-2;
    L = size(sigma,2)-2;
    M = size(rho,2)-2;

    % Inefficient, but can calculate some additional summary statistics
    % such as percentiles
    
% % Generates random sample
%     P = U*I;
%     ps = randperm(U*I,P);
% 
%     phi = zeros(P,K+L+M);
% 
%     for p = 1:P
%         u = mod(ps(p),U)+1;
%         i = ceil(ps(p)/U);
%     %         disp(ps(p))
%     %     disp(u)
%     %     disp(i)
%         phi(p,1:K) = theta(u,3:end).*beta(i,3:end);
%         phi(p,K+1:K+L) = sigma(u,3:end).*obsItem(i,2:end);
%         phi(p,K+L+1:K+L+M) = obsUser(u,2:end).*rho(i,3:end);
% 
% %         phi(p,:) = phi(p,:)/sum(phi(p,:));
%     end
% 
%     phimean1 = mean(phi,1);
%     phistd1 = std(phi,0,1);
%     upper = phimean1+phistd1;
%     lower = phimean1-phistd1;
%     p50 = prctile(phi,50,1);
%     p95 = prctile(phi,95,1);
%     p75 = prctile(phi,75,1);
%     p25 = prctile(phi,25,1);
%     p5 = prctile(phi,5,1);   


    % Means and standard deviations for population
    thetameans = zeros(K,1);
    betameans = zeros(K,1);
    sigmameans = zeros(L,1);
    rhomeans = zeros(M,1);
    obsItemMeans = zeros(L,1);
    obsUserMeans = zeros(M,1);
    thetass = zeros(K,1);
    betass = zeros(K,1);
    sigmass = zeros(L,1);
    rhoss = zeros(M,1);
    obsItemSs = zeros(L,1);
    obsUserSs = zeros(M,1);
    
    phimean = zeros(K+L+M,1);
    phistd = zeros(K+L+M,1);
    
    for k = 1:K
        for i = 1:I
            betameans(k) = betameans(k)+1/I*beta(i,k+2);
            betass(k) = betass(k)+beta(i,k+2)^2;
        end
            
        for u = 1:U
            thetameans(k) = thetameans(k)+1/U*theta(u,k+2);
            thetass(k) = thetass(k)+theta(u,k+2)^2;
        end
        
        phimean(k) = thetameans(k)*betameans(k);
        phistd(k) = sqrt(1/(I*U)*thetass(k)*betass(k)-thetameans(k)^2*betameans(k)^2);
    end
    
    for l = 1:L
        for u = 1:U
            sigmameans(l) = sigmameans(l)+1/U*sigma(u,l+2);
            sigmass(l) = sigmass(l)+sigma(u,l+2)^2;
        end
        
        for i = 1:I
            obsItemMeans(l) = obsItemMeans(l)+1/I*obsItem(i,l+1);
            obsItemSs(l) = obsItemSs(l)+obsItem(i,l+1)^2;
        end
        
        phimean(K+l) = sigmameans(l)*obsItemMeans(l);
        phistd(K+l) = sqrt(1/(I*U)*sigmass(l)*obsItemSs(l)-sigmameans(l)^2*obsItemMeans(l)^2);
    end
    
    for m = 1:M
        for u = 1:U
            obsUserMeans(m) = obsUserMeans(m)+1/U*obsUser(u,m+1);
            obsUserSs(m) = obsUserSs(m)+obsUser(u,m+1)^2;
        end
        
        for i = 1:I
            rhomeans(m) = rhomeans(m)+1/I*rho(i,m+2);
            rhoss(m) = rhoss(m)+rho(i,m+2)^2;
        end
        
        phimean(K+L+m) = obsUserMeans(m)*rhomeans(m);
        phistd(K+L+m) = sqrt(1/(I*U)*obsUserSs(m)*rhoss(m)-obsUserMeans(m)^2*rhomeans(m)^2);
    end

%     subplot(2,1,1)
%     plot(1:(K+L+M),phimean1,1:(K+L+M),phistd1,25.5*[1,1],[0,max(phistd)],75.5*[1,1],[0,max(phistd)],...
%         'r',26.5*[1,1],[0,max(phistd)],'k',28.5*[1,1],[0,max(phistd)],'k',46.5*[1,1],[0,max(phistd)],'k',...
%         76.5*[1,1],[0,max(phistd)],'k',78.5*[1,1],[0,max(phistd)],'k',...
%         84.5*[1,1],[0,max(phistd)],'k',90.5*[1,1],[0,max(phistd)],'k',97.5*[1,1],[0,max(phistd)],'k',...
%         103.5*[1,1],[0,max(phistd)],'k',111.5*[1,1],[0,max(phistd)],'k')
%     legend('mean of phi','standard deviation of phi')
%     xlabel('factor')
%     subplot(2,1,2)
    plot(1:(K+L+M),phimean,1:(K+L+M),phistd,25.5*[1,1],[0,max(phistd)],75.5*[1,1],[0,max(phistd)],...
        'r',26.5*[1,1],[0,max(phistd)],'k',28.5*[1,1],[0,max(phistd)],'k',46.5*[1,1],[0,max(phistd)],'k',...
        76.5*[1,1],[0,max(phistd)],'k',78.5*[1,1],[0,max(phistd)],'k',...
        84.5*[1,1],[0,max(phistd)],'k',90.5*[1,1],[0,max(phistd)],'k',97.5*[1,1],[0,max(phistd)],'k',...
        103.5*[1,1],[0,max(phistd)],'k',111.5*[1,1],[0,max(phistd)],'k')
    legend('mean of phi','standard deviation of phi')
    xlabel('factor')

    print(strcat(dirOut,'plots/',name),'-dpdf')
end