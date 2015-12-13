clear
clc

% Loads the test files
train = dlmread('/afs/.ir/users/m/r/mridulk/GSB_RA/RecommendSysHPF/PoissFactObs/Yogurt/observables/train.tsv','\t');
test = dlmread('/afs/.ir/users/m/r/mridulk/GSB_RA/RecommendSysHPF/PoissFactObs/Yogurt/observables/test.tsv','\t');
validation = dlmread('/afs/.ir/users/m/r/mridulk/GSB_RA/RecommendSysHPF/PoissFactObs/Yogurt/observables/validation.tsv','\t');
obsItem = dlmread('/afs/.ir/users/m/r/mridulk/GSB_RA/RecommendSysHPF/PoissFactObs/Yogurt/observables/obsItem.tsv','\t');
obsUser = dlmread('/afs/.ir/users/m/r/mridulk/GSB_RA/RecommendSysHPF/PoissFactObs/Yogurt/observables/obsUser.tsv','\t');

% % Loads the files
% train = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/train.tsv','\t');
% test = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/test.tsv','\t');
% validation = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/validation.tsv','\t');
% obsItem = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/obsItem.tsv','\t');
% obsUser = dlmread('~/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt/observables/obsUser.tsv','\t');

% obsItemB = obsItem(:,2:end);
% obsItemID = num2str(obsItem(:,1));

% Keys from item/user code to the order in which they will be stored in
% matrices
% Item/user code to number or number to code
i_c2n = containers.Map('KeyType','uint64','ValueType','uint64');
i_n2c = containers.Map('KeyType','uint64','ValueType','uint64');
u_c2n = containers.Map('KeyType','uint64','ValueType','uint64');
u_n2c = containers.Map('KeyType','uint64','ValueType','uint64');

% Reads files and saves the keys
i = 1;
u = 1;
for j = 1:size(train,1)
    if ~isKey(u_c2n,train(j,1)) 
        u_c2n(train(j,1))=u;
        u_n2c(u)=train(j,1);
        u = u+1;
    end
    
    if ~isKey(i_c2n,train(j,2)) 
        i_c2n(train(j,2))=i;
        i_n2c(i)=train(j,2);
        i = i+1;
    end
end
for j = 1:size(validation,1)
    if ~isKey(u_c2n,validation(j,1)) 
        u_c2n(validation(j,1))=u;
        u_n2c(u)=validation(j,1);
        u = u+1;
    end
    
    if ~isKey(i_c2n,validation(j,2)) 
        i_c2n(validation(j,2))=i;
        i_n2c(i)=validation(j,2);
        i = i+1;
    end
end
for j = 1:size(test,1)
    if ~isKey(u_c2n,test(j,1)) 
        u_c2n(test(j,1))=u;
        u_n2c(u)=test(j,1);
        u = u+1;
    end
    
    if ~isKey(i_c2n,test(j,2)) 
        i_c2n(test(j,2))=i;
        i_n2c(i)=test(j,2);
        i = i+1;
    end
end

% Number of users and items
U = u-1;
I = i-1;

% Saves train as a sparse array
trainS = sparse(U,I);
trainM = zeros(U,I);

for j = 1:size(train,1)
    trainS(u_c2n(train(j,1)),i_c2n(train(j,2))) = train(j,3);
    trainM(u_c2n(train(j,1)),i_c2n(train(j,2))) = train(j,3);
end

%trainS(u_c2n(30),i_c2n(70))

% Number of observed and unobserved characteristics
K = 7;
L = size(obsItem,2)-1;  % Number of observed item characteristics
M = size(obsUser,2)-1;  % Number of observed user characteristics
% L=0;
% M=0;

% Saves item and user observables as arrays
obsU = zeros(U,M);
obsI = zeros(I,L);

if L ~= 0
    for j = 1:size(obsUser,1)
        if obsU(u_c2n(obsUser(j,1)),1) ~= 0
            disp('Duplicate users')
        end

        obsU(u_c2n(obsUser(j,1)),:) = obsUser(j,2:end);
    end
end

if M ~= 0
    for j = 1:size(obsItem,1)
        if obsI(i_c2n(obsItem(j,1)),1) ~= 0
            disp('Duplicate items')
        end

        obsI(i_c2n(obsItem(j,1)),:) = obsItem(j,2:end);
    end
end

% 1: ones
% 2: mean
% 3: standard deviation
scale = 2;
scaleFactor = 1;

if scale == 1
    scaleU = scaleFactor*ones(1,M);
    scaleI = scaleFactor*ones(1,L);
elseif scale == 2
    scaleU = scaleFactor*mean(obsU,1);
    scaleI = scaleFactor*mean(obsI,1);
elseif scale == 3
    scaleU = scaleFactor*std(obsU,0,1);
    scaleI = scaleFactor*std(obsI,0,1);
end

% Hyperparameters
a = 0.3;
ap = 0.4;
bp = 0.5;
c = 0.1;
cp = 0.6;
dp = 0.7;
e = 0.8;
f = 0.9;

totRanking = sum(train(:,3));

bp = sqrt(ap*cp*a*c/(ap-1)/(cp-1)*(K+L+M)*U*I/totRanking);
dp = bp;

offset = 0.00;

% Initializes parameters
gamma_s = a*ones(U,K) + offset*0.01*randn(U,K);
gamma_r = bp*ones(U,K) + offset*0.01*randn(U,K);
mu_s = e*ones(U,L) + offset*0.01*randn(U,L);
% mu_r = bp*ones(U,L) + offset*0.01*randn(U,L);
mu_r = e*bp/(c*a)*ones(U,1)*scaleI + offset*0.01*randn(U,L);
kappa_s = (ap)*ones(U,1);
% kappa_s = (ap+K*a+L*e)*ones(U,1);
kappa_r = ap/bp*ones(U,1) + offset*0.01*randn(U,1);

lambda_s = c*ones(I,K) + offset*0.01*randn(I,K);
lambda_r = dp*ones(I,K) + offset*0.01*randn(I,K);
nu_s = f*ones(I,M) + offset*0.01*randn(I,M);
% nu_r = dp*ones(I,M) + offset*0.01*randn(I,M);
nu_r = f*dp/(c*a)*ones(I,1)*scaleU + offset*0.01*randn(I,M);
tau_s = (cp)*ones(I,1);
% tau_s = (cp+K*c+M*f)*ones(I,1);
tau_r = cp/dp*ones(I,1) + offset*0.01*randn(I,1);

phi = zeros(U,I,K+L+M);
yphi = zeros(U,I,K+L+M);
yphik = zeros(U,I,K);
yphil = zeros(U,I,L);
yphim = zeros(U,I,M);

% Main loop
r=1;
while (r <= 60)
    % Update multinomial parameters
    for u = 1:U
        for i = 1:I
            for k = 1:K
                phi(u,i,k) = exp( psi(gamma_s(u,k))-log(gamma_r(u,k)) + ...
                    psi(lambda_s(i,k))-log(lambda_r(i,k)) );
                
%                 disp(psi(gamma_s(u,k))-log(gamma_r(u,k))+ psi(lambda_s(i,k))-log(lambda_r(i,k)))
            end
            
            for l = 1:L
                phi(u,i,K+l) = obsI(i,l)*exp( psi(mu_s(u,l))-log(mu_r(u,l))-psi(tau_s(i))+log(tau_r(i)) );
                
%                 disp(log(obsI(i,l)*exp( psi(mu_s(u,l))-log(mu_r(u,l))-psi(tau_s(i))+log(tau_r(i)) )));
            end
            
            for m = 1:M
                phi(u,i,K+L+m) = obsU(u,m)*exp( psi(nu_s(i,m))-log(nu_r(i,m))-psi(kappa_s(u))+log(kappa_r(u)) );
%                 disp( log( obsU(u,m)*exp( psi(nu_s(i,m))-log(nu_r(i,m))-psi(kappa_s(u))+log(kappa_r(u)) ) ) )
%                 disp( log( obsU(u,m) ));
%                 disp( psi(nu_s(i,m))-log(nu_r(i,m)) );
%                 disp( psi(kappa_s(u))-log(kappa_r(u)) );
%                 disp( psi(tau_s(u))-log(tau_r(u)) );
            end
            
%             squeeze(log(phi(u,i,:)))'
            
            phi(u,i,:) = phi(u,i,:)/(sum(phi(u,i,:)));
            
            yphi(u,i,:) = trainM(u,i)*phi(u,i,:);
           
%             squeeze(log(yphi(u,i,:)))'
        end
    end
    
    yphik = yphi(:,:,1:K);
    yphil = yphi(:,:,K+1:K+L);
    yphim = yphi(:,:,K+L+1:K+L+M);

    % First cycle: updates user parameters
    Ngamma_s = a*ones(U,K);
    Nmu_s = e*ones(U,L);

    Ngamma_r = zeros(U,K);
    Nmu_r = zeros(U,L);
    for u = 1:U
        Ngamma_r(u,:) = kappa_s(u)/kappa_r(u);
        Nmu_r(u,:) = e/(c*a)*scaleI*kappa_s(u)/kappa_r(u);
    end
    
    for u = 1:U
        for i = 1:I
            for k = 1:K
                Ngamma_s(u,k) = Ngamma_s(u,k) + yphik(u,i,k);
                Ngamma_r(u,k) = Ngamma_r(u,k) + lambda_s(i,k)/lambda_r(i,k);
            end
            
            for l = 1:L 
                Nmu_s(u,l) = Nmu_s(u,l) + yphil(u,i,l);
                Nmu_r(u,l) = Nmu_r(u,l) + tau_r(i)/(tau_s(i)-1)*obsI(i,l);
            end
        end
    end
   
    % Second cycle: updates item parameters
    
    Nlambda_s = c*ones(I,K);
    Nnu_s = f*ones(I,M);
    
    Nlambda_r = zeros(I,K);
    Nnu_r = zeros(I,M);
    for i = 1:I
        Nlambda_r(i,:) = tau_s(i)/tau_r(i);
        Nnu_r(i,:) = f/(c*a)*scaleU*tau_s(i)/tau_r(i);
    end

    for u = 1:U
        for i = 1:I
            for k = 1:K
                Nlambda_s(i,k) = Nlambda_s(i,k) + yphik(u,i,k);
                Nlambda_r(i,k) = Nlambda_r(i,k) + Ngamma_s(u,k)/Ngamma_r(u,k);
            end

            for m = 1:M 
                Nnu_s(i,m) = Nnu_s(i,m) + yphim(u,i,m);
                Nnu_r(i,m) = Nnu_r(i,m) + kappa_r(u)/(kappa_s(u)-1)*obsU(u,m);
            end
        end
    end
    
    Nlambda_s;
           
    % Final cycle: updates popularity and activity parameters
    
    Nkappa_s = (ap+K*a+L*e)*ones(U,1);
    Nkappa_r = ap/bp*ones(U,1);
    Ntau_s = (cp+K*c+M*f)*ones(I,1);
    Ntau_r = cp/dp*ones(I,1);
    
    colsum = zeros(u,1);
    colsum2 = zeros(u,1);
    
    for u = 1:U
        for k = 1:K
            Nkappa_r(u) = Nkappa_r(u) + Ngamma_s(u,k)/Ngamma_r(u,k);
%             colsum(u) = colsum(u) + Ngamma_s(u,k)/Ngamma_r(u,k);
        end
        
        for l = 1:L
            Nkappa_r(u) = Nkappa_r(u) + e/(c*a)*scaleI(l)*Nmu_s(u,l)/Nmu_r(u,l);
            %colsum2(u) = colsum2(u) + e/(c*a)*scaleI(l)*Nmu_s(u,l)/Nmu_r(u,l);
%             colsum2(u) = colsum2(u) + scaleI(l)*Nmu_s(u,l)/Nmu_r(u,l);
        end
    end
    
    for i = 1:I
        for k = 1:K
            Ntau_r(i) = Ntau_r(i) + Nlambda_s(i,k)/Nlambda_r(i,k);
        end
        
        for m = 1:M
            Ntau_r(i) = Ntau_r(i) + f/(c*a)*scaleU(m)*Nnu_s(i,m)/Nnu_r(i,m);
        end
    end
    
    % Save new parameters as current
    gamma_s = Ngamma_s;
    gamma_r = Ngamma_r;
    lambda_s = Nlambda_s;
    lambda_r = Nlambda_r;
    mu_s = Nmu_s;
    mu_r = Nmu_r;
    nu_s = Nnu_s;
    nu_r = Nnu_r;
    kappa_r = Nkappa_r;
    kappa_s = Nkappa_s;
    tau_r = Ntau_r;
    tau_s = Ntau_s;
    
    % Compute log likelihood of validation set
    
    llValidation = 0;
    
    for i = 1:size(validation,1)
        theta = gamma_s(u_c2n(validation(i,1)),:)./gamma_r(u_c2n(validation(i,1)),:);
        beta = lambda_s(i_c2n(validation(i,2)),:)./lambda_r(i_c2n(validation(i,2)),:);
        sigma = mu_s(u_c2n(validation(i,1)),:)./mu_r(u_c2n(validation(i,1)),:);
        rho = nu_s(i_c2n(validation(i,2)),:)./nu_r(i_c2n(validation(i,2)),:);
        x = obsI(i_c2n(validation(i,2)),:);
        w = obsU(u_c2n(validation(i,1)),:);
        
        rate = theta*beta'+sigma*x'+rho*w';
        
        rating = validation(i,3);
        
        indLL = rating*log(rate)-rate-log(factorial(rating));
        
        llValidation = llValidation + indLL;
    end
    
%     disp(llValidation)
    
    % Compute log likelihood of test set
    
    llTest = 0;
    
    for i = 1:size(test,1)
        theta = gamma_s(u_c2n(test(i,1)),:)./gamma_r(u_c2n(test(i,1)),:);
        beta = lambda_s(i_c2n(test(i,2)),:)./lambda_r(i_c2n(test(i,2)),:);
        sigma = mu_s(u_c2n(test(i,1)),:)./mu_r(u_c2n(test(i,1)),:);
        rho = nu_s(i_c2n(test(i,2)),:)./nu_r(i_c2n(test(i,2)),:);
        x = obsI(i_c2n(test(i,2)),:);
        w = obsU(u_c2n(test(i,1)),:);
        
        rate = theta*beta'+sigma*x'+rho*w';
        
        rating = test(i,3);
        
        indLL = rating*log(rate)-rate-log(factorial(rating));
        
        llTest = llTest + indLL;
    end
    
     disp(llTest)
    
    r = r+1;
end
