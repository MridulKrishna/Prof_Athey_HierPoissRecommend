clear all

tic

% Directory where files will be saved
dir = '/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData';

% Hyperparameters
a = 0.3;
c = 0.3;
e = 0.3;
f = 0.3;
ap = 1.5;
bp = 30;
cp = 1.5;
dp = 30;

% Size of database
U = 12000;
I = 10000;

% Number of factors
K = 4;
L = 3;
M = 2;

lambda = (K+L+M)*a*c*ap*cp/(bp*dp*(ap-1)*(cp-1));
disp(lambda)

% Set random seed
rng(10);

% Generates item and user observables randomly
X = (10+10*rand(I,L));
W = (1+2*rand(U,M));

% Makes each variable have a different scale
X = X*diag(rand(L,1));
W = W*diag(rand(M,1));

% Computes the scale of observables
Xscale = mean(X);
Wscale = mean(W);

% Generates popularity and activity
xi = gamrnd(ap,bp/ap,U,1);
eta = gamrnd(cp,dp/cp,I,1);

% Generates variables for users
theta = zeros(U,K);
sigma = zeros(U,L);

for u = 1:U
    for k = 1:K
        theta(u,k) = gamrnd(a,1/xi(u));
    end
    for l = 1:L
        sigma(u,l) = gamrnd(e,1/(e/(c*a)*Xscale(l)*xi(u)));
    end
end

% Generates variables for items
beta = zeros(I,K);
rho = zeros(I,M);

for i = 1:I
    for k = 1:K
        beta(i,k) = gamrnd(c,1/eta(i));
    end
    for m = 1:M
        rho(i,m) = gamrnd(f,1/(f/(c*a)*Wscale(m)*eta(i)));
    end
end

% Computes observables divided by popularity/activity
XoverEta = X;
for i = 1:I
    for l = 1:L
        XoverEta(i,l) = X(i,l)/eta(i);
    end
end

WoverXi = W;
for u = 1:U
    for m = 1:M
        WoverXi(u,m) = W(u,m)/xi(u);
    end
end

rates = theta*beta'+sigma*XoverEta'+WoverXi*rho';

toc
% mean(mean(rates))

tic
rating = poissrnd(rates);
toc

%max(max(rating));

tic
sprating = sparse(rating);
[iu,ii] = find(sprating);

ratcol = iu;

for i = 1:length(iu)
    ratcol(i) = sprating(iu(i),ii(i));
end

obsUserMat = [10*(1:U)',W];
obsItemMat = [10*(1:I)',X];

ratingMat = [10*iu,10*ii,ratcol];

N = size(ratingMat,1);
perm = randperm(N);
train = perm(1:floor(0.7999999*N));
validation = perm(ceil(0.8*N):floor(0.81*N));
test = perm(ceil(0.81*N):N);
toc

%length(train)+length(validation)+length(test)

tic
dlmwrite(strcat(dir,'/gen',...
    num2str(U),num2str(I),'/obsUser.tsv'),obsUserMat,'\t');
dlmwrite(strcat(dir,'/gen',...
    num2str(U),num2str(I),'/obsItem.tsv'),obsItemMat,'\t');
dlmwrite(strcat(dir,'/gen',num2str(U),num2str(I),'/train.tsv'),...
    ratingMat(train,:),'delimiter','\t','precision','%.0f');

dlmwrite(strcat(dir,'/gen',num2str(U),num2str(I),'/validation.tsv'),...
    ratingMat(validation,:),'delimiter','\t','precision','%.0f');

dlmwrite(strcat(dir,'/gen',num2str(U),num2str(I),'/test.tsv'),...
    ratingMat(test,:),'delimiter','\t','precision','%.0f');
toc
