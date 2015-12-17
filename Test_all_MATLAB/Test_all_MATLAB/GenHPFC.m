% R = gamrnd(shape,1/rate,dimensions);

a = 0.3;
c = 0.3;
e = 0.3;
f = 0.3;
ap = 0.3;
bp = 3;
cp = 0.3;
dp = 3;

U = 10000;
I = 10000;
K = 4;
L = 3;
M = 2;

rng(10);

X = 1/10*(10+10*rand(I,L));
W = 1/10*(5+20*rand(U,M));

xi = gamrnd(ap,ap/bp,U,1);
eta = gamrnd(cp,cp/dp,I,1);

theta = zeros(U,K);
sigma = zeros(U,L);

for u = 1:U
    for k = 1:K
        theta(u,k) = gamrnd(a,xi(u));
    end
    for l = 1:L
        sigma(u,l) = gamrnd(e,xi(u));
    end
end

beta = zeros(I,K);
rho = zeros(I,M);

for i = 1:I
    for k = 1:K
        beta(i,k) = gamrnd(c,eta(i));
    end
    for m = 1:M
        rho(i,m) = gamrnd(f,eta(i));
    end
end

rates = theta*beta'+sigma*X'+W*rho';

%mean(mean(rates));

tic
rating = poissrnd(rates);
toc

%max(max(rating));

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

%length(train)+length(validation)+length(test)

tic
dlmwrite(strcat('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData/gen',...
    num2str(U),num2str(I),'/obsUser.tsv'),obsUserMat,'\t');
dlmwrite(strcat('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData/gen',...
    num2str(U),num2str(I),'/obsItem.tsv'),obsItemMat,'\t');
dlmwrite(strcat('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData/gen',...
    num2str(U),num2str(I),'/train.tsv'),ratingMat(train,:),'\t');
dlmwrite(strcat('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData/gen',...
    num2str(U),num2str(I),'/validation.tsv'),ratingMat(validation,:),'\t');
dlmwrite(strcat('/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/HPF/generatedData/gen',...
    num2str(U),num2str(I),'/test.tsv'),ratingMat(test,:),'\t');
toc
