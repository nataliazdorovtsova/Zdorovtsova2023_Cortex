%% Run two-block PLS for degree, clustering coef., and communicability across 383 participants

%% Load data

clear; clc;

% set path
cd('/home/nz01/Nodal Analyses/Structural Analyses'); %***** set directory to data folder

% 300 x 383 file for degree, clustering coef, and communicability
graphs = csvread('graphmeasures383.csv');

% specify degree, clust coef, and comm positions in file
degree = graphs(1:100,:);
clustcoef = graphs(101:200,:);
comm = graphs(201:end,:);

% load behavioural data (factor scores)
factorscores = csvread('factorscores_383_new.csv');

% load age data
age = csvread('age383.csv');

%% PLS hyperparameters

% set sample size
samp    = [1:383];


% set hyperparameters *****
nsamp   = length(samp);
nnode   = 100; 
ncomp   = 1; % number of components *****
nperm   = 1000; 
nboot   = 1000;
kfold   = 10;

% set predictors
X       = factorscores;

% set graph metric
Y       = degree;

%%%%%%%%%% can comment this section out if you don't want to age regress %%%%%%%%%%
age_regressed_Y = zeros(size(Y));

for nodeu = 1:nnode
    for nodev = 1:nnode
        node = Y;
        % compute residuals
        b = regress(node, [ones(length(age),1) age]);
        Yhat = [ones(length(age),1) age]*b;
        resid = node - Yhat;
        % keep residuals
        age_regressed_Y(:,nodeu,nodev) = resid;
    end
end

Y = age_regressed_Y;

%% PLS between factor scores and graph measures

% nodal PLS
[Xl1,Yl1,Xs1,Ys1,beta1,var1,mse1] = plsregress(normalize(X,'range'),Y1,ncomp,'CV',kfold);

% A. permutations

% initialise
Xl1_perm  = zeros(nperm,size(X,2),ncomp); 
Yl1_perm  = zeros(nperm,nedge,ncomp);
Xs1_perm  = zeros(nperm,nsamp,ncomp);
Ys1_perm  = zeros(nperm,nsamp,ncomp);
var1_perm = zeros(nperm,2,ncomp);

for perm = 1:nperm
    
    % permute (reorder) the response
    Y1_rand = Y1(randperm(nsamp),:);
    
    [Xl1_perm(perm,:),Yl1_perm(perm,:),...
        Xs1_perm(perm,:),Ys1_perm(perm,:),~,...
        var1_perm(perm,:)] = plsregress(...
        normalize(X,'range'),Y1_rand,ncomp);
    
    disp(sprintf('Permutation %g of %g complete.',perm,nperm));
end

% B. bootstrapping

bootsamp = repmat(Y1,nboot,1); % make an inflated sample
Xl1_boot = zeros(nboot,size(X,2),ncomp);
Yl1_boot = zeros(nboot,nedge,ncomp);

for boot = 1:nboot;
    
    % take a random nsamp sample from bootsamp
    Y1_samp = bootsamp(randperm(nsamp*nboot,nsamp),:);
    
    [Xl1_boot(boot,:),Yl1_boot(boot,:)] = plsregress(...
        normalize(X,'range'),Y1_samp,ncomp);
    
    disp(sprintf('Bootstrap %g of %g complete.',boot,nboot));
end

% then procrusties rotation for bootstrapped loadings
Xl1_boot_Z = [];
Yl1_boot_Z = [];
for boot = 1:nboot
    % predictor
    [~,Xl1_boot_Z(boot,:),~] = procrustes(Xl1,squeeze(Xl1_boot(boot,:)));
    % response
    [~,Yl1_boot_Z(boot,:),~] = procrustes(Yl1,squeeze(Yl1_boot(boot,:)));
    disp(sprintf('Procrusties rotation %g of %g complete.',boot,nboot));
end

%% PLS between environmental and brain measures, edge-wise - ANALYSE

% for each component, variance explained
var1_pred = var1(1,:);
var1_resp = var1(2,:);
cov1      = sum(var1);

% p value in covariance for each component, initialise
var1_null_pred = squeeze(var1_perm(:,1,:));
var1_null_resp = squeeze(var1_perm(:,2,:));
cov1_null      = squeeze(sum(var1_perm,2));

p_pred         = [];
p_resp         = [];
p_cov          = [];

r = []; p_s = [];
for comp = 1:ncomp
    
    % covariance significance
    null_p    = var1_null_pred(:,comp);
    null_r    = var1_null_resp(:,comp);
    null_c    = cov1_null(:,comp);
    
    p_pred(comp) = sum(var1_pred(comp) < null_p) / nperm;
    p_resp(comp) = sum(var1_resp(comp) < null_r) / nperm;
    p_cov(comp)  = sum(cov1(comp) < null_c) / nperm;
    
    % observed score correlations
    [Observed_r(comp) Observed_p_s(comp)] = corr(Xs1(:,comp),Ys1(:,comp));
    
    % permuted score correlations
    for perm = 1:nperm
    [Permuted_r(perm,comp),Permuted_r_s(perm,comp)] = corr(...
        Xs1_perm(perm,:,comp)',Ys1_perm(perm,:,comp)');
    end
end

% calculate score significance against random permutations
p_score = [];
for comp = 1:ncomp
    perm = Permuted_r(:,comp);
    obs  = Observed_r(comp);
    p_score(comp) = sum(obs>perm)/nperm;
end

% tabulate
T1 = table([1:ncomp]',p_pred',p_resp',cov1'*100,p_cov',Observed_r',Observed_p_s',p_score','VariableNames',...
    {'Component','Ppred','Presp',...
    'Covariance explained (%)','Pcov',...
    'Observed score R','Observed score P','Pscore'});

%% PLS between environmental and brain measures, edge-wise - BOOTSTRAPPING CI

% set z for confidence interval
z  = 1.96;

% initialise
CI_Xl1_upper = [];
CI_Xl1_lower = [];
Xl1_diff     = [];

CI_Yl1_upper = [];
CI_Yl1_lower = [];
Yl1_diff     = [];

% the standard deviation of all the bootstrapped values is the standard error of the bootstrap: 
% https://www.youtube.com/watch?v=O_Fj4q8lgmc&t=82s @ 14m

for comp = 1:ncomp
    
    % predictor
    for pred = 1:size(X,2)
    CI_Xl1_upper(pred,comp)     = Xl1(pred,comp) + z * (std(Xl1_boot_Z(:,pred,comp)));
    CI_Xl1_lower(pred,comp)     = Xl1(pred,comp) - z * (std(Xl1_boot_Z(:,pred,comp)));
    Xl1_diff(pred,comp)         = z * (std(Xl1_boot_Z(:,pred,comp)));
    end
    
    % response
    for response = 1:size(Y1,2)
    CI_Yl1_upper(response,comp) = Yl1(response,comp) + z * (std(Yl1_boot_Z(:,response,comp)));
    CI_Yl1_lower(response,comp) = Yl1(response,comp) - z * (std(Yl1_boot_Z(:,response,comp)));
    Yl1_diff(response,comp)     = z * (std(Yl1_boot_Z(:,response,comp)));
    end
    
end

% find predictors which cross 0
notSigPred = [];
for comp = 1:ncomp
    for pred = 1:size(X,2);
        a = CI_Xl1_lower(pred,comp);
        b = CI_Xl1_upper(pred,comp);
        if a < 0 && b > 0;
            notSigPred(pred,comp) = 1;
        end
    end
end

% find edges which cross 0
notSigResp = [];
for comp = 1:ncomp
    for resp = 1:size(Y1,2);
        a = CI_Yl1_lower(resp,comp);
        b = CI_Yl1_upper(resp,comp);
        if a < 0 && b > 0;
            notSigResp(resp,comp) = 1;
        end
    end
end
