%% Clustering analysis for inattentive/hyperactive subsample (n = 232)

% This script will consensus-threshold connectome data from participants,
% calculate graph metrics for each participant, and perform various
% different types of clustering.

% Load data
ConnersIDs = csvread('ConnersIDs.csv');
IA=ConnersIDs(ConnersIDs(:,2)==2,1);%select those with the elevated range
Conts=ConnersIDs(ConnersIDs(:,2)~=2,1);%form the control group
IDs383 = csvread('384_CALM_IDs.csv');
load('/imaging/astle/qsiprep_analysis_da04/qsiprep_analysis_processed/calm_qsiprep.mat')%%all the connectomes from CALM
IDs392 = calm_qsiprep.sample.id.calm;
age = calm_qsiprep.sample.age.scan;

%% CONSENSUS THRESHOLDING

clear Y
C = squeeze(calm_qsiprep.schaefer100x7.connectivity(:,1,:,:));
%C = reshape(C,[384,100,100]); %get them into connectome shape
set     = 0.6; % [0 1] ***** this is the thresholding number to change
nsamp   = size(C,1);
t       = floor(nsamp * set);
% find non zero elements
k       = C~=0;
u       = squeeze(sum(k,1));
% keep the indices
ind     = u<=t;
indkeep = u>t;
for sub = 1:nsamp
    A = squeeze(C(sub,:,:));
    % remove edges
    A(ind) = 0;
    Y(sub,:,:) = A;
end

check = [];
for sub = 1:nsamp
    A = squeeze(Y(sub,:,:));
    check(sub) = nnz(A)/2;
end

%% CALCULATE CONNECTOME DENSITY

brainconnectivity = '/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT';
addpath(brainconnectivity);
thr  = 39;
density               = [];
binarised_connectomes = [];
% loop through subjects
W_thr = [];
nsub=size(Y,1);
for sub = 1:nsub;
    W     = squeeze(Y(sub,:,:));          % take unthresholded connectome
    W_thr = threshold_absolute(W,thr);                              % absolute threshold
    density(sub) = density_und(W_thr);                              % calculate density
    B = zeros(size(W_thr));
    B(find(W_thr)) = 1;                                             % binarise
    binarised_connectomes(sub,:,:) = B;
end
% display
disp(sprintf('At a threshold of 39 streamlines, a mean density of %g%% is produced across the sample.',100*mean(density)));


%% CALCULATE GRAPH MEASURES FOR EACH PARTICIPANT

clear Graph
for i = 1:size(Y,1);
%Graph(:,1,i)=betweenness_bin(squeeze(binarised_connectomes(i,:,:)));
Graph(:,1,i)=((mean(expm(squeeze(binarised_connectomes(i,:,:))),1))/2); %mean communicability
Graph(:,2,i)=clustering_coef_bu(squeeze(binarised_connectomes(i,:,:))); %clustering_coef
Graph(:,3,i)=degrees_und(squeeze(binarised_connectomes(i,:,:))); %mean node degree
%Graph(:,4,i)=strengths_und(squeeze(binarised_connectomes(i,:,:)));
end

%% Z SCORING

clear Reshaped_Graph
for i = 1:size(Graph,3);
tmp=reshape(Graph(:,:,i),[300,1]);
Reshaped_Graph(:,i)=tmp;
end

clear Reshaped_Graph_z_tmp
Reshaped_Graph_comm = Reshaped_Graph(1:100,:);
Reshaped_Graph_clust = Reshaped_Graph(101:200,:);
Reshaped_Graph_degree = Reshaped_Graph(201:300,:);

%comm
for i = 1:100
    tmp=(Reshaped_Graph_comm(i,:));
    tmp=zscore(tmp);
    Reshaped_Graph_z_comm(i,:)=tmp;
end

clear tmp
%clust coef
for i = 1:100
    tmp=(Reshaped_Graph_clust(i,:));
    tmp=zscore(tmp);
    Reshaped_Graph_z_clust(i,:)=tmp;
end

clear tmp
%degree
for i = 1:100
    tmp=(Reshaped_Graph_degree(i,:));
    tmp=zscore(tmp);
    Reshaped_Graph_z_degree(i,:)=tmp;
end

% These are the graph metrics for each participant in the sample (n = 383)
Reshaped_Graph_z = cat(1,Reshaped_Graph_z_comm,Reshaped_Graph_z_clust,Reshaped_Graph_z_degree);

%% FIND PARTICIPANTS WITH ELEVATED INATTENTION/HYPERACTIVITY (Conners Questionnaire >60)

clear IAconnect IAids IAdensity
k=1;
for i=1:392;
   if find(IA==IDs392(i))>0, IAconnect(:,k)=Reshaped_Graph_z(:,i);
       IAids(k)=IDs392(i);
       IAdensity(k)=density(i);
       k=k+1;
   end  
end

%% IDENTIFY CONTROLS

clear Controls Controlids Contconnect
k=1;
for i=1:392;
   if find(Conts==IDs392(i))>0, Contconnect(:,k)=Reshaped_Graph_z(:,i);
       Controlids(k)=IDs392(i);
       k=k+1;
   end 
end

%% MDS SCALING

diss = pdist(transpose(IAconnect));
Z = squareform(diss);
figure()
imagesc(Z)
MDS=cmdscale(diss,2); 
maxrelerr = max(abs(diss-pdist(MDS(:,1:2))))./max(diss)
[sig, mu, mah, outliers]=robustcov(MDS); %calculates outliers
MDS=MDS(outliers<1,:); %removes multivariate outliers
eucD = pdist(MDS);

%% HIERARCHICAL CLUSTERING

clustTreeEuc = linkage(eucD,'weighted');
cophenet(clustTreeEuc,eucD)

% calculate silhouettes
clear SilTree
k=1;
for i=1:1:10 % you need to adjust these values based upon your cluster tree
hidx = cluster(clustTreeEuc,'MaxClust',i);
[silh2]= silhouette(MDS,hidx)
SilTree(k)=mean(silh2)
k=1+k;
end

hidx = cluster(clustTreeEuc,'MaxClust',2); % put the value for the best Silhouette in here
figure();
[silh2,h]= silhouette(MDS,hidx)
mean(silh2)

% plot in MDS 2 dimensions
T = cluster(clustTreeEuc,'MaxClust',2);
figure()
gscatter(MDS(:,1),MDS(:,2),T)

%% K-MEANS CLUSTERING

for i=1:10
hidx=kmeans(MDS,i)
[silh2,h]= silhouette(MDS,hidx)
KmeansClusts(i)=mean(silh2) %this gives you the mean silhouette values
end

hidx=kmeans(MDS,2) %%put your preferred cluster membership in here
[silh2,h]= silhouette(MDS,hidx)
mean(silh2) 

Grouping = hidx; % Save this grouping variable for later analyses

%% FUZZY CLUSTERING

[centres,hidx]=fcm(MDS,2)
for i = 1:length(hidx);
    Fuzzyhidx(i)=find(hidx(:,i)==max(hidx(:,i)))
end

[silh2,h]= silhouette(MDS,Fuzzyhidx)
mean(silh2)
