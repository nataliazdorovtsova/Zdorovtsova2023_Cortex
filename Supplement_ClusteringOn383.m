%% Clustering 383 connectomes (full sample)

%% FIND PARTICIPANTS WITH ELEVATED INATTENTION/HYPERACTIVITY

clear IAconnect IAids IAdensity
k=1;
for i=1:392;
   if find(IA==IDs392(i))>0, IAconnect(:,k)=Reshaped_Graph_z(:,i);
       IAids(k)=IDs392(i);
       IAdensity(k)=density(i);
       k=k+1;
   end %%is this subject on the elevated list   
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

fullconnect = cat(2,IAconnect,Contconnect);

%% MDS SCALING

diss = pdist(transpose(fullconnect));
Z = squareform(diss);
figure()
imagesc(Z)
MDS=cmdscale(diss,2); 
maxrelerr = max(abs(diss-pdist(MDS(:,1:2))))./max(diss)
[sig, mu, mah, outliers]=robustcov(MDS); %calculates outliers
MDS=MDS(outliers<1,:); %removes multivariate outliers
eucD = pdist(MDS);


%% KMEANS CLUSTERING
for i=1:10
hidx=kmeans(MDS,i)
[silh2,h]= silhouette(MDS,hidx)
KmeansClusts(i)=mean(silh2) %this gives you the mean Silhouette values
end

hidx=kmeans(MDS,2) %%put your preferred cluster membership in here
[silh2,h]= silhouette(MDS,hidx)
mean(silh2) 
%figure()
%gscatter(MDS(:,1),MDS(:,2),hidx)
Grouping = hidx; 


