%% Graph metric comparisons between inattentive and hyperactive clusters

IAconnectTrim=IAconnect(:,outliers<1); %%remember to trim out the outliers
Grouping = csvread('kmeans_hidx.csv');

% Retain variables AgeScan, AgeTest, IADensity, and Motion from scripts 3-4
clear NodeTvals NodePvals
for i=1:300;
[b, dev, stats]=glmfit([Grouping, transpose(AgeScan), IAdensity(outliers<1)', transpose(AgeTest), transpose(Motion)],IAconnectTrim(i,:));
NodePvals(i)=stats.p(2);
NodeTvals(i)=stats.t(2);
end

[padj,alpha] = multicmp (NodePvals,'fdr',0.05)
ThreshNodeVals = NodeTvals;
ThreshNodeVals(padj>0.05)=0;
ThreshNodeVals(isnan(ThreshNodeVals))=0;

Comm=ThreshNodeVals([1:100]); %Thresholded communicability t-values
ClustCoef=ThreshNodeVals([101:200]); %Thresholded Clustering t-values
Degree=ThreshNodeVals([201:300]);%Thresholded Degree t-values