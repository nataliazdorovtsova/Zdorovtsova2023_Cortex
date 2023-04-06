%% Statistical comparisons on behavioural and cognitive data between the two inattentive/hyperactive clusters

% Run this script directly after the third script in the pipeline (graph
% metric calculations and clustering) in order to preserve variable names.
% This script will compare the two clusters on twenty measures of behaviour
% and seven measures of cognitive ability.

%% Extract data for inattentive/hyperactive clusters

clear Pheno Motion AgeTest AgeScan
tmp=IAids(outliers<1) %%remember to take out the MDS outliers
test=xlsread('CALM_Data_Full_Standardised.xlsx');
for i=1:length(tmp);
Pheno(i,:)=test(find(test(:,1)==tmp(i)),:); %%this should now be in the right order
Motion(i)=calm_qsiprep.qc.summary(find(calm_qsiprep.sample.id.calm==tmp(i)),1);
if str2num(cell2mat(calm_qsiprep.sample.age.test(find(calm_qsiprep.sample.id.calm==tmp(i)))))>1,
AgeTest(i)=str2num(cell2mat(calm_qsiprep.sample.age.test(find(calm_qsiprep.sample.id.calm==tmp(i)))));
end
AgeScan(i)=calm_qsiprep.sample.age.scan(find(calm_qsiprep.sample.id.calm==tmp(i)));
end

[transpose(AgeTest*12),Pheno(:,2)] %this is a sanity check... they should be identical if we have successfully lined up the pheno and connectomes

%% Extract data for 'controls'

%This creates similar fields for the controls (for later comparisons, if
%desired)
clear ContPheno ContMotion ContAgeTest ContAgeScan
tmp=Controlids; %%remember to take out the outliers
test=xlsread('CALM_Data_Full_Standardised.xlsx');
for i=1:length(tmp);
ContPheno(i,:)=test(find(test(:,1)==tmp(i)),:); %%this should now be in the right order
ContMotion(i)=calm_qsiprep.qc.summary(find(calm_qsiprep.sample.id.calm==tmp(i)),1);
if str2num(cell2mat(calm_qsiprep.sample.age.test(find(calm_qsiprep.sample.id.calm==tmp(i)))))>1,
ContAgeTest(i)=str2num(cell2mat(calm_qsiprep.sample.age.test(find(calm_qsiprep.sample.id.calm==tmp(i)))));
end
ContAgeScan(i)=calm_qsiprep.sample.age.scan(find(calm_qsiprep.sample.id.calm==tmp(i)));
end

[transpose(ContAgeTest*12),ContPheno(:,2)]
Gender = csvread('Gender BOTH CLUSTERS.csv');

%% GLMs

% Here, hidx is a grouping variable which corresponds to the groups found
% during the k-means clustering procedure. 

clear stats.p
for i= 1:33 %cognition = [4 5 8 10 11 12 13] and behaviour = 14:33
[b, dev, stats]=glmfit([hidx(ismember(hidx,[1,2])),transpose(AgeScan(ismember(hidx,[1,2]))), transpose(AgeTest(ismember(hidx,[1,2]))),transpose(Motion(ismember(hidx,[1,2])))],Pheno((ismember(hidx,[1,2])),i));
Pvals(i)=stats.p(2);
Tvals(i)=stats.t(2);
end

cognition=Pvals([4 5 8 10 11 12 13])
[padj,alpha] = multicmp(cognition,'fdr',0.05)

behaviour=Pvals([14:33])
[padj,alpha] = multicmp(behaviour,'fdr',0.05)
