%% Comparisons of average communicability, clustering coefficient, and degree between clusters

data = csvread('GraphMeasures232.csv');
    Comm = data(1:100,:);
    ClustCoef = data(101:200,:);
    Degree = data(201:300,:);
    
grouping = csvread('Grouping.csv');

Comm1 = Comm(:,grouping == 1); 
Comm2 = Comm(:,grouping == 2);

ClustCoef1 = ClustCoef(:,grouping == 1);
ClustCoef2 = ClustCoef(:,grouping == 2);

Degree1 = Degree(:,grouping == 1);
Degree2 = Degree(:,grouping == 2);

%% Communicability

Comm1avg = zeros(50,1);
for x = 1:length(Comm1avg)
    Comm1avg(x,:) = mean(Comm1(:,x));
end

Comm2avg = zeros(60,1);
for x = 1:length(Comm2avg)
    Comm2avg(x,:) = mean(Comm2(:,x));
end

[h,p,ci,stats] = ttest2(Comm1avg,Comm2avg)

%% Clustering Coefficient

ClustCoef1avg = zeros(50,1);
for x = 1:length(ClustCoef1avg)
    ClustCoef1avg(x,:) = mean(ClustCoef1(:,x));
end

ClustCoef2avg = zeros(60,1);
for x = 1:length(ClustCoef2avg)
    ClustCoef2avg(x,:) = mean(ClustCoef2(:,x));
end

[h,p,ci,stats] = ttest2(ClustCoef1avg,ClustCoef2avg)

%% Degree

Degree1avg = zeros(50,1);
for x = 1:length(Degree1avg)
    Degree1avg(x,:) = mean(Degree1(:,x));
end

Degree2avg = zeros(60,1);
for x = 1:length(Degree2avg)
    Degree2avg(x,:) = mean(Degree2(:,x));
end

[h,p,ci,stats] = ttest2(Degree1avg,Degree2avg)