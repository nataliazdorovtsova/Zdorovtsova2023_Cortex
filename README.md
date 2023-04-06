This repository holds analysis scripts for 'Exploring Neural Heterogeneity in Inattention and Hyperactivity' (Zdorovtsova et al., 2023). Each script deals with a different component of the analysis procedure detailed in the paper:

Exploratory factor analysis: '1_EFA.py' (Dependencies: Pandas, FactorAnalyzer, Numpy, and Matplotlib).
Partial least squares regression on nodal graph measure data and EFA factor scores: '2_PLS.m'
Calculation of graph measures and k-means clustering: '3_GraphMeasures_Clustering.m' 
Behavioural and cognitive comparisons (GLMs) between inattentive and hyperactive clusters: '4_BehaviouralAndCognitiveComparisons.m'
Nodal graph measure comparisons (GLMs) between inattentive and hyperactive clusters: '5_NodalMeasureComparisons.m'
(Supplementary materials) K-means clustering on full sample (n = 383): 'Supplement_ClusteringOn383.m'

For some of these analyses, you will need to download a version of the Brain Connectivity Toolbox (link: https://sites.google.com/site/bctnet/)

The MATLAB scripts should be run in order so that variable names are retained and easily accessed. For data availability, please see the statement provided in our published manuscript; as of April 2023, CALM data is being prepared for publication on a public repository. As our data is derived from NHS patients, access will be granted on the basis of an application to the CALM ethics committee.

If you have any further questions, please reach out to the lead author, Natalia Zdorovtsova, at:
natalia.zdorovtsova@mrc-cbu.cam.ac.uk
