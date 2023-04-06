#import libraries as required
import pandas as pd
from factor_analyzer import FactorAnalyzer
from factor_analyzer import utils
import matplotlib.pyplot as plt
import numpy as np

#import the dataset
df = pd.read_csv("BehaviouralData_struct_clean_labels.csv")
df.info()
df.head()

#Normalise the data
from scipy.stats import zscore
df.apply(zscore)

# Create factor analysis object and perform factor analysis
# Try: with and without BRIEF: one factor, no rotation vs oblimin, two factors, no rotation vs varimax. 
fa = FactorAnalyzer(n_factors=1,rotation=None,method='minres',
               use_smc=True,is_corr_matrix=False,bounds=(0.005,1),impute='median',
               svd_method='randomized', rotation_kwargs=None)
fa.fit(df)
loadings = fa.loadings_ # factor loadings
communalities = fa.get_communalities()
variance = fa.get_factor_variance() # variance, prop. variance**, and cumulative variance for each factor
uniquenesses = fa.get_uniquenesses()
factorscores = fa.transform(df) #factor scores for each participant

# Check Eigenvalues
ev, v = fa.get_eigenvalues()
print(ev)

# Create scree plot using matplotlib
plt.scatter(range(1,df.shape[1]+1),ev)
plt.plot(range(1,df.shape[1]+1),ev)
plt.title('Scree Plot')
plt.xlabel('Factors')
plt.ylabel('Eigenvalue')
plt.grid()
plt.show()

#Bartlett's Test of Sphericity
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
chi_square_value,p_value=calculate_bartlett_sphericity(df)
print()
print(chi_square_value, p_value)

#Kaiser-Meyer-Olkin (KMO) Test
from factor_analyzer import calculate_kmo
kmo_all,kmo_model=calculate_kmo(df)
print(kmo_model)

# Correlation matrix between questionnaires
matrix = df.corr()



