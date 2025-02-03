import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import zscore

# Load the dataset
data = pd.read_csv("static/examples/matrix_myeloma.csv")

# Transpose for PCA (samples as rows, genes as columns)
data_T = data.set_index("ID_REF").transpose()

# Standardize the data
scaler = StandardScaler()
data_standardized = scaler.fit_transform(data_T)

# Perform PCA
pca = PCA(n_components=2)  # Start with 2 components
principal_components = pca.fit_transform(data_standardized)

# Compute Z-scores for the first two principal components
pc_df = pd.DataFrame(principal_components, columns=['PC1', 'PC2'])
pc_df['Z_PC1'] = zscore(pc_df['PC1'])
pc_df['Z_PC2'] = zscore(pc_df['PC2'])

# Identify outliers using a Z-score threshold (e.g., >3 standard deviations)
threshold = 1.5
outliers = pc_df[(abs(pc_df['Z_PC1']) > threshold) | (abs(pc_df['Z_PC2']) > threshold)]

# Get GSM identifiers of outliers
outlier_gsms = data_T.index[outliers.index]

print("Outliers:", list(outlier_gsms))

# Remove outliers
data_filtered = data.drop(columns=outlier_gsms)

# Save or proceed with filtered data
output_file = "filtered_data.csv"
data_filtered.to_csv(output_file, index=False)
