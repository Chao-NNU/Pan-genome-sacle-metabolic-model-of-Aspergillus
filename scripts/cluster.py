import os
import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Path settings
input_fba = "/data/gan_1/add_sbml_models/results/FBA.xlsx"
input_meta = "/data/gan_1/add_sbml_models/results/objective.xlsx"
output_dir = "/data/gan_1/add_sbml_models/py/cluster_analysis"
os.makedirs(output_dir, exist_ok=True)

# Step 1: Read data
fba_df = pd.read_excel(input_fba)
meta_df = pd.read_excel(input_meta)

# Extract raw flux data, excluding model name and growth rate
flux_df = fba_df.drop(columns=["model", "BIOMASS_flux"], errors="ignore")
model_names = fba_df["model"]

# Merge metabolite names
flux_df.columns.name = None
col_rename = {row["abbreviaion"]: f"{row['description']} ({row['abbreviaion']})" for _, row in meta_df.iterrows()}
flux_df = flux_df.rename(columns=col_rename)

# Standardization
scaler = StandardScaler()
flux_scaled = scaler.fit_transform(flux_df)
flux_scaled_df = pd.DataFrame(flux_scaled, columns=flux_df.columns)

# Step 2: PCA dimensionality reduction and clustering evaluation
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(flux_scaled)

# Determine the optimal k using inertia and silhouette score
inertia = []
sil_scores = []
K_range = range(5, 11)

for k in K_range:
    km = KMeans(n_clusters=k, random_state=42)
    labels = km.fit_predict(flux_scaled)
    inertia.append(km.inertia_)
    sil_scores.append(silhouette_score(flux_scaled, labels))

# Clustering with the optimal number of clusters
best_k = K_range[np.argmax(sil_scores)]
kmeans = KMeans(n_clusters=best_k, random_state=42)
cluster_labels = kmeans.fit_predict(flux_scaled)

# Save clustering results
fba_df["cluster"] = cluster_labels
fba_df.to_excel(os.path.join(output_dir, "clustered_fba_results.xlsx"), index=False)
