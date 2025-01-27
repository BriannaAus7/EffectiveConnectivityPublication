# This code has 2 main purposes. 
# 1. To reduce the high dimensionality of the input dataset, as in out case, it was 94, 45,796. 
# 2. Perform a cluster analysis on the reduced dimensional data, and assess whether the EC patterns are able form clusters that are characteristic of healthy controls and patients.
pwd

file_path ='/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)'

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score
import umap
import hdbscan
import matplotlib.pyplot as plt

def read_and_flatten(file_path):
    df = pd.read_csv(file_path, header=None)
    affinity_matrix = df.iloc[1:, 1:].values.astype(float)
    return affinity_matrix.flatten()

def process_directory(directory, label):
    X, y = [], []
    for filename in os.listdir(directory):
        if filename.endswith("_Aff.csv"):
            filepath = os.path.join(directory, filename)
            X.append(read_and_flatten(filepath))
            y.append(label)
    return X, y

def load_data():
    cue_hc_dir = "/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/CueData/HC_Cue"
    cue_patients_dir = "/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/CueData/Patients_Cue"
    mid_hc_dir = "/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData/HC_MID"
    mid_p_dir = "/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)/MIDData/Patients_MID"

    X, y = [], []
    for directory, label in [(cue_hc_dir, 0), (cue_patients_dir, 1), (mid_hc_dir, 0), (mid_p_dir, 1)]:
        X_dir, y_dir = process_directory(directory, label)
        X.extend(X_dir)
        y.extend(y_dir)
    return np.array(X), np.array(y)


print("Loading data...")
X, y = load_data()
print(f"Data shape: {X.shape}, Labels shape: {y.shape}")


scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)


X_scaled = np.nan_to_num(X_scaled)
#UMAP parameters 
umap_n_neighbors = 10
umap_min_dist = 0.1
umap_n_components = 10

#Applying UMAP
print("Applying UMAP...")
reducer = umap.UMAP(
    n_neighbors=umap_n_neighbors,
    min_dist=umap_min_dist,
    n_components=umap_n_components,
    random_state=42
)
X_reduced = reducer.fit_transform(X_scaled)
print(f"UMAP reduced data shape: {X_reduced.shape}")


#Clustering, specifically chosen HDBSCAN.
hdbscan_min_cluster_size = 5
hdbscan_min_samples = 10

# Apply HDBSCAN
print("Applying HDBSCAN...")
hdbscan_model = hdbscan.HDBSCAN(
    min_samples=hdbscan_min_samples,
    min_cluster_size=hdbscan_min_cluster_size,
    metric='euclidean'
)
clusters = hdbscan_model.fit_predict(X_reduced)

#Cluster eval metrics.
# SS, cluster purity
print("Evaluating clustering...")
non_noise_mask = clusters != -1
silhouette_avg = silhouette_score(X_reduced[non_noise_mask], clusters[non_noise_mask])
print(f"Silhouette Score: {silhouette_avg:.4f}")

# Cluster Purity
def calculate_purity(clusters, labels):
    unique_clusters = np.unique(clusters)
    purity = []
    for cluster_id in unique_clusters:
        if cluster_id == -1:  # Skip noise points
            continue
        cluster_mask = clusters == cluster_id
        cluster_labels = labels[cluster_mask]
        dominant_group = max(np.sum(cluster_labels == 0), np.sum(cluster_labels == 1))
        purity.append(dominant_group / len(cluster_labels))
    return np.mean(purity)

cluster_purity = calculate_purity(clusters, y)
print(f"Cluster Purity: {cluster_purity:.4f}")


#Cluster profiling, also telling me about the noise points.
def cluster_composition_analysis(clusters, labels):
    unique_clusters = np.unique(clusters)
    results = []

    for cluster_id in unique_clusters:
        if cluster_id == -1:  # Skip noise points
            continue
        cluster_mask = clusters == cluster_id
        hc_count = np.sum((labels == 0) & cluster_mask)
        patient_count = np.sum((labels == 1) & cluster_mask)
        total_count = hc_count + patient_count

        results.append({
            'Cluster': cluster_id,
            'HC Count': hc_count,
            'Patient Count': patient_count,
            'Total': total_count,
            'HC %': (hc_count / total_count) * 100 if total_count > 0 else 0,
            'Patient %': (patient_count / total_count) * 100 if total_count > 0 else 0,
        })

    return pd.DataFrame(results)

cluster_composition = cluster_composition_analysis(clusters, y)
print("\nCluster Composition:")
print(cluster_composition)


noise_mask = clusters == -1  #Noise points are labelled as -1, so this way it is only them being identified
num_noise_points = np.sum(noise_mask) 
total_points = len(clusters) 


noise_proportion = (num_noise_points / total_points) * 100

# Patient and HC composition in noise
num_hc_in_noise = np.sum((y == 0) & noise_mask)
num_patients_in_noise = np.sum((y == 1) & noise_mask)


hc_proportion_in_noise = (num_hc_in_noise / num_noise_points) * 100 if num_noise_points > 0 else 0
patient_proportion_in_noise = (num_patients_in_noise / num_noise_points) * 100 if num_noise_points > 0 else 0


print("\nNoise Point Analysis:")
print(f"Total Noise Points: {num_noise_points} ({noise_proportion:.2f}%)")
print(f"HC in Noise: {num_hc_in_noise} ({hc_proportion_in_noise:.2f}%)")
print(f"Patients in Noise: {num_patients_in_noise} ({patient_proportion_in_noise:.2f}%)")


# Visualize clusters including noise

plt.figure(figsize=(10, 8))
plt.scatter(
    X_reduced[clusters != -1, 0],  # UMAP dimension 1 for clustered points
    X_reduced[clusters != -1, 1],  # UMAP dimension 2 for clustered points
    c=clusters[clusters != -1],  # Use cluster labels for color
    cmap='viridis',  # Colormap for clusters
    label='Clusters',
    alpha=0.7,  # Transparency
    s=30  # Marker size
)
plt.scatter(
    X_reduced[clusters == -1, 0],  # UMAP dimension 1 for noise points
    X_reduced[clusters == -1, 1],  # UMAP dimension 2 for noise points
    c='red',  # Color for noise points
    label='Noise',
    alpha=0.5,  # Transparency for noise points
    s=30  # Marker size
)
plt.title("UMAP Embedding with Clusters and Noise Points", fontsize=16)
plt.xlabel("UMAP Dimension 1", fontsize=12)
plt.ylabel("UMAP Dimension 2", fontsize=12)
plt.legend(loc="best", fontsize=10)
plt.grid(True, alpha=0.3)
plt.show()

from sklearn.manifold import trustworthiness

# Calculate trustworthiness
trust = trustworthiness(X_scaled, X_reduced, n_neighbors=umap_n_neighbors)
print(f"Trustworthiness Score: {trust:.4f}")

from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting

# 3D Visualization of UMAP embedding with clusters and noise points
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot clustered points
scatter_clusters = ax.scatter(
    X_reduced[clusters != -1, 0],  # UMAP dimension 1
    X_reduced[clusters != -1, 1],  # UMAP dimension 2
    X_reduced[clusters != -1, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[clusters != -1].shape[0]),  # UMAP dimension 3 (or zeros if only 2D)
    c=clusters[clusters != -1],  # Color based on cluster ID
    cmap='viridis',
    label='Clusters',
    alpha=0.7,
    s=30
)

# Plot noise points
scatter_noise = ax.scatter(
    X_reduced[clusters == -1, 0],  # UMAP dimension 1
    X_reduced[clusters == -1, 1],  # UMAP dimension 2
    X_reduced[clusters == -1, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[clusters == -1].shape[0]),  # UMAP dimension 3 (or zeros)
    c='red',
    label='Noise',
    alpha=0.5,
    s=30
)

# Add labels and title
ax.set_title("3D UMAP Embedding with Clusters and Noise Points", fontsize=16)
ax.set_xlabel("UMAP Dimension 1", fontsize=12)
ax.set_ylabel("UMAP Dimension 2", fontsize=12)
ax.set_zlabel("UMAP Dimension 3", fontsize=12)
plt.legend(loc="best", fontsize=10)
plt.show()

# Visualization by HC and Patient labels
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# HC points
scatter_hc = ax.scatter(
    X_reduced[y == 0, 0],  # UMAP dimension 1 for HC
    X_reduced[y == 0, 1],  # UMAP dimension 2 for HC
    X_reduced[y == 0, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[y == 0].shape[0]),  # UMAP dimension 3 (or zeros)
    c='blue',
    label='HC',
    alpha=0.7,
    s=30
)

# Patient points
scatter_patients = ax.scatter(
    X_reduced[y == 1, 0],  # UMAP dimension 1 for Patients
    X_reduced[y == 1, 1],  # UMAP dimension 2 for Patients
    X_reduced[y == 1, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[y == 1].shape[0]),  # UMAP dimension 3 (or zeros)
    c='orange',
    label='Patients',
    alpha=0.7,
    s=30
)

# Add labels and title
ax.set_title("3D UMAP Embedding with HC and Patient Labels", fontsize=16)
ax.set_xlabel("UMAP Dimension 1", fontsize=12)
ax.set_ylabel("UMAP Dimension 2", fontsize=12)
ax.set_zlabel("UMAP Dimension 3", fontsize=12)
plt.legend(loc="best", fontsize=10)
plt.show()

# 1. Feature Importance
#2 ways to assess feature Importance.

#1. Within Each Cluster (Grouping Clusters)
#Which features most influence the grouping of points in each cluster. Specifcially chose the Kruskal-Wallis as non-parametric data.

#2. Distinguishing Clusters (Seperating Clusters)

def calculate_feature_importance(X, clusters):
    """
    Calculate feature importance using Kruskal-Wallis test for non-parametric data.
    Skips features with identical values within clusters.
    """
    feature_importance = []
    unique_clusters = np.unique(clusters)

    for feature_idx in range(X.shape[1]):
        feature_values = X[:, feature_idx]

        # Group feature values by cluster
        cluster_values = [feature_values[clusters == cluster] for cluster in unique_clusters if cluster != -1]

        # Check if all values are identical across all clusters
        if all(len(values) > 1 and np.any(np.diff(values)) for values in cluster_values):
            # Perform Kruskal-Wallis test
            stat, p_value = kruskal(*cluster_values)
            feature_importance.append((feature_idx, stat, p_value))

    # Convert to DataFrame
    importance_df = pd.DataFrame(feature_importance, columns=['Feature', 'K-Statistic', 'P-Value'])
    return importance_df.sort_values(by='K-Statistic', ascending=False)


pip install tqdm
from tqdm import tqdm


# Permutation Testing functions

def permutation_test(X, labels, feature_index, num_permutations=1000):
    group_0 = X[labels == 0, feature_index]
    group_1 = X[labels == 1, feature_index]
    observed_stat = np.mean(group_1) - np.mean(group_0)

    null_dist = []
    for _ in range(num_permutations):
        permuted_labels = np.random.permutation(labels)
        perm_group_0 = X[permuted_labels == 0, feature_index]
        perm_group_1 = X[permuted_labels == 1, feature_index]
        null_stat = np.mean(perm_group_1) - np.mean(perm_group_0)
        null_dist.append(null_stat)

    p_value = np.mean(np.abs(null_dist) >= np.abs(observed_stat))
    return observed_stat, p_value

def run_permutation_tests(X, labels, top_features, num_permutations=1000):
    results = []
    for feature in tqdm(top_features, desc="Running Permutation Tests"):
        obs_stat, p_val = permutation_test(X, labels, feature, num_permutations)
        results.append({'Feature': feature, 'Observed Statistic': obs_stat, 'P-Value': p_val})
    return pd.DataFrame(results)
