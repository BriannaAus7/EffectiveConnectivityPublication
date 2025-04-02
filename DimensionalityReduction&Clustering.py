# This code has 2 main purposes.
# 1. To reduce the high dimensionality of the input dataset, as in our case, it was 94, 45,796. This is derived from the pairwise EC matrices for each individual.
# 2. Perform a cluster analysis on the reduced dimensional data, and assess whether the EC patterns are able form clusters that are characteristic of healthy controls and patients.
# 3. Assess the noise cluster, and investigate whether the noise points are boundaried or transitional compared to the main clusters formed.

file_path ='/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)'

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score
import umap
import hdbscan
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu

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

#UMAP being applied
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

#Cluster evaluation metrics (Silhouette Score and Cluster Purity).
# Silhouette Score
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

# MD patient and HC composition in noise
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
    X_reduced[clusters != -1, 0],  # UMAP 1-D for clustered points
    X_reduced[clusters != -1, 1],  # UMAP 2-D for clustered points
    c=clusters[clusters != -1],  # Use cluster labels for color
    cmap='viridis',  # Colourmap for clusters
    label='Clusters',
    alpha=0.7,  # Transparency
    s=30  # Marker size
)
plt.scatter(
    X_reduced[clusters == -1, 0],  # UMAP 1-D for noise points
    X_reduced[clusters == -1, 1],  # UMAP 2-D for noise points
    c='red',  # Colour for noise points
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

# Calculate trustworthiness. (Scroll further down in the code for explanation of trustworthiness.)
trust = trustworthiness(X_scaled, X_reduced, n_neighbors=umap_n_neighbors)
print(f"Trustworthiness Score: {trust:.4f}")

from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting as want to visualise the UMAP embedding in 3D

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

ax.set_title("3D UMAP Embedding with Clusters and Noise Points", fontsize=16)
ax.set_xlabel("UMAP Dimension 1", fontsize=12)
ax.set_ylabel("UMAP Dimension 2", fontsize=12)
ax.set_zlabel("UMAP Dimension 3", fontsize=12)
plt.legend(loc="best", fontsize=10)
plt.show()

# Visualization with HC and Patient labels now, to see the distribution
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# HC points
scatter_hc = ax.scatter(
    X_reduced[y == 0, 0],  # UMAP 1-D for HC
    X_reduced[y == 0, 1],  # UMAP 2-D for HC
    X_reduced[y == 0, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[y == 0].shape[0]),  # UMAP 3-D (or zeros)
    c='blue',
    label='HC',
    alpha=0.7,
    s=30
)

# Patient (MD) points
scatter_patients = ax.scatter(
    X_reduced[y == 1, 0],  # UMAP 1-D for Patients
    X_reduced[y == 1, 1],  # UMAP 2-D for Patients
    X_reduced[y == 1, 2] if X_reduced.shape[1] > 2 else np.zeros(X_reduced[y == 1].shape[0]),  # UMAP 3-D (or zeros)
    c='orange',
    label='Patients',
    alpha=0.7,
    s=30
)


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
    Skips features with identical values within clusters as Kruskal cannot be done with identical values.
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

#Plotting functions

def plot_feature_distribution(X, labels, feature_idx, group_labels):
    control_values = X[labels == 0, feature_idx]
    patient_values = X[labels == 1, feature_idx]
    plt.boxplot([control_values, patient_values], labels=group_labels)
    plt.title(f"Feature {feature_idx} Distribution")
    plt.ylabel("Value")
    plt.show()

def map_features_to_region_pairs(feature_indices, region_names):
    num_regions = len(region_names)
    mapped_features = []

    for feature_idx in feature_indices:
        region1 = feature_idx // num_regions
        region2 = feature_idx % num_regions
        region_pair = (region_names[region1], region_names[region2])
        mapped_features.append(region_pair)

    return mapped_features


    # Feature Importance
    print("Calculating feature importance...")
    feature_importance_df = calculate_feature_importance(X_scaled, clusters)
    print("\nTop 10 Important Features:")
    print(feature_importance_df.head(10))

    # Permutation Testing
    top_features = feature_importance_df['Feature'].values[:10]
    print("\nRunning permutation tests for top features...")
    permutation_results_df = run_permutation_tests(X_scaled, y, top_features)
    print("\nPermutation Test Results:")
    print(permutation_results_df)

    # Plot the significant figures
    significant_features = permutation_results_df[permutation_results_df['P-Value'] < 0.05]['Feature'].values
    for feature in significant_features:
        plot_feature_distribution(X_scaled, y, feature, ["Control", "Patient"])

    # Map the features to the brain regions, so can see which networks the coming up.
    region_names_df = pd.read_excel("ComboNames.xlsx", header=None)
    region_names = region_names_df[0].tolist()
    mapped_features = map_features_to_region_pairs(significant_features, region_names)

    print("\nMapped Features:")
    for feature, region_pair in zip(significant_features, mapped_features):
        print(f"Feature {feature}: {region_pair}")


kruskal_significance_level = 0.05

significant_kruskal_features = feature_importance_df[
    feature_importance_df['P-Value'] < kruskal_significance_level
]['Feature'].values

print(f"\nSignificant Features from Kruskal-Wallis Test (p < {kruskal_significance_level}):")
print(significant_kruskal_features)

# PermTest on sig kw features

print("\nRunning permutation tests for significant Kruskal-Wallis features...")
permutation_results_kruskal_df = run_permutation_tests(X_scaled, y, significant_kruskal_features)

print("\nPermutation Test Results for Kruskal-Wallis Features:")
print(permutation_results_kruskal_df)


# Now this is to show you the significant features after PermTesting. (Expecting to see significant drop).

significant_permutation_features = permutation_results_kruskal_df[
    permutation_results_kruskal_df['P-Value'] < 0.05
]['Feature'].values

print(f"\nSignificant Features from Permutation Tests (p < 0.05): {significant_permutation_features}")

for feature in significant_permutation_features:
    plot_feature_distribution(X_scaled, y, feature, ["Control", "Patient"])

# Mapping features after KS
mapped_features_kruskal = map_features_to_region_pairs(significant_permutation_features, region_names)

print("\nMapped Features from Kruskal-Wallis + Permutation Tests:")
for feature, region_pair in zip(significant_permutation_features, mapped_features_kruskal):
    print(f"Feature {feature}: {region_pair}")




def cluster_feature_importance(X, clusters, cluster_id):
    cluster_mask = clusters == cluster_id  # samples in the target cluster siuch as cluster 0 or 1
    other_mask = clusters != cluster_id   # Samples not in the target cluster  0 or 1 
    
    results = []
    for feature_idx in range(X.shape[1]):
        # Numbers for the feature in the target cluster vs. others
        cluster_values = X[cluster_mask, feature_idx]
        other_values = X[other_mask, feature_idx]
        
        #MannUWHIT test being performed here specifically.
        stat, p_value = mannwhitneyu(cluster_values, other_values, alternative='two-sided')
        
        # Then the mean difference, as a more positive is more associated with HC and a more negative is associated with pateint
        mean_diff = cluster_values.mean() - other_values.mean()
        
        results.append({
            'Feature': feature_idx,
            'P-Value': p_value,
            'Mean Difference': mean_diff
        })
    
    # Now make the results into a dataframe for easier access going forward.
    importance_df = pd.DataFrame(results)
    return importance_df.sort_values(by='P-Value', ascending=True)


cluster_0_importance_df = cluster_feature_importance(X_scaled, clusters, cluster_id=0)
cluster_1_importance_df = cluster_feature_importance(X_scaled, clusters, cluster_id=1)

print("Top Features for Cluster 0:")
print(cluster_0_importance_df.head(10))

print("Top Features for Cluster 1:")
print(cluster_1_importance_df.head(10))

# What does it mean to be a "top feature" within a cluster?
# Interpretation: Has a statistically significant difference in its distribution within the cluster compared to outside the cluster.
# Likely plays a critical role in defining the characteristics or behavior of that cluster relative to the rest of the dataset.


mapped_features_cluster_0 = map_features_to_region_pairs(
    cluster_0_importance_df['Feature'].values[:10],
    region_names
)
print("\nMapped Features for Cluster 0:")
for feature, region_pair in zip(cluster_0_importance_df['Feature'].values[:10], mapped_features_cluster_0):
    print(f"Feature {feature}: {region_pair}")

mapped_features_cluster_1 = map_features_to_region_pairs(
    cluster_1_importance_df['Feature'].values[:10],
    region_names
)
print("\nMapped Features for Cluster 1:")
for feature, region_pair in zip(cluster_1_importance_df['Feature'].values[:10], mapped_features_cluster_1):
    print(f"Feature {feature}: {region_pair}")

# Noise Point Analysis

# Create a mask for noise points
noise_mask = clusters == -1

# Subset the data for noise points
X_noise = X_scaled[noise_mask]
y_noise = y[noise_mask]

print(f"Number of noise points: {X_noise.shape[0]}")


from scipy.stats import mannwhitneyu

def noise_feature_importance(X, clusters, noise_label=-1):
    """
    Identify features most important for separating noise points from other points.

    Parameters:
    - X: Scaled feature matrix.
    - clusters: Cluster labels for each sample.
    - noise_label: Label representing noise points (default: -1).

    Returns:
    - A DataFrame with features ranked by p-value and mean difference.
    """
    noise_mask = clusters == noise_label
    non_noise_mask = clusters != noise_label

    results = []
    for feature_idx in range(X.shape[1]):
        # Feature values for noise and non-noise points
        noise_values = X[noise_mask, feature_idx]
        non_noise_values = X[non_noise_mask, feature_idx]

        # Perform Mann-Whitney U test
        stat, p_value = mannwhitneyu(noise_values, non_noise_values, alternative='two-sided')

        # Calculate mean difference
        mean_diff = noise_values.mean() - non_noise_values.mean()

        results.append({
            'Feature': feature_idx,
            'P-Value': p_value,
            'Mean Difference': mean_diff
        })

    # Convert to DataFrame and sort by p-value
    importance_df = pd.DataFrame(results)
    return importance_df.sort_values(by='P-Value', ascending=True)



noise_importance_df = noise_feature_importance(X_scaled, clusters)

print("Top Features for Noise Points:")
print(noise_importance_df.head(10))


mapped_noise_features = map_features_to_region_pairs(
    noise_importance_df['Feature'].values[:10],
    region_names
)

print("\nMapped Features for Noise Points:")
for feature, region_pair in zip(noise_importance_df['Feature'].values[:10], mapped_noise_features):
    print(f"Feature {feature}: {region_pair}")


num_patients_in_noise = np.sum(y_noise == 1)
num_controls_in_noise = np.sum(y_noise == 0)
print(f"Noise: {num_patients_in_noise} Patients, {num_controls_in_noise} Controls")


# Masks for non-noise patients and controls
non_noise_patients_mask = (clusters != -1) & (y == 1)
non_noise_controls_mask = (clusters != -1) & (y == 0)

# Compare noise points to non-noise patients
noise_vs_patients_df = noise_feature_importance(
    X_scaled, 
    clusters, 
    noise_label=-1
)
print("Top Features Comparing Noise to Patients:")
print(noise_vs_patients_df.head(10))

# I want to seee if a comparison for the noise points to the non-noise so controls and then patients.
noise_vs_controls_df = noise_feature_importance(
    X_scaled, 
    clusters, 
    noise_label=-1
)
print("Top Features Comparing Noise to Controls:")
print(noise_vs_controls_df.head(10))


from scipy.spatial.distance import cdist

# Cluster core samples.
core_samples_mask = (clusters != -1)  # Not the non-noise
X_core = X_scaled[core_samples_mask]
cluster_labels_core = clusters[core_samples_mask]

# Calculate distance of noise points to the nearest core point
distances = cdist(X_noise, X_core, metric='euclidean')
nearest_cluster_distances = distances.min(axis=1)  # Minimum distance to any cluster core

# Analyze the distances
print(f"Average distance of noise points to nearest cluster core: {nearest_cluster_distances.mean():.4f}")
print(f"Max distance: {nearest_cluster_distances.max():.4f}, Min distance: {nearest_cluster_distances.min():.4f}")

from sklearn.neighbors import NearestNeighbors

# Fitting the nearestneighbors to all points that were used, noise and not noise.
nbrs = NearestNeighbors(n_neighbors=10).fit(X_scaled)
distances, indices = nbrs.kneighbors(X_noise)

# Average local density for noise points
noise_density = 1 / distances.mean(axis=1)  # Inverse of mean distance to neighbors

print(f"Average density of noise points: {noise_density.mean():.4f}")
print(f"Density range: {noise_density.min():.4f} - {noise_density.max():.4f}")

# Reclustering the noise points

# UMAP red for noise
reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, random_state=42)
X_noise_reduced = reducer.fit_transform(X_noise)

# HDBSCAN on noise
hdbscan_model = hdbscan.HDBSCAN(min_cluster_size=3, min_samples=2)
noise_subclusters = hdbscan_model.fit_predict(X_noise_reduced)

print(f"Noise subclusters: {np.unique(noise_subclusters)}")

plt.figure(figsize=(10, 8))
plt.scatter(
    X_noise_reduced[noise_subclusters == -1, 0],
    X_noise_reduced[noise_subclusters == -1, 1],
    c='gray',
    label='Unclustered Noise',
    alpha=0.6
)
plt.scatter(
    X_noise_reduced[noise_subclusters == 0, 0],
    X_noise_reduced[noise_subclusters == 0, 1],
    c='blue',
    label='Noise Subcluster 0',
    alpha=0.8
)
plt.scatter(
    X_noise_reduced[noise_subclusters == 1, 0],
    X_noise_reduced[noise_subclusters == 1, 1],
    c='orange',
    label='Noise Subcluster 1',
    alpha=0.8
)
plt.title("UMAP Visualization of Noise Subclusters", fontsize=16)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.show()

subcluster_0_mask = noise_subclusters == 0
subcluster_1_mask = noise_subclusters == 1

# subclusters comparison
noise_subcluster_importance_df = noise_feature_importance(
    X_noise, 
    noise_subclusters
)
print(noise_subcluster_importance_df.head(10))

subcluster_0_patients = np.sum((y_noise == 1) & (noise_subclusters == 0))
subcluster_0_controls = np.sum((y_noise == 0) & (noise_subclusters == 0))

subcluster_1_patients = np.sum((y_noise == 1) & (noise_subclusters == 1))
subcluster_1_controls = np.sum((y_noise == 0) & (noise_subclusters == 1))

print(f"Subcluster 0: {subcluster_0_patients} Patients, {subcluster_0_controls} Controls")
print(f"Subcluster 1: {subcluster_1_patients} Patients, {subcluster_1_controls}")

mapped_features_subclusters = map_features_to_region_pairs(
    [2992, 17450, 42780, 11043, 42685, 39789, 41715, 42608, 31780, 16463],
    region_names
)

print("\nMapped Features for Noise Subclusters:")
for feature, region_pair in zip([2992, 17450, 42780, 11043, 42685, 39789, 41715, 42608, 31780, 16463], mapped_features_subclusters):
    print(f"Feature {feature}: {region_pair}")

# Patient/Control composition in Cluster 0 and Cluster 1
cluster_0_patients = np.sum((y == 1) & (clusters == 0))
cluster_0_controls = np.sum((y == 0) & (clusters == 0))

cluster_1_patients = np.sum((y == 1) & (clusters == 1))
cluster_1_controls = np.sum((y == 0) & (clusters == 1))

print("Cluster 0 Composition:")
print(f"Patients: {cluster_0_patients}, Controls: {cluster_0_controls}")

print("\nCluster 1 Composition:")
print(f"Patients: {cluster_1_patients}, Controls: {cluster_1_controls}")

# Compare to Subcluster composition
print("\nSubcluster 0 Composition:")
print(f"Patients: {np.sum((y_noise == 1) & (noise_subclusters == 0))}, Controls: {np.sum((y_noise == 0) & (noise_subclusters == 0))}")

print("\nSubcluster 1 Composition:")
print(f"Patients: {np.sum((y_noise == 1) & (noise_subclusters == 1))}, Controls: {np.sum((y_noise == 0) & (noise_subclusters == 1))}")

plt.figure(figsize=(10, 8))
plt.scatter(X_reduced[clusters == 0, 0], X_reduced[clusters == 0, 1], c='blue', label='Cluster 0', alpha=0.6)
plt.scatter(X_reduced[clusters == 1, 0], X_reduced[clusters == 1, 1], c='green', label='Cluster 1', alpha=0.6)
plt.scatter(X_noise_reduced[noise_subclusters == 0, 0], X_noise_reduced[noise_subclusters == 0, 1], c='orange', label='Subcluster 0', edgecolor='k', s=100)
plt.scatter(X_noise_reduced[noise_subclusters == 1, 0], X_noise_reduced[noise_subclusters == 1, 1], c='red', label='Subcluster 1', edgecolor='k', s=100)
plt.title("UMAP Visualization with Subclusters and Main Clusters", fontsize=16)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.show()

from scipy.stats import kruskal

features_to_test2 = [2992, 17450, 42780, 11043, 42685, 39789, 41715, 42608, 31780, 16463]

# Now I need to do a KW for each feature that is the top
kruskal_results = []
for feature in features_to_test2:
    stat, p_value = kruskal(
        X_scaled[clusters == 0, feature],  # Cluster 0
        X_scaled[clusters == 1, feature],  # Cluster 1
        X_noise[noise_subclusters == 0, feature],  # Subcluster 0
        X_noise[noise_subclusters == 1, feature]   # Subcluster 1
    )
    kruskal_results.append({'Feature': feature, 'K-Statistic': stat, 'P-Value': p_value})

# Convert results to a DataFrame for readability
import pandas as pd
kruskal_results_df = pd.DataFrame(kruskal_results)

# Display significant results
print("Kruskal-Wallis Test Results:")
print(kruskal_results_df[kruskal_results_df['P-Value'] < 0.05])


import numpy as np

def permutation_test(feature_values, group_labels, num_permutations=1000):
    """
    Perform a permutation test for a given feature across groups.
    
    Parameters:
    - feature_values: Values of the feature to test.
    - group_labels: Group labels for each value.
    - num_permutations: Number of permutations.
    
    Returns:
    - Observed statistic and p-value.
    """
    # observed statistic (difference in means across groups)
    group_means = [feature_values[group_labels == group].mean() for group in np.unique(group_labels)]
    observed_stat = max(group_means) - min(group_means)

    # PermTest
    null_distribution = []
    for _ in range(num_permutations):
        permuted_labels = np.random.permutation(group_labels)
        permuted_means = [feature_values[permuted_labels == group].mean() for group in np.unique(group_labels)]
        null_stat = max(permuted_means) - min(permuted_means)
        null_distribution.append(null_stat)

    #P-value
    p_value = np.mean(np.array(null_distribution) >= observed_stat)
    return observed_stat, p_value

# PermTests for each feature
permutation_results = []
group_labels = np.concatenate([
    np.full(np.sum(clusters == 0), 0),  # Cluster 0
    np.full(np.sum(clusters == 1), 1),  # Cluster 1
    np.full(np.sum(noise_subclusters == 0), 2),  # Subcluster 0
    np.full(np.sum(noise_subclusters == 1), 3)   # Subcluster 1
])

for feature in features_to_test2:
    feature_values = np.concatenate([
        X_scaled[clusters == 0, feature],
        X_scaled[clusters == 1, feature],
        X_noise[noise_subclusters == 0, feature],
        X_noise[noise_subclusters == 1, feature]
    ])
    observed_stat, p_value = permutation_test(feature_values, group_labels)
    permutation_results.append({'Feature': feature, 'Observed Statistic': observed_stat, 'P-Value': p_value})

# results--> DataFrame
permutation_results_df = pd.DataFrame(permutation_results)

#SigFigs
print("Permutation Test Results:")
print(permutation_results_df[permutation_results_df['P-Value'] < 0.05])


import matplotlib.pyplot as plt
import numpy as np


significant_features = [42780, 42685, 39789, 41715, 42608, 31780, 16463]

n_features = len(significant_features)
n_cols = 2 
n_rows = -(-n_features // n_cols)

# Create a figure with subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows), sharex=False, sharey=False)

# Flatten axes for easier iteration
axes = axes.flatten()

# Loop through features and plot in subplots
for i, feature in enumerate(significant_features):
    ax = axes[i]
    
    # Prepare data for the boxplot
    data = [
        X_scaled[clusters == 0, feature],  # Cluster 0
        X_scaled[clusters == 1, feature],  # Cluster 1
        X_noise[noise_subclusters == 0, feature],  # Subcluster 0
        X_noise[noise_subclusters == 1, feature]   # Subcluster 1
    ]
    
    # Create the boxplot
    box = ax.boxplot(data, patch_artist=True, labels=["Cluster 0", "Cluster 1", "Subcluster 0", "Subcluster 1"])
    
    # Customize boxplot colors
    for patch in box['boxes']:
        patch.set_facecolor("lightblue")
    
    # Overlay scatterplot for patients (x) and controls (o) with jitter
    jitter = 0.1
    for j, (group, label) in enumerate([
        (clusters == 0, "Cluster 0"),
        (clusters == 1, "Cluster 1"),
        (noise_subclusters == 0, "Subcluster 0"),
        (noise_subclusters == 1, "Subcluster 1")
    ]):
        # Get x-position for scatter points
        x_pos = j + 1
        
        # Get patient and control points for this group
        if "Cluster" in label:
            group_patients = X_scaled[(group) & (y == 1), feature]
            group_controls = X_scaled[(group) & (y == 0), feature]
        else:
            group_patients = X_noise[(group) & (y_noise == 1), feature]
            group_controls = X_noise[(group) & (y_noise == 0), feature]
        
        # Scatter points for patients (with jitter)
        ax.scatter(
            x_pos + np.random.uniform(-jitter, jitter, size=len(group_patients)), group_patients, 
            color='red', marker='x', label="Patient" if j == 0 else "", zorder=3
        )
        
        # Scatter points for controls (with jitter)
        ax.scatter(
            x_pos + np.random.uniform(-jitter, jitter, size=len(group_controls)), group_controls, 
            color='blue', marker='o', label="HC" if j == 0 else "", zorder=3
        )
    
    # Add title for each subplot
    ax.set_title(f"Feature {feature} Distribution", fontsize=12)
    ax.set_ylabel("Feature Value")
    ax.set_xlabel("Groups")
    ax.grid(alpha=0.5)
    ax.legend(loc="upper right", fontsize=8, title="Markers")

# Remove empty subplots if n_features < n_rows * n_cols
for ax in axes[n_features:]:
    ax.remove()

# Add overall spacing
fig.suptitle("Feature Distributions Across Groups", fontsize=16, y=0.95)
plt.tight_layout(rect=[0, 0.03, 1, 0.93])
plt.show()


# Include scetion about the feature importance statistical tests

