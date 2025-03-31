### Investigating hierarchal control amongst functional networks disrupted by Opioid use disorder using effective connectivity analysis.
 ---
 
This repository contains the code for the analyses described in the paper "Investigating hierarchal control amongst functional networks disrupted by opioid use disorder using effective connectivity.
Please cite this work should you use any part of this code. 

 ---
 ### Dependencies

In order to run the statistical experiments and algorithm as done in the paper, your environment will need to contains at minimum the following main packages:

- pandas 2.0.1
- numpy 1.24.4
- scipy 1.8.0
- python 3.12.3
- MATLAB_R2023b

 ---

### Large scale non-linear granger causality (LsNGC)

Large scale non-linear granger causality (LsNGC) was used for the effective connectivity analysis.
LsNGC was developed by Wismüller et al (2021), with a full repository of the original analysis and code for lsNGC [here](https://github.com/Large-scale-causality-inference/Large-scale-nonlinear-causality).

Additionally, public access to the full paper of [Wismüller et al (2021)](https://www.nature.com/articles/s41598-021-87316-6) can be found when clicking on the hyperlink.

---
After, lsNGC, a clustering model was created to investigate whether effective connectivity patterns are able to characterise healthy controls and people with opioid use disorder.

This pipeline utilizes data from 22 healthy controls (HC) and 25 methadone-dependent (MD) patients from the [Neural Correlates of Reward and Emotion in opioid dependence (NCORE) study](https://www.imperial.ac.uk/brain-sciences/research/psychiatry/ncore/). Each subject performed a series of neurocognitive tasks, and these analyses focus on the Monetary Incentive Delay (MID) and Cue Reactivity (CR) task during functional magnetic resonance imaging acquistion. 

Timeseries of task activity were extracted from each brain region from the Schaefer 2018 atlas for each subject. Following timeseries extraction, each subject had a CSV file. For these analyses, the MID task had 214 brain regions, and 292 timepoints. In the CR task, there were 214 brain regions and 322 timepoints. Many thanks to Dr Danielle Kurtin for completing the parcellation and timeseries extraction.

 ---

### Pipeline Steps

1. Make sure your timeseries are in a CSV format. For these analyses, each subject had a 214-by-214 CSV. For the MID task, the CSV were 214,292 and for the CR task, 214,322.
2. Ensure all the CSVs used for effective connectivity analysis are in the correct directory.
3. Run the Analysis_LsNGC.py script.
4. Save the outputted affinity matricies and f-statistics to your chosen directory.
5. Run the ECPermutation.m script.
6. Run the EdgeCountsAndDigraph.m script.
   * Evaluate the dominant network of influence.
   * Visualise the edge counts across each network for both conditions, HC>MD and MD>HC.
8. Run the DimensionalityReduction&Clustering.py script.
   * Reduce the dimensionality of the dataset using Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP).
   * Perform Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN).
   * Evaluate the cluster profiling, and seperation.
   * Evaluate the features important for distinguishing and seperating clusters.
   * Evaluate the noise point cluster, and investigate whether the noise points are boundaried or transitional in relation to the main clusters.

This infographic further summarises the analysis pipeline used:
<img width="588" alt="Screenshot 2025-01-28 at 14 08 59" src="https://github.com/user-attachments/assets/d1b25ec7-d7a9-4fd3-a712-0bbb4a982960" />

If you have any questions about the code, analyses, or results, please don't hesitate to email Brianna Austin (ba223@ic.ac.uk).
