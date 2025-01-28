### <ins>Investigating hierarchal control amongst functional networks disrupted by opioid use disorder using effective connectivity analysis</ins> .


This repository contains the code for the analyses described in the paper "Investigating hierarchal control amongst functional networks disrupted by opioid use disorder using effective connectivity".
Please cite this work should you use any part of this code.

 ### Dependencies

In order to run the statistical experiments and algorithm as done in the paper, your environment will need to contains at minimum the following main packages:

- pandas 2.0.1
- numpy 1.24.4
- scipy 1.8.0
- python 3.12.3
- MATLAB_R2023b


### Large scale non-linear granger causality (LsNGC)

Large scale non-linear granger causality (LsNGC) was used for the effective connectivity analysis.
LsNGC was developed by Wismüller et al (2021), with a full repository of the original analysis and code for lsNGC here: https://github.com/Large-scale-causality-inference/Large-scale-nonlinear-causality.
Additionally, public access to the full paper of [Wismüller et al (2021)](https://www.nature.com/articles/s41598-021-87316-6) can be found when clicking on the hyperlink.
After, a clustering model was created to investigate whether effective connectivity patterns are able to characterise healthy controls and people with opioid use disorder.

This pipeline utilizes data from 22 healthy controls and 25 methadone-dependent patients from the Neural Correlates of Reward and Emotion in opioid dependence (NCORE) study [1]. Each subject performed a series of neurocognitive tasks, and these analyses focus on the Monetary Incentive Delay and Cue Reactivity task during functional magnetic resonance imaging acquistion. 
Timeseries of task activity were extracted from each region of the Schaefer 2018 atlas from each subject. Following timeseries extraction, each subject has a CSV file. For these analyses, the MID task had 214 brain regions, and 292 timepoints. In the CR task, there were 214 brain regions and 322 timepoints. Many thanks to Dr Danielle Kurtin for completing the parcellation and timeseries extraction.


### Pipeline Steps

1. Make sure your timeseries are in a CSV format. For these analyses, each subject had a 214-by-214 CSV.
2. Ensure all the CSVs used for effective connectivity analysis are in the correct directory.
3. Run the Analysis_LsNGC script.
4. Save the outputted affinity matricies and f-statistics to your chosen directory.
5. Run the ECPermutation.m script.
6. Run the EdgeCountsAndDigraph script.
   * Evaluate the dominant network of influence. 6.1
   * Visualise the edge counts across each network for both conditions, HC>MD and MD>HC 6.2.
8. Run the DimensionalityReduction&Clustering script.
   * Reduce the dimensionality of the dataset using Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP).
   * Perform Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN).
   * Evaluate the cluster profiling, and seperation.
   * Evaluate the features important for distinguishing and seperating clusters.
   * Evaluate the noise point cluster, and investigate whether the noise points are boundaried or transitional in relation to the main clusters.

This infographic summarises the analysis pipeline, and the flow of work used:

<img width="391" alt="Screenshot 2024-11-12 at 19 23 43" src="https://github.com/user-attachments/assets/e346d82d-efec-47c6-812a-29382cbfdc09" />







If you have any questions about the code, analyses, or results, please don't hesitate to email Brianna Austin (ba223@ic.ac.uk).
