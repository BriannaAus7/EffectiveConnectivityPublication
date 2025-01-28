### Investigating hierarchal control amongst functional networks disrupted by opioid use disorder using effective connectivity analysis.

This repository contains the code for the analyses described in the paper "Investigating hierarchal control amongst functional networks disrupted by opioid use disorder using effective connectivity".
Please cite this work should you use any part of this code.

 **Dependencies**

In order to run the statistical experiments and algorithm as done in the paper, your environment will need to contains at minimum the following main packages:

- pandas 2.0.1
- numpy 1.24.4
- scipy 1.8.0
- python 3.12.3
- MATLAB_R2023b


**Large scale non-linear granger causality (LsNGC)**

Large scale non-linear granger causality (LsNGC) was used for the effective connectivity analysis.
LsNGC was developed by Wismüller et al (2021), with a full repository of the original analysis and code for lsNGC here: https://github.com/Large-scale-causality-inference/Large-scale-nonlinear-causality.
Additionally, public access to the full paper of [Wismüller et al (2021)[https://www.nature.com/articles/s41598-021-87316-6] can be found when clicking on the hyperlink.
After, a clustering model was created to investigate whether effective connectivity patterns are able to characterise healthy controls and people with opioid use disorder.


This pipeline utilizes data from 22 healthy controls and 25 methadone-dependent patients from the Neural Correlates of Reward and Emotion in opioid dependence (NCORE) study [1]. Each subject performed a series of neurocognitive tasks, and these analyses focus on the Monetary Incentive Delay and Cue Reactivity task during functional magnetic resonance imaging acquistion. 
Timeseries of task activity were extracted from each region of the Schaefer 2018 atlas from each subject. Following timeseries extraction, each subject has a CSV file. For these analyses, the MID task had 214 brain regions, and 292 timepoints. In the CR task, there were 214 brain regions and 322 timepoints.


1. Make sure your timeseries are in a CSV format. For these analyses, each subject had a 214-by-214 CSV.
2. Ensure all the CSVs used for the EC are in the correct directory
3. Run the LsNGC analysis script
4. Save the Affinity matricies and f-statistic
5. Run the EC Permutation test.
6. Run the EdgeCountsAndDigraph script
   i. Evaluate the dominant network of influence
   ii. Visualise the edge counts across each network for both conditions, HC>MD and MD>HC.
7. Run the UMAPandHDBSCAN script
   i. Reduce the dimensionality of the dataset using Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP).
   ii. Perform Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN).
   iii. Evaluate the cluster profiling, and seperation.
   iv. Evaluate the features important for distinguishing and seperating clusters.
   v. Evaluate the noise point cluster, whether the noise points are boundaries or transitional in relation to the main clusters.


**For lsNGC:**
- The functions for lsNGC can be found in the lsNGC functions folder, and these can also be found by the original authors here: https://github.com/Large-scale-causality-inference/Large-scale-nonlinear-causality
- After, the lsNGC analysis script and the functions downloaded, can be run locally.

**For the EC Permutation test: **



This infographic summarises the analysis pipeline, and the flow of work used.

If you have any questions about the code, analyses, or results, please don't hesitate to email ba223@ic.ac.uk
