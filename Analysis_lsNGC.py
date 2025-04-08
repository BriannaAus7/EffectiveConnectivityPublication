# This script is the first for the analysis, and is for the effective connectivity analysis.
# Python 3.12.3 was used.
# This code includes parts of the large scale nonlinear granger causality function described in Wismuller et al (2021), which the full link to the repository can be found in the READme.md but also in the next line::
# Wismuller et al (2021): https://www.nature.com/articles/s41598-021-87316-6.
# Headings have been provided within the code for clear instruction and details of what each code section does.
# The outputs are an affinity matrix and an f-statistic for each subject, across the timeseries, and this will be saved in your chosen directory.
# Important to note that the affinity matricies are not to be symmertised after lsNGC unlike what is described in Wismuller et al (2021) as doing so would remove the directional information cardinal to EC.

# 1. Importing libraries and preparing the input data (CSV files) to be compatabile for LsNGC analysis downstream.
import os
import numpy as np
import pandas as pd

directory_path = '/Users/briannaaustin/Desktop/lsngc/EC_Brianna(3)/MIDData/HC_MID'

def process_csv(file_path):
    df = pd.read_csv(file_path)
    print(f"Loaded {file_path}")

    # Transposing the dataframes as for my analysis, they did not have the correct orientation. If you dataframes are already in the correct format, then you can comment this section and any references to transposing out.
    df_transposed = df.T

    # 2. Save the transposed DataFrame as a new CSV file, appending 'T' to the filename just to make identifying clearer, as you will have multiple CSVs in the directory.
    new_file_path = file_path.replace('.csv', 'T.csv')
    df_transposed.to_csv(new_file_path)
    print(f"Saved transposed CSV to {new_file_path}")

    return df_transposed

# 3. Create a dictionary to store each subject's data
subjects_data = {}
for filename in os.listdir(directory_path):
    if filename.endswith('.csv') and not filename.endswith('T.csv'): 
        full_file_path = os.path.join(directory_path, filename)
        process_csv(full_file_path)

# 4. Now only loop through the transposed csvs
for filename in os.listdir(directory_path):
    if filename.endswith('T.csv'):  # Checking for transposed csvs
        full_file_path = os.path.join(directory_path, filename)
        df_transposed = pd.read_csv(full_file_path, index_col=0)
        print(f"Loaded transposed data from {full_file_path}")
        
# 5. Convert the transposed dataframe to NumPy array, and then save the brain regions extraction for later when plotting and labelling indicies
        data_array = df_transposed.to_numpy()
        brain_regions = df_transposed.index.tolist()

        
        inp_series = data_array

# 6. Store the data array, brain regions, and inp_series by subject
        subjects_data[filename] = {
            'data_array': data_array,
            'brain_regions': brain_regions,
            'inp_series': inp_series
        }

for subject, data in subjects_data.items():
    print(f"Processing data for {subject}")
    print(f"inp_series shape: {data['inp_series'].shape}")


# 7. This section of code is the lsNGC function, which can be found in full here: 
# https://github.com/Large-scale-causality-inference/Large-scale-nonlinear-causality?tab=readme-ov-file

def lsNGC(inp_series, ar_order=1, k_f=3, k_g=2, normalize=0):
    if normalize:
        X_normalized=normalize_0_mean_1_std(inp_series)
    else:
        X_normalized=inp_series.copy()
    
    X_train, Y_train , X_test, Y_test=multivariate_split(X=X_normalized,ar_order=ar_order)
    
    X_train=torch.flatten(X_train, start_dim=1)

    km= KMeans(n_clusters= k_f, max_iter= 100, random_state=123)
    km.fit(X_train)
    cent= km.cluster_centers_


    max=0 

    for i in range(k_f):
        for j in range(k_f):
            d= np.linalg.norm(cent[i]-cent[j])
            if(d> max):
                max= d
    d= max

    sigma= d/math.sqrt(2*k_f)

    sig_d=np.zeros((np.shape(X_normalized)[0],np.shape(X_normalized)[0]));
    sig=np.zeros((np.shape(X_normalized)[0],np.shape(X_normalized)[0]));

#    Z_train_label=Y_train
    for i in range(X_normalized.shape[0]):
        Z_temp=X_normalized.copy()
        Z_train, Z_train_label , _ , _=multivariate_split(X=Z_temp,ar_order=ar_order)
        Z_train=torch.flatten(Z_train, start_dim=1)
        Z_train_label=torch.flatten(Z_train_label, start_dim=1)

        # Obtain phase space Z_s by exclusing time series of of x_s
        Z_s_train, Z_s_train_label , _ , _=multivariate_split(X=np.delete(Z_temp,[i],axis=0),ar_order=ar_order)
        # Obtain phase space reconstruction of x_s
        W_s_train, W_s_train_label , _ , _=multivariate_split(X=np.array([Z_temp[i]]),ar_order=ar_order)

        # Flatten data
        Z_s_train=torch.flatten(Z_s_train, start_dim=1)
        Z_s_train_label=torch.flatten(Z_s_train_label, start_dim=1)

        W_s_train=torch.flatten(W_s_train, start_dim=1)
        W_s_train_label=torch.flatten(W_s_train_label, start_dim=1)
        # Obtain k_g number of cluster centers in the phase space W_s with k-means clustering, will have dim=(k_g * d)
        kmg= KMeans(n_clusters= k_g, max_iter= 100, random_state=123)
        kmg.fit(W_s_train)
        cent_W_s= kmg.cluster_centers_
        # Calculate activations for each of the k_g neurons
        shape= W_s_train.shape
        row= shape[0]
        column= k_g
        G= np.empty((row,column), dtype= float)
        maxg=0 

        for ii in range(k_g):
            for jj in range(k_g):
                dg= np.linalg.norm(cent_W_s[ii]-cent_W_s[jj])
                if(dg> maxg):
                    maxg= dg
        dg= maxg

        sigmag= dg/math.sqrt(2*k_g)
        if sigmag==0:
            sigmag=1
        for ii in range(row):
            for jj in range(column):
                dist= np.linalg.norm(W_s_train[ii]-cent_W_s[jj])
                G[ii][jj]= math.exp(-math.pow(dist,2)/math.pow(2*sigmag,2))
        # Generalized radial basis function
        g_ws=np.array([G[ii]/sum(G[ii]) for ii in range(len(G))])
        # Calculate activations for each of the k_f neurons 
        shape= Z_s_train.shape
        row= shape[0]
        column= k_f
        F= np.empty((row,column), dtype= float)
        for ii in range(row):
            for jj in range(column):
                cent_temp=cent.copy()
                cent_temp=np.delete(cent_temp,np.arange(jj,jj+ar_order),axis=1)
                dist= np.linalg.norm(Z_s_train[ii]-cent_temp)
                F[ii][jj]= math.exp(-math.pow(dist,2)/math.pow(2*sigma,2))
        # Generalized radial basis function
        f_zs=np.array([F[ii]/sum(F[ii]) for ii in range(len(F))])

        # Prediction in the presence of x_s
        num_samples=f_zs.shape[0]

        f_new=np.concatenate((0.5*f_zs,0.5*g_ws),axis=1)
        GTG= np.dot(f_new.T,f_new)
        GTG_inv= np.linalg.pinv(GTG)
        fac= np.dot(GTG_inv,f_new.T)
        W_presence= np.dot(fac,Z_train_label)
        
        prediction_presence= np.dot(f_new,W_presence)
        error_presence=prediction_presence-np.array(Z_train_label)
        sig[i,:]=np.diag(np.cov(error_presence.T))

        # Prediction without x_s
        GTG= np.dot(f_zs.T,f_zs)
        GTG_inv= np.linalg.pinv(GTG)
        fac= np.dot(GTG_inv,f_zs.T)
        W_absence= np.dot(fac,Z_train_label)

        prediction_absence= np.dot(f_zs,W_absence)
        error_absence=prediction_absence-np.array(Z_train_label)
        sig_d[i,:]=np.diag(np.cov(error_absence.T))
    # Comupte the Granger causality index

    Aff=np.log(np.divide(sig_d,sig))
    Aff=(Aff>0)*Aff
    np.fill_diagonal(Aff,0)
    f_stat=calc_f_stat(sig_d, sig, n=num_samples+1, pu=k_f+k_g, pr=k_f)
    np.fill_diagonal(f_stat,0)
    
    return Aff, f_stat

# 8. Create a results dictionary to store each subject's Affinity matrix and F-statistic.
results = {}

for subject, data in subjects_data.items():
        Aff, f_stat = lsNGC(data['inp_series'])
        results[subject] = {
            'Aff': Aff,
            'f_stat': f_stat
        }


# 9. This is for the results dictionary and subjects_data with brain_regions, again in your chosen directory.
base_directory = '/Users/briannaaustin/Desktop/lsngc/EC_Brianna(3)/MIDData/HC_MID'
if not os.path.exists(base_directory):
    os.makedirs(base_directory)  # Double check that the directory exists.

for subject, result in results.items():

# 10. Retrieve the brain regions for each subject (The brain regions were saved earlier in the script before DataFrame --> NumPy array conversion).
    brain_regions = subjects_data[subject]['brain_regions']

# 11. Convert the affinity matrices to DataFrames, setting brain regions as index and columns
    aff_df = pd.DataFrame(result['Aff'], index=brain_regions, columns=brain_regions)
    f_stat_df = pd.DataFrame(result['f_stat'], index=brain_regions, columns=brain_regions)

# 12. Create file paths for both the affinity matricies and f-statistic using your base directory.
    aff_file_path = f'{base_directory}/{subject}_Aff.csv'
    f_stat_file_path = f'{base_directory}/{subject}_F_stat.csv'

# 13. Save the affinity matricies and f-statistic as CSV files to the paths.
    aff_df.to_csv(aff_file_path)
    f_stat_df.to_csv(f_stat_file_path)

# 14. A print statement acting as a sanity check to ensure the CSV files have been saved in the correct place.
    print(f"Saved Aff matrix to {aff_file_path}")
    print(f"Saved F-stat matrix to {f_stat_file_path}")

# Now, the next step is to run the ECPermutation.m script to identify which edges had significantly different EC between MD vs HC participants. 

