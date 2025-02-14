def multivariate_split(X,ar_order, valid_percent=0.15):
    X=X.copy() #input, multidimension array where each row represents variable and column.
    TS=np.shape(X)[1] #extracting the number of time series (TS) in the dataset (the dimension is the 2)
    n_vars=np.shape(X)[0] #n_vars = number of variables: extracting the number of variables  from other dimension
    val_num=int(valid_percent*TS) #the no' of time series used for validation based on the specificed percentage. Below it is stated as 0
    my_data_train=torch.zeros((TS-ar_order-val_num,ar_order,n_vars)) #creating a tensor for the training data outputs, with each output correspoidning to the next time series following the input.
    my_data_y_train=torch.zeros((TS-ar_order-val_num,1,n_vars))#creating a tensor to 
    my_data_val=torch.zeros((val_num,ar_order,n_vars))#creating a tensor as above
    my_data_y_val=torch.zeros((val_num,1,n_vars)) #storing outputs of the validation
    for i in range(TS-ar_order-val_num):#Fills my_data_train with sequences of length based on ar_order for training. The data is transposed to switch rows and columns, making columns represent variables and rows time series.
        my_data_train[i]=torch.from_numpy(X.transpose()[i:i+ar_order,:]) #Sets the target output for each training sequence, corresponding to the time series immediately following.
        my_data_y_train[i]=torch.from_numpy(X.transpose()[i+ar_order,:])
    return my_data_train, my_data_y_train, my_data_val, my_data_y_val   
#Returns the prepared training and validation datasets (both inputs and targets).
