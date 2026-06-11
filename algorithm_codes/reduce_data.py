import torch
import numpy as np
import pandas as pd

def load_data(folder = '/Users/psturm/Desktop/Mass conservation on manifolds/', file = "experiments_11e5_1hour_5mins_falsecombinatoricratelaws.csv"):
    
    n_points_per_experiment = 13 # 1 hours * 1 step/5 minutes + 1 for the first step
    n_steps_per_experiment = n_points_per_experiment - 1
    
    # Load reference model output:
    #   C -- concentration values
    
    # Additionally, calculate and output
    #   D -- tendency values (units of concentration)
    print("loading concentration data")
    # read dataframe from csv
    df = pd.read_csv(folder + file,  dtype='float64')

    # Get concentration values C, which are all columns except the first two
    C = df.iloc[:,2:].values

    # Get tendency values D, which are the difference between consecutive concentration values
    D = np.delete(np.diff(C,axis = 0), list(range(n_steps_per_experiment, C.shape[0]-1, n_points_per_experiment)), axis=0)

    # Get C of active species, ignoring H2O, O2, and buildup HNO3, CO, H2
    ignore_species = ['H2O', 'O2', 'HNO3', 'CO', 'H2']
    active_species_columns = [col for col in df.columns if col.split(' ')[0] not in ignore_species and col.endswith('[ppb]')]
    C_active = df[active_species_columns].values

    return df,C,D,C_active

def createIO(C,D,C_active):
    
    n_points_per_experiment = 13 # 1 hours * 1 step/5 minutes + 1 for the first step
    n_steps_per_experiment = n_points_per_experiment - 1
    
    # This function creates input X and output Y train and test sets, as well as test C and J
    # Inputs:
    #   C -- concentration values (ppb)
    #   D -- tendency values (units of concentration) over a 5 minute time step
    
    # Outputs:
    #   X_train --  scaled input values for training and validating NN   
    #   Y_train --  target values for training and validating NN  
    #   X_test  --  unscaled input values for testing/evaluating NN   
    #   Y_test  --  target values for testing/evaluating NN 
    #   C_test  --  concentrations of test data
    #   scalerX  --  scaling function used to scale NN inputs X_train and X_test
    
    print("creating input and output data")
    X = C_active   

    # delete the last step of each experiment -> Because the last step is the first step of the next experiment! TK
    X = np.delete(C_active, 
                  list(range(n_steps_per_experiment, 
                             C_active.shape[0], n_points_per_experiment)), 
                             axis=0)

    X_all = np.delete(C, 
                  list(range(n_steps_per_experiment, 
                             C_active.shape[0], n_points_per_experiment)), 
                             axis=0)

    Y = D

    #12 consecutive steps are 1 experiment, then the a new series is started!indeces 11,23,35...
    idx_withoutExp_shift = np.array(range(n_steps_per_experiment-1, X.shape[0], n_steps_per_experiment))
    diff = X_all[:-1]-X_all[1:]+D[:-1]
    diff = np.delete(diff, idx_withoutExp_shift[:-1], axis=0)
    print("max diff all:", np.max(np.abs(diff)))
    print("mean diff all:", np.mean(np.abs(diff)))

    # Create a train/test split
    split = 0.90
    trainsplit = int(split*X.shape[0])
    print("train size:", trainsplit)
    testsplit = int(round(1-split,2)*X.shape[0])
    print("test size:", testsplit)
    num_test_exps = int(testsplit/(n_steps_per_experiment))
    print("number of test experiments:", num_test_exps)

    X_train_raw = X[0:trainsplit,:]
    Y_train_raw = Y[0:trainsplit,:]
    X_test_raw = X[trainsplit:,:]
    Y_test_raw = Y[trainsplit:,:]
    C_test_raw = C[int(trainsplit*13/12):,:]
    C_test = np.reshape(C_test_raw, [num_test_exps, n_points_per_experiment, C.shape[1]])

    #to torch tensors
    X_train_raw = torch.tensor(X_train_raw)
    Y_train_raw = torch.tensor(Y_train_raw)
    X_test_raw = torch.tensor(X_test_raw)
    Y_test_raw = torch.tensor(Y_test_raw)

    return X_train_raw, Y_train_raw, X_test_raw, Y_test_raw, C_test