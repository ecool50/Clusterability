#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 22:57:07 2020

@author: elijah
"""
import randomly
import pandas as pd
import numpy as np

def Process_Data(DataFrame):
    
    DataFrame = pd.DataFrame(DataFrame)
    
    # initialize the model
    model = randomly.Rm()
    
    model.preprocess(DataFrame, min_tp=0, 
                            min_genes_per_cell=5, 
                            min_cells_per_gene=200,
                        refined=True)
    model.refining()
    
    # fit the model
    model.fit()
    
    # get the number of significant eigenvalues 
    L = model.L
    
    # get he cleaned data
    Denoised_Data = pd.DataFrame(model.return_cleaned(fdr=0.05))
    
    if model.n_components >= 2:
        # get the critical lambda
        lambda_c = model.lambda_c
    
        # get the indices of the signal eigenvectors
        indices = [i for i, x in enumerate(L) if x > lambda_c]
    
        # get the eigenvectors
        V = model.V
        V_sig = pd.DataFrame([V[:, i] for i in indices])
        V_sig = V_sig.transpose()
    
        # project the cells onto the significant eigenvectors
        Cells_Projected = model._project_cells(DataFrame, V_sig)
        
        
        # retun the projected data and the denoised data
        return[Cells_Projected, Denoised_Data]
        
    # return the dataframes
    return(Denoised_Data)