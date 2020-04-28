#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 14:58:49 2020

@author: elijah
"""

import modality as md
import knn_density_test as kdt
import numpy as np
import pandas as pd


def TestModality(data, dim, alpha=0.05, null="normal"):
    print(data.head())
    #data = np.hstack([data['X1'].values.reshape(-1, 1), data['X2'].values.reshape(-1, 1)])
    data = data[:].values.reshape(-1, dim)
    knn_dip = kdt.pointwise_test(data, plot = True)
    #dip_test = md.calibrated_diptest(data, alpha=alpha, null=null)
    #silver_test = md.calibrated_bwtest(data, alpha=alpha, null=null)
    
    return(knn_dip)
    