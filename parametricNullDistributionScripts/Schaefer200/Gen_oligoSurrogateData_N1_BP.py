#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
from scipy import stats
import gseapy as gp
from gseapy import barplot, dotplot
import time 
import matplotlib.pyplot as plt
from collections import Counter  
import SpatialAutocorrelationGSEA
import sys
import pandas as pd



# running the spatial autocorrelation correction functions 

# For reference, use the gene_sets below as last argument in spatialAutocorrCorrectedPval function
gene_expression_file = 'LH_Schaefer200GeneExpression.csv'
brain_map_file = "formatted_segstatsFilesSchaefer200/formatted_BP-Schaefer200.txt"
dist_mat_file = "Schaefer200_L_GeodesicParcelDistMat_midThickness_Left.txt"

npermutations = 1
geneOntology = "geneset_LAKE.gmt"
nh = 10
pv = 70

BPSurrMaps = SpatialAutocorrelationGSEA.generateSurrogateMaps(brain_map_file,dist_mat_file,npermutations,nh,pv)
BPSurrMaps_df = pd.DataFrame(BPSurrMaps)
BPSurrMaps_df.to_csv('BPSurrMap_forOligoSurrMapWorkflowDiag.csv',index=False,header=False)
