#!/usr/bin/env python
# coding: utf-8


from brainsmash.mapgen.base import Base 
from brainsmash.workbench.geo import cortex
from brainsmash.workbench.geo import parcellate 
from brainsmash.mapgen.eval import base_fit
from brainsmash.mapgen.stats import spearmanr, pairwise_r, pearsonr
import numpy as np
from scipy import stats
from brainsmash.mapgen.stats import nonparp
import pandas as pd
import gseapy as gp
from gseapy import barplot, dotplot



def generateSurrogateMaps(PhenotypeTxtFile,ParcelDistTxtFile,npermutations,
                          number_uniformly_spaced_distance_intervals,
                          pairwise_distance_distribution_truncation_percent):
    gen = Base(PhenotypeTxtFile, ParcelDistTxtFile, resample = True, 
               nh = number_uniformly_spaced_distance_intervals,
               pv = pairwise_distance_distribution_truncation_percent)
    nperm = npermutations
    surrogate_maps = gen(n=nperm) 
    return surrogate_maps



def geneWiseNull(input_geneCSV,generatedSurrogateMaps):

    surrogate_lactate_corrs_temp = []
    surrogate_lactate_corrs = []

    geneData = pd.read_csv(input_geneCSV)
    
    iterGenes = iter(geneData.items())
    next(iterGenes)
    
    for series_name, series in iterGenes:

        surrogate_lactate_corrs_temp = np.transpose(spearmanr(geneData[series_name].values,
                                                             generatedSurrogateMaps).flatten())
        
        surrogate_lactate_corrs.append(surrogate_lactate_corrs_temp)

        
    df_surrogate_lactate_corrs = pd.DataFrame(surrogate_lactate_corrs)
    
    return df_surrogate_lactate_corrs
    
    

def spatialAutocorrCorrectedPval(geneWiseNull_output,GeneExpression,PhenotypeTxtFile):

    empiricalSpearmanCorrelation_temp = []
    empiricalSpearmanCorrelation = [] 
    correctedPvalues_temp = []
    correctedPvalues = []
    
    geneData = pd.read_csv(GeneExpression)
    iterGenes = iter(geneData.items())
    next(iterGenes)

    brainMapFile = np.loadtxt(PhenotypeTxtFile)
    for series_name, series in iterGenes:
        
        empiricalSpearmanCorrelation_temp = stats.spearmanr(brainMapFile,geneData[series_name])[0]
        empiricalSpearmanCorrelation.append(empiricalSpearmanCorrelation_temp)
        
    temp = pd.DataFrame(empiricalSpearmanCorrelation)    
    df_empiricalSpearmanCorrelation = temp.set_axis(['SpearmanR'],axis = 1)
    
    
    for rows in range(len(geneWiseNull_output)):

        n = float(len(geneWiseNull_output.loc[rows]))
        correctedPvalues_temp = np.sum(np.abs(geneWiseNull_output.loc[rows]) > abs(df_empiricalSpearmanCorrelation['SpearmanR'].loc[rows]))/n
        correctedPvalues.append(correctedPvalues_temp)
        
    df_correctedPvalues = pd.DataFrame(correctedPvalues)
    geneLabels = pd.DataFrame(list(geneData)[1:len(list(geneData))])
    
    df_correctedLabelledPvalues = pd.concat([geneLabels,df_empiricalSpearmanCorrelation,df_correctedPvalues],ignore_index=True,axis=1)
    
    df_correctedLabelledPvalues = df_correctedLabelledPvalues.set_axis(['GeneLabel','EmpiricalSpearman','CorrectedPval'],
                                                                       axis = 1)
    return df_correctedLabelledPvalues



def empiricalSpearman(brain_phenotype_file,geneExpressionCSV):
    phenotype = np.loadtxt(brain_phenotype_file)
    geneData = pd.read_csv(geneExpressionCSV)
    
    iterGenes = iter(geneData.items())

    next(iterGenes)
    SpearmanCorrelations_temp = []
    SpearmanCorrelations = []
    for series_name, series in iterGenes:

        SpearmanCorrelations_temp = np.transpose(stats.spearmanr(phenotype,geneData[series_name]))[0]
                                           
        SpearmanCorrelations.append(SpearmanCorrelations_temp)
    
    
    temp = pd.DataFrame(SpearmanCorrelations)    
    
    df_empiricalSpearmanCorrelation = temp.set_axis(['EmpiricalSpearman'],axis = 1)
    geneNames = pd.DataFrame(geneData.columns.values[1:geneData.shape[1]],columns = ['Genes'])
    df_rnkdEmpiricalSpearmanCorrelation = pd.concat([geneNames,df_empiricalSpearmanCorrelation],axis = 1,)      
    
    print(df_rnkdEmpiricalSpearmanCorrelation)

    return df_rnkdEmpiricalSpearmanCorrelation 

