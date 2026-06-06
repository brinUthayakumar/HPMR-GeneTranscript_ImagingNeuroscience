#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


# 1000 spatial autocorrelation preserved lactate brain maps for left hand side of Yeo atlas
# The inputs to this function are two text files: the segstats of the phenotype for parcellation of interest
# and a distance text file between each of the parcels in the parcellation of interest, output by brainSMASH. Example in 20240220 lab book.
# Also need to format the phenotype segstats (i.e. lactate segstats) to be only the left hand side or right hand side hemispheres.
def generateSurrogateMaps(PhenotypeTxtFile,ParcelDistTxtFile,npermutations,
                          number_uniformly_spaced_distance_intervals,
                          pairwise_distance_distribution_truncation_percent):
    gen = Base(PhenotypeTxtFile, ParcelDistTxtFile, resample = True, 
               nh = number_uniformly_spaced_distance_intervals,
               pv = pairwise_distance_distribution_truncation_percent)
    nperm = npermutations
    surrogate_maps = gen(n=nperm) # Increase this to 10000 to avoid ties 

    return surrogate_maps



# For every column in the left hand side Yeo abagen expression data (i.e. every gene), use the 1000 surrogate brain maps 
def geneWiseNull(generatedSurrogateMaps,npermutations,geneOntology,gene_expression_file):
    global geneData
    surrogateCorrs_temp = []
    surrogateCorrs = []
    surrogatePrerankES_temp = []
    surrogatePrerankES = []
    surrogatePrerankTerm = []
    surrogatePrerankTerm_temp = []
    geneName = []
    #gene_expression_file = "Yeo17Split_lh_expression.csv"
    geneData = pd.read_csv(gene_expression_file)
    rnk = pd.DataFrame()

    nperm = npermutations
    # Create iterator on the csv of gene expression. This function keeps the iterator open so that the next function below simply
    # moves down the columns by one before being used to loop through. Done this way so that the first column can be skipped
    iterGenes = iter(geneData.items())
    # Skip first iteration (the 'labels' column)
    next(iterGenes)
    for series_name, series in iterGenes:
        #print(series_name)
        surrogateCorrs_temp = np.transpose(spearmanr(geneData[series_name].values,
                                                             generatedSurrogateMaps).flatten())
        
        surrogateCorrs.append(surrogateCorrs_temp)

    geneNames = pd.DataFrame(geneData.columns.values[1:geneData.shape[1]],columns = ['Genes'])
    rnk_temp = pd.DataFrame(surrogateCorrs,columns = (range(nperm)))
    rnk = pd.concat([geneNames,rnk_temp],axis = 1,)        
    
    progress_indicator = 0
    for i in range(1,np.shape(rnk)[1]):
        
        # Convert list to dataframe. More computationally efficient to go from list to dataframe then to do the above with dataframes
        
        surrogatePrerank = gp.prerank(rnk = rnk.iloc[:,[0,i]],
                                      gene_sets = geneOntology,
                                      min_size = 10,
                                      max_size = 500,
                                      ascending = False)
        surrogatePrerankES_temp = surrogatePrerank.res2d.sort_values(by = 'Term').ES
        surrogatePrerankTerm = surrogatePrerank.res2d.sort_values(by = 'Term').Term
        surrogatePrerankES.append(surrogatePrerankES_temp)
        
        
    df_surrogatePrerankES = pd.DataFrame(surrogatePrerankES)
    df_surrogatePrerankES.columns = surrogatePrerankTerm
    
    
    return df_surrogatePrerankES
    


# In[6]:


# Started March 07, 2024
#Function that calculates a spatial autocorrelation corrected pvalue using lactate hemisphere .txt file
# the geneWiseNull_output results and a gene expression text file from abagen of the appropriate hemisphere.

def spatialAutocorrCorrectedPval(geneWiseNull_output,gene_expression_file,PhenotypeTxtFile,geneOntology):

    empiricalSpearmanCorrelation_temp = []
    empiricalSpearmanCorrelation = [] 
    correctedPvalues_temp = []
    correctedPvalues = []
    
    geneData = pd.read_csv(gene_expression_file)
    iterGenes = iter(geneData.items())
    next(iterGenes)

    brainMapFile = np.loadtxt(PhenotypeTxtFile)
    for series_name, series in iterGenes:
        
        #print(GeneExpression[series_name])
        empiricalSpearmanCorrelation_temp = stats.spearmanr(brainMapFile,geneData[series_name])[0]
        empiricalSpearmanCorrelation.append(empiricalSpearmanCorrelation_temp)
        
    temp = pd.DataFrame(empiricalSpearmanCorrelation)    
    # need to create new dataframe for df.set_axis changes to be retained .... I don't get why
    df_empiricalSpearmanCorrelation = temp.set_axis(['SpearmanR'],axis = 1)
    # Create dataframe with the empirical correlations and corresponding gene label to run prerank on for empirical enrichment score
    geneNames = pd.DataFrame(geneData.columns.values[1:geneData.shape[1]],columns = ['Genes'])
    df_rnkdEmpiricalSpearmanCorrelation = pd.concat([geneNames,df_empiricalSpearmanCorrelation],axis = 1,)        
    # Empirical GSEA
    EmpiricalPrerank = gp.prerank(rnk = df_rnkdEmpiricalSpearmanCorrelation,
                                  gene_sets = geneOntology,
                                  min_size = 10,
                                  max_size = 500,
                                  ascending = False)
    
    
    #return EmpiricalPrerank.res2d.sort_values(by='Term'), EmpiricalPrerank.res2d.sort_values(by='Term').ES
    
    for geneSetID in range(len(geneWiseNull_output.columns)):

       
        n = len(geneWiseNull_output)
        correctedPvalues_temp = np.sum(np.abs(geneWiseNull_output.iloc[:,geneSetID]) > abs(
            EmpiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[geneSetID,'ES']))/n
        correctedPvalues.append(correctedPvalues_temp)
        
    geneSetTerm = EmpiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'Term']
    enrichmentScore = EmpiricalPrerank.res2d.sort_values(by = 'Term').reset_index().loc[:,'ES']
    df_correctedPvalues = pd.DataFrame(correctedPvalues)
    df_correctedLabelledPvalues = pd.concat([geneSetTerm,enrichmentScore,df_correctedPvalues],
                                             ignore_index=True,
                                             axis=1)
    
    df_correctedLabelledPvalues = df_correctedLabelledPvalues.set_axis(['geneSetTerm','EmpiricalES','CorrectedPval'],axis = 1)
    return df_correctedLabelledPvalues




# In[8]:


# Function that calculates empirical Spearman correlations and creates dataframe of results 
# 

def empiricalSpearman(brain_phenotype_file,geneExpressionCSV):
    phenotype = np.loadtxt(brain_phenotype_file)
    geneData = pd.read_csv(geneExpressionCSV)
    
    iterGenes = iter(geneData.items())
    # Skip first iteration (the 'labels' column)
    next(iterGenes)
    SpearmanCorrelations_temp = []
    SpearmanCorrelations = []
    for series_name, series in iterGenes:
        #print(series_name)
        SpearmanCorrelations_temp = np.transpose(stats.spearmanr(phenotype,geneData[series_name]))[0]
                                           
        SpearmanCorrelations.append(SpearmanCorrelations_temp)
    
    #Concatenate Spearman correlations with gene name for use as input into prerank
    temp = pd.DataFrame(SpearmanCorrelations)    
    # need to create new dataframe for df.set_axis changes to be retained .... I don't get why
    df_empiricalSpearmanCorrelation = temp.set_axis(['EmpiricalSpearman'],axis = 1)
    geneNames = pd.DataFrame(geneData.columns.values[1:geneData.shape[1]],columns = ['Genes'])
    df_rnkdEmpiricalSpearmanCorrelation = pd.concat([geneNames,df_empiricalSpearmanCorrelation],axis = 1,)      
    
    print(df_rnkdEmpiricalSpearmanCorrelation)

    return df_rnkdEmpiricalSpearmanCorrelation 

