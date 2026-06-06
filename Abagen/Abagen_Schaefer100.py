#!/usr/bin/env python
# coding: utf-8

# In[1]:


import abagen as ab
from abagen import reporting
from abagen import images
import pandas as pd


# In[7]:


atlas_Schaefer100 = ('Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii.gz')
images.check_atlas(atlas_Schaefer100)


# In[5]:


# Using downgraded version of pandas in brainSMASH environment
pd.__version__


# In[8]:


Schaefer100_expression = ab.get_expression_data(atlas_Schaefer100,
                                               probe_selection='average',
                                               donor_probes='aggregate',
                                              missing='centroids')


# In[9]:


print(Schaefer100_expression)
Schaefer100_expression.to_csv('Schaefer100_expression.csv')


# In[14]:


Schaefer100_lh_expression = Schaefer100_expression.drop(Schaefer100_expression.index[50:100],axis=0)
print(Schaefer100_lh_expression)
Schaefer100_lh_expression.to_csv('Schaefer100_lh_expression.csv')


# In[15]:


Schaefer100_lh_expression.shape

