#!/usr/bin/env python
# coding: utf-8

# In[1]:


import abagen as ab
from abagen import reporting
from abagen import images
import pandas as pd


# In[2]:


atlas_Schaefer200Yeo17Split = ('Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.nii.gz')
images.check_atlas(atlas_Schaefer200Yeo17Split)


# In[3]:


pd.__version__


# In[4]:


Schaefer200Yeo17Split_expression = ab.get_expression_data(atlas_Schaefer200Yeo17Split,
                                               probe_selection='average',
                                               donor_probes='aggregate')


# In[8]:


print(Schaefer200Yeo17Split_expression)
Schaefer200Yeo17Split_expression.to_csv('Schaefer200Yeo17Split_expression.csv')

