#!/usr/bin/env python
# coding: utf-8

# In[125]:


import allel 
import numpy as np

import pandas as pd
from Bio import SeqIO

import os


# In[126]:


minmaf = 0.05
samplefile = "samples.txt"
gapfile="{}_consensus_gaps.txt"

idir = "../../pipelines/rsveillance/results/ivar/"

targets = ["RSVA","RSVB"]
#targetfile = ""

outsum = "ivar_maf_summary_maf0.05.txt"
outmaf = "ivar_maf_maf0.05.txt"


# In[127]:



samples = [line.strip() for line in open(samplefile, 'r')]

#targets = [line.strip() for line in open(targetfile, 'r')]


# In[128]:


allmaf = pd.DataFrame({'MAF':list(),
             'sample':list(),
             'target':list()})


# In[130]:



for target in targets:
    gaps = pd.read_csv("{}_consensus_gaps.txt".format(target),sep="\t",index_col=False)

    for sample in samples:
        ivfile = "{}/{}_{}_ivariants.tsv".format(idir,sample,target)
        if not os.path.exists(ivfile): continue
        
        ivars = pd.read_csv(ivfile,sep="\t")
        
        ivars['MAF'] = ivars['ALT_FREQ']
        ivars.loc[(ivars['ALT_FREQ'] > 0.5),'MAF'] = 1-ivars.loc[(ivars['ALT_FREQ'] > 0.5),'ALT_FREQ']

        for i in range(0,len(gaps['start'])):
            c = "NC_001781.1"
            s = gaps.loc[i,'start']
            e = gaps.loc[i,'end']
            ivars = ivars[(ivars['POS']<s) | (ivars['POS']>e)]

        ivars = ivars[(ivars['MAF'] > minmaf)]

        if (len(ivars['MAF'])>0):
            smaf = pd.DataFrame({'MAF':ivars['MAF'],
                     'sample':[sample] * len(ivars['MAF']),
                     'target':[target] * len(ivars['MAF'])})

            allmaf = pd.concat([allmaf,smaf])


# In[131]:


sumstats = pd.merge(allmaf.groupby(['sample','target'])['MAF'].count().reset_index().rename(columns={'MAF':"n"}),
                   pd.merge(
                       allmaf.groupby(['sample','target'])['MAF'].mean().reset_index().rename(columns={'MAF':"MAFmean"}),
                       allmaf.groupby(['sample','target'])['MAF'].std().reset_index().rename(columns={'MAF':"MAFstd"}))
)




# In[124]:


sumstats.to_csv(outsum, index=False,sep="\t")
allmaf.to_csv(outmaf, index=False,sep="\t")


# In[ ]:




