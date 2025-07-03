#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('--meta', '-m', help='metadata file')
parser.add_argument('--mindate', '-D', help='minimum date (YYYY-MM-DD)')
parser.add_argument('--samples', '-s', help='samples file')
parser.add_argument('--out', '-o', help='outfile')
parser.add_argument('--datecol', '-d', help='data column (metadata)',default="date")
parser.add_argument('--namecol', '-n', help='name column (metadata)',default="name")


args = parser.parse_args()
metafile = args.meta
samplefile = args.samples
outfile = args.out
datecol=args.datecol
namecol=args.namecol

print(datecol,namecol)


#read in meta file
meta = pd.read_csv(metafile,sep="\t")
meta[datecol] = pd.to_datetime(meta[datecol],format="%Y-%m-%d")

#limit to samples in samplefile
samples = open(samplefile).read().split()
#remove reseq marker if necessary
samples = [re.sub("-1$","",s) for s in samples]

meta = meta[meta[namecol].isin(samples)]
print(samples)

#covert to days since patient 0
print(meta[datecol])
mindate = min(meta[datecol])
meta['epidays'] = pd.to_timedelta(meta[datecol]-mindate).dt.days
meta[[namecol,'epidays']]

#write to outfile
meta[[namecol,'epidays']].to_csv(outfile,index=False,header=False)
