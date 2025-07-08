#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--depth', '-d', help='depth file format')
parser.add_argument('--samples', '-s', help='samples file')
parser.add_argument('--target', '-t', help='reference name')
parser.add_argument('--out', '-o', help='outfile')

depthfile=None

args = parser.parse_args()
samplefile = args.samples
depthfile = args.depth
target = args.target 
outfile = args.out

if depthfile is None:
    depthfile = "results/ivar/{}_{}_depth.txt"

samples = open(samplefile).read().split()


alldepth = None
si=-1
for sample in samples:
    si = si+1
    sdepth = pd.read_csv(depthfile.format(sample,target),sep="\t").iloc[:,2]
    if alldepth is None:
        alldepth = pd.DataFrame(index=range(1,len(sdepth)),columns=samples)
    alldepth[sample] = sdepth

alldepth.to_csv(outfile,index=False, )

