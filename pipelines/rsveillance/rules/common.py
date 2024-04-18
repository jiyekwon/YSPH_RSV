import os

def get_value_col(search,fromcol,tocol):
        #
        pass

def get_amplicon_panel_id(wildcards):
        ampset=get_value_col(wildcards.sample, "sample","ampset")
	return ampset

def get_subsample_factor(fqsize,maxfqsize):
        factor=1/floor(fqsize/maxfqsize)
        if factor > 0.5: return 1
        else: return factor
