import numpy as np
import pandas as pd
import glob
import os

# Read all TSV files from the /data/ folder

rsvatab = pd.read_csv('data/rsv-a_metadata_2025-06-19T1926.tsv', sep='\t')
rsvbtab = pd.read_csv('data/rsv-b_metadata_2025-06-19T1927.tsv', sep='\t')
rsvatab['genotype'] = 'RSVA'
rsvbtab['genotype'] = 'RSVB'

# Concatenate all case tables
cases_df = pd.concat([rsvatab, rsvbtab], ignore_index=True)

# Ensure date column is datetime
cases_df['date'] = pd.to_datetime(cases_df['date'])

# Extract year-month for grouping
cases_df['year_month'] = cases_df['date'].dt.to_period('M')

# Summarise counts per month for each location
summary = cases_df.groupby(['location', 'year_month']).size().reset_index(name='case_count')

print(summary)