# YSPH / PGCoE RSV Surveillance

This repository contains analysis and data processing code, panel designs and supporting files for Respiratory Syncytial Virus (RSV) surveillance conducted at the Yale School of Public Health (YSPH) as part of the CDC's Pathogen Genomics Center of Excellence (PGCoE).

Contents
- ampschemes/       — primer sequences and amplicon locations
- analysis/         - ad-hoc analyses / QC
- background/       - background phylo data
- containers/       - docker/singularity containers
- pipelines/        - Snakemake data processing pipelines
- meta/             - clean sample metadata
- refs/             - fasta / gff references
- datasets          - processed analysis-ready data

Data access & privacy
- No individually identifiable patient data is stored in this public repo.
- All specimens were de-identified, remnant specimens used previously for diagnostic testing or IRB-approved human subjects research, in accordance with Yale University IRB-exempt protocol #2000033281.
