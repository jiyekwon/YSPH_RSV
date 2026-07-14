#!/bin/bash
#SBATCH --job-name=hmpv_smk
#SBATCH --partition=day
#SBATCH --account=pi_dmw63
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jiye.kwon@yale.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=hmpv_smk_%j.log

# Single fat allocation + LOCAL snakemake execution.
# YCRC rate-limits sbatch to 200 jobs/hour, which makes the per-job SLURM
# executor unusable for this ~4000-job workflow. Instead we take one big node
# and let snakemake run all rules as local processes across its cores
# (no per-rule sbatch -> no rate limit). apptainer works on a compute node.

set -euo pipefail
module load snakemake/8.27.0-foss-2024a git-lfs
cd /nfs/roberts/project/pi_dmw63/jk2666/YSPH_RSV/pipelines/hmpv
# CONFIG can be overridden at submit time, e.g.:  CONFIG=config_HMPV012.yaml sbatch run_hmpv.sh
ln -sf "${CONFIG:-config_hmpv.yaml}" config/config.yaml

snakemake \
  --cores 16 \
  --resources mem_mb=60000 \
  --default-resources disk_mb=1000 \
  --sdm apptainer \
  --apptainer-args "-B /nfs/roberts -B /home/jk2666" \
  --group-components mashcall=10 \
  --latency-wait 180 \
  --retries 1 \
  --keep-going \
  --rerun-trigger mtime --rerun-incomplete
