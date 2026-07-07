#!/usr/bin/env bash
#
# fetch_batch.sh — download a YCGA sequencing batch and stage it for the
#                  rsveillance pipeline.
#
# Handles the three ways YCGA delivers data:
#   (a) S3 presigned URL   https://ycgasequencearchive.s3.amazonaws.com/....tar?X-Amz-...
#   (b) S3 object path      s3://ycgasequencearchive/archive/.../Project_XXXX.tar
#   (c) HTTP directory      http://fcb.ycga.yale.edu:3010/<token>/sample_dir_00000000/
#       (the recursive "wget -r -np -nH -A '*.fastq.gz' ..." style link)
#
# Two ways to stage:
#   * per-batch (default)      -> its own reads/ + <batch>_samples.txt + config
#   * accumulate  (-c NAME)    -> add this batch into a shared collection NAME,
#                                 merging into one reads/ + one NAME_samples.txt.
#                                 Run several batches with the same -c NAME, then
#                                 run the pipeline once on the collection.
#
# For a batch it will:
#   1. download into  $SCRATCH_BASE/<batch>/           (tarball, or recursive fastq pull)
#   2. extract if it's a tarball
#   3. build  <target>/reads/<SAMPLE>/  dirs of R1/R2 symlinks, sample id taken
#      from the fastq filename and made pipeline-legal (no underscores):
#         Yale-RSV-1213_S265_L007_R1_001.fastq.gz -> "Yale-RSV-1213"
#         RV024_cDNA_S12_L007_R1_001.fastq.gz     -> "RV024-cDNA"
#   4. write / update  $SAMPLE_DIR/<target>_samples.txt   (sorted, de-duplicated)
#   5. validate exactly one R1 and one R2 per sample
#   6. write  config/config_<target>.yaml  (copy of the template with readdir /
#      samples / runname swapped) so the collection is ready to run
#   ( <target> = the collection name with -c, otherwise the batch name )
#
# Usage:
#   ./fetch_batch.sh [-c COLLECTION] <batch_name> "<link>"
#   ./fetch_batch.sh -s [-c COLLECTION] <batch_name>          # stage-only: skip download,
#                                                             # reuse data already in <batch>/
#
# Examples:
#   # one-off batch, its own run:
#   ./fetch_batch.sh RSV_batch01 "https://ycgasequencearchive.s3.amazonaws.com/.../Project_XXXX.tar?X-Amz-..."
#
#   # accumulate several batches into collection "RSV_all", then run once:
#   ./fetch_batch.sh -c RSV_all RSV_batch00 "http://fcb.ycga.yale.edu:3010/tok/sample_dir_000021410/"
#   ./fetch_batch.sh -c RSV_all RSV_batch01 "https://ycgasequencearchive.s3.amazonaws.com/.../P1.tar?X-Amz-..."
#   ./fetch_batch.sh -c RSV_all RSV_batch02 s3://ycgasequencearchive/archive/.../P2.tar
#   # ...then:  cd ..;  CONFIG=config_RSV_all.yaml sbatch run_rsv025.sh
#
# Pass just the LINK (in quotes). For the wget style, don't paste the whole
# "wget -r -np ..." command — only the URL. A stray leading backslash is stripped.
#
# Notes:
#   - Run in an interactive allocation, not the login node:
#       salloc -p day -c 4 --mem 16G -t 3:00:00
#   - s3:// needs the AWS CLI ( module load awscli ) and, if not public, creds:
#       export AWS_ACCESS_KEY_ID=...  AWS_SECRET_ACCESS_KEY=...  [AWS_SESSION_TOKEN=...]
#   - Override defaults without editing the script:
#       SCRATCH_BASE=/path  SAMPLE_DIR=/path  CONFIG_TEMPLATE=/path/config_x.yaml  ./fetch_batch.sh ...

set -euo pipefail

usage() { sed -n '2,60p' "$0" | sed 's/^#//; s/^ //'; }

# ---------------------------------------------------------------- config ----
SCRATCH_BASE="${SCRATCH_BASE:-/home/jk2666/scratch_pi_dmw63/jk2666}"
SAMPLE_DIR="${SAMPLE_DIR:-/home/jk2666/project_pi_dmw63/jk2666/rsv_yale}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PIPELINE_DIR/config"
CONFIG_TEMPLATE="${CONFIG_TEMPLATE:-$CONFIG_DIR/config_RSV025.yaml}"

# ---------------------------------------------------------------- args ------
COLLECTION=""
STAGE_ONLY=0
while getopts "c:sh" opt; do
    case "$opt" in
        c) COLLECTION="$OPTARG" ;;
        s) STAGE_ONLY=1 ;;
        h) usage; exit 0 ;;
        *) echo "usage: $0 [-c COLLECTION] [-s] <batch_name> [\"<link>\"]" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [[ "$STAGE_ONLY" -eq 1 ]]; then
    if [[ $# -lt 1 ]]; then
        echo "usage (stage-only): $0 -s [-c COLLECTION] <batch_name>" >&2; exit 1
    fi
    BATCH="$1"; SRC="(stage-only)"
else
    if [[ $# -ne 2 ]]; then
        echo "usage: $0 [-c COLLECTION] <batch_name> \"<link>\"" >&2; exit 1
    fi
    BATCH="$1"
    SRC="$2"
    SRC="${SRC#\\}"                         # tolerate a stray leading backslash from a pasted command
    SRC="${SRC#"${SRC%%[![:space:]]*}"}"    # trim leading whitespace
fi

TARGET="${COLLECTION:-$BATCH}"          # where reads/, samples list and config live
BATCHDIR="$SCRATCH_BASE/$BATCH"         # each batch downloads to its own dir (provenance)
READS="$SCRATCH_BASE/$TARGET/reads"     # staging dir (shared in accumulate mode)
TARBALL="$BATCHDIR/_download.tar"        # fixed name: sidesteps presigned "?..." suffixes
SAMP="$SAMPLE_DIR/${TARGET}_samples.txt"
OUTCONFIG="$CONFIG_DIR/config_${TARGET}.yaml"

mkdir -p "$BATCHDIR" "$SAMPLE_DIR" "$SCRATCH_BASE/$TARGET"

if [[ -n "$COLLECTION" ]]; then
    echo ">> accumulate mode: batch '$BATCH' -> collection '$COLLECTION'"
else
    echo ">> per-batch mode: '$BATCH'"
fi

# ---------------------------------------------------------------- 1. detect + download
if [[ "$STAGE_ONLY" -eq 1 ]]; then
    echo ">> [1/6] stage-only: skipping download; using existing data in $BATCHDIR"
    [[ -d "$BATCHDIR" ]] || { echo "!! $BATCHDIR not found — nothing to stage" >&2; exit 1; }
    MODE=stageonly
else
    if   [[ "$SRC" == s3://* ]];       then MODE=s3tar
    elif [[ "$SRC" =~ \.tar($|\?) ]];  then MODE=httptar
    else                                    MODE=wgetdir      # recursive HTTP fastq directory
    fi
    echo ">> [1/6] download mode: $MODE"

    case "$MODE" in
      s3tar)
        command -v aws >/dev/null || { echo "!! aws CLI not found — run: module load awscli" >&2; exit 1; }
        aws s3 cp "$SRC" "$TARBALL"
        ;;
      httptar)
        curl -f -L -C - -o "$TARBALL" "$SRC"     # -f fail on errors, -L redirects, -C - resume
        ;;
      wgetdir)
        # Reject index/robots pages rather than -A '*.fastq.gz': YCGA directory
        # listings have no .html extension, so -A discards them before wget can
        # follow links into the Sample_*/ subdirs -> 0 fastqs. -R lets it traverse.
        wget -e robots=off -r -np -nH -R 'index.html*,robots.txt,*.tmp' -P "$BATCHDIR/download" "$SRC"
        ;;
    esac
fi

# ---------------------------------------------------------------- 2. extract (tar modes only)
if [[ "$MODE" == s3tar || "$MODE" == httptar ]]; then
    if ! file "$TARBALL" | grep -qi 'tar archive'; then
        echo "!! downloaded file is not a tar archive — likely an expired/invalid link:" >&2
        head -c 400 "$TARBALL" >&2; echo >&2
        exit 1
    fi
    echo ">> [2/6] extracting tarball"
    tar -xf "$TARBALL" -C "$BATCHDIR"
else
    echo ">> [2/6] (no tar to extract)"
fi

# ---------------------------------------------------------------- 3. locate fastqs
echo ">> [3/6] locating fastq.gz files"
fq_list="$(mktemp)"
find "$BATCHDIR" -type d -name reads -prune -o -type f -name '*.fastq.gz' -print | sort > "$fq_list"
n_fq=$(wc -l < "$fq_list")
if [[ "$n_fq" -eq 0 ]]; then
    echo "!! no *.fastq.gz found under $BATCHDIR" >&2
    rm -f "$fq_list"; exit 1
fi
echo "   found $n_fq fastq files in this batch"

# ---------------------------------------------------------------- 4. stage reads/ + list
echo ">> [4/6] staging reads/<sample>/ symlinks and updating sample list"
if [[ -z "$COLLECTION" ]]; then
    rm -rf "$READS"                     # per-batch: clean rebuild (reads/ holds only our symlinks)
    : > "$SAMP"
fi
mkdir -p "$READS"; touch "$SAMP"

batch_tmp="$(mktemp)"
while IFS= read -r fq; do
    base=$(basename "$fq")
    # strip extension, then everything from _S<n>_ on (Illumina), else the _R1/_R2 suffix
    s=$(printf '%s' "$base" | sed -E 's/\.fastq\.gz$//; s/_S[0-9]+_.*//; s/_R[12]([._].*)?$//')
    s=${s//_/-}                         # RV024_cDNA -> RV024-cDNA  (no underscores)
    mkdir -p "$READS/$s"
    ln -sfn "$fq" "$READS/$s/$base"
    echo "$s" >> "$batch_tmp"
done < "$fq_list"
rm -f "$fq_list"
n_new=$(sort -u "$batch_tmp" | wc -l)
cat "$batch_tmp" >> "$SAMP"
sort -u "$SAMP" -o "$SAMP"              # merge + de-dup (idempotent if a batch is re-run)
rm -f "$batch_tmp"

# provenance: record what went into this collection/batch
printf '%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M')" "$BATCH" "$n_new" "$SRC" \
    >> "$SCRATCH_BASE/$TARGET/batches.log"

# ---------------------------------------------------------------- 5. validate
echo ">> [5/6] validating one R1 + one R2 per sample"
bad=0
while read -r s; do
    r1=$(find -L "$READS/$s" -name '*_R1*.fastq.gz' | wc -l)
    r2=$(find -L "$READS/$s" -name '*_R2*.fastq.gz' | wc -l)
    if [[ "$r1" -ne 1 || "$r2" -ne 1 ]]; then
        echo "   WARN  $s  R1=$r1 R2=$r2"
        bad=1
    fi
done < "$SAMP"

# ---------------------------------------------------------------- 6. (re)generate config
echo ">> [6/6] writing config for target '$TARGET'"
if [[ -f "$CONFIG_TEMPLATE" ]]; then
    mkdir -p "$CONFIG_DIR"
    sed -E "s|^readdir:.*|readdir: '$READS'|; s|^samples:.*|samples: '$SAMP'|; s|^runname:.*|runname: '$TARGET'|" \
        "$CONFIG_TEMPLATE" > "$OUTCONFIG"
    echo "   wrote $OUTCONFIG"
else
    echo "   !! template $CONFIG_TEMPLATE not found — skipped config generation" >&2
    OUTCONFIG=""
fi

# ---------------------------------------------------------------- summary
echo
echo "================ batch staged ================"
echo "batch        : $BATCH  (+$n_new samples)"
echo "target       : $TARGET"
echo "reads dir    : $READS  ($(wc -l < "$SAMP") samples total)"
echo "sample list  : $SAMP"
[[ -n "$OUTCONFIG" ]] && echo "config       : $OUTCONFIG"
echo
if [[ -n "$COLLECTION" ]]; then
    echo "accumulating — add more batches with:  $0 -c $COLLECTION <batch> \"<link>\""
    echo "when the collection is complete, run the pipeline:"
else
    echo "to run the pipeline:"
fi
[[ -n "$OUTCONFIG" ]] && {
    echo "    cd $PIPELINE_DIR"
    echo "    CONFIG=config_${TARGET}.yaml sbatch run_rsv025.sh"
}
if [[ "$bad" -eq 0 ]]; then
    echo "all samples have exactly one R1/R2  [OK]"
else
    echo "SOME SAMPLES NEED ATTENTION (see WARN above) — e.g. a multi-lane sample"
    echo "with >1 R1 will break cp_local_fq; those lanes need concatenating first."
fi
