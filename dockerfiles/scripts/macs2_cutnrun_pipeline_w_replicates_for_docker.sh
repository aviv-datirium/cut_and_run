#!/bin/bash

set -o pipefail

cat <<'BANNER'

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ CUT&RUN PIPELINE (Paired-End) - Replicates · E. coli Spike-in · MACS2 Peaks ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
 This Bash workflow trims, aligns, deduplicates, filters, and peak-calls
 replicated CUT&RUN libraries. It supports an optional IgG/empty-vector control
 set, scales coverage tracks by E. coli spike-in, and writes three MACS2 peak
 tiers:

   • replicate/   – one peak set per treatment replicate
   • pooled/      – each treatment replicate vs pooled control
   • merged/      – treatment-merged vs control-merged

 OUTPUT TREE
 ───────────
   output/
     fastqc_reports/           ⇠   FastQC HTML per FASTQ
     bigwig_bedgraphs/         ⇠   *.bw + (optionally) *.bedgraph
     macs2_peaks/
        replicate/*.narrowPeak *.xls *.bed
        pooled/*.narrowPeak *.xls *.bed
        merged/*.narrowPeak *.xls *.bed
     annotated_peaks/          ⇠   *.annotated.{bed,tsv}
     preseq/
     logs/                     ⇠   STAR, Picard, Trim, etc.

 MAJOR STEPS
 ───────────
  0  Read paths & parameters (config.json)
  1  FastQC on raw FASTQs
  2  Trimming (Trimmomatic) – all treatment & control replicates
  3  STAR   : host genome (hg38/mm10/…) → per-replicate BAM + index
  4  STAR   : E. coli (spike-in) → per-replicate BAM + index
  5  Picard : add-RG + duplicate removal
  6  Fragment-size filtering (histone/TF/≤1 kb; user selection)
  7  MACS2  : replicate, merged, pooled peak calling (BAMPE mode)
  8  Spike-in scale factors – host/spike read ratio per replicate
  9  BedGraph + BigWig generation for viewing (scaled if factors exist)
 10  Peak-to-gene annotation (bedtools intersect)
 11  Preseq complexity estimation and plotting

 USAGE
 ─────
   • Edit the config.json with absolute paths and replicate lists.
   • Run: bash macs2_cut_and_run_pipeline_w_replicates.sh (no CLI arguments)
   • Logs stream to stdout *and* to output/logs/pipeline.log

 REQUIREMENTS
 ────────────
   bash ≥4 · samtools ≥1.10 · bedtools ≥2.28 · STAR ≥2.7 · Java ≥17
   Trim Galore ≥0.6.10 · Picard ≥2.18 · MACS2 ≥2.2 · cutadapt ≥4.1
   R 4.x (for plotting) · GNU coreutils/awk · PIGZ · GNU Parallel

 CITATIONS
 ────────────
   O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
     ;login: The USENIX Magazine, February 2011:42-47.
  
BANNER

###############################################################################
# 0  CONFIG + PATHS                                                           #
###############################################################################
# Make script directory discoverable (portable execution from any location)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

CONFIG_FILE="${1:-config_for_docker.json}"
shift || true

# Move to the directory containing the config file
cd "$(dirname "$CONFIG_FILE")" || exit 1

# Use the basename so jq works with the local file
CONFIG_FILE="$(basename "$CONFIG_FILE")"

# Directory containing the config file — assumed to be project root
CONFIG_DIR=$(pwd)

# Use it to construct all other relative paths
ALIGNMENT_DIR="$CONFIG_DIR/$(jq -r '.alignment_dir' "$CONFIG_FILE")"
OUTPUT_DIR="$CONFIG_DIR/$(jq -r '.output_dir' "$CONFIG_FILE")"
LOG_DIR="$CONFIG_DIR/$(jq -r '.log_dir' "$CONFIG_FILE")"
REFERENCE_GENOME="$CONFIG_DIR/$(jq -r '.reference_genome' "$CONFIG_FILE")"
ECOLI_INDEX="$CONFIG_DIR/$(jq -r '.ecoli_index' "$CONFIG_FILE")"
ANNOTATION_GENES="$CONFIG_DIR/$(jq -r '.annotation_genes' "$CONFIG_FILE")"
CHROM_SIZE="$CONFIG_DIR/$(jq -r '.chrom_sizes' "$CONFIG_FILE")"

GENOME_SIZE_STRING=$(jq -r '.genome_size'      "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' "$CONFIG_FILE")
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
NUM_THREADS=$(jq    -r '.num_threads'          "$CONFIG_FILE")
NUM_PARALLEL_THREADS=$(jq    -r '.num_parallel_threads'          "$CONFIG_FILE")

TREAT_R1=($(jq -r '.samples.treatment[]?.r1' "$CONFIG_FILE" | sed -E "s|^(/)|\1|; t; s|^|$CONFIG_DIR/|"))
TREAT_R2=($(jq -r '.samples.treatment[]?.r2' "$CONFIG_FILE" | sed -E "s|^(/)|\1|; t; s|^|$CONFIG_DIR/|"))
CTRL_R1=($(jq -r '.samples.control[]?.r1 // empty' "$CONFIG_FILE" | sed -E "s|^(/)|\1|; t; s|^|$CONFIG_DIR/|"))
CTRL_R2=($(jq -r '.samples.control[]?.r2 // empty' "$CONFIG_FILE" | sed -E "s|^(/)|\1|; t; s|^|$CONFIG_DIR/|"))


###############################################################################
# 1  TOOL LOCATIONS                                                           #
###############################################################################
PICARD_CMD="$(command -v picard)"
if [[ -z $PICARD_CMD ]]; then #Fallback in case command -v picard fails
  echo "Error: Picard not found in PATH." >&2
  exit 1
fi
export PICARD_CMD ALIGNMENT_DIR LOG_DIR NUM_THREADS
FASTQC_BIN="$(command -v fastqc)"
STAR_BIN="$(command -v STAR)"
MACS2="macs2"

###############################################################################
# 2  FOLDERS + LOGGER                                                         #
###############################################################################
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR"

# timestamp the old log (if any) before we start the new one (macOS/Linux)
RUN_TS=$(date -u +'%Y%m%d_%H%M%S')   # e.g. 20250624_104832 
if [[ -f $LOG_DIR/pipeline.log ]]; then
    mv "$LOG_DIR/pipeline.log" \
       "$LOG_DIR/pipeline_${RUN_TS}.log"
fi
# create a fresh empty logfile for this run
: > "$LOG_DIR/pipeline.log"
                               
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports"
SPIKE_DIR="$ALIGNMENT_DIR/spikein"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks"
export PEAK_DIR
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs"
ANN_DIR="$OUTPUT_DIR/annotated_peaks"
for d in "$FASTQC_DIR" "$SPIKE_DIR" "$PEAK_DIR"/{replicate,merged,pooled} \
         "$BW_DIR" "$ANN_DIR"; do mkdir -p "$d"; done

log(){ printf '[%(%F %T)T] %-10s %-15s %s\n' -1 "$1" "$2" "${*:3}" \
      | tee -a "$LOG_DIR/pipeline.log"; }

log START Paramaters "Config=$CONFIG_FILE"

###############################################################################
# 3  HELPER FUNCTIONS                                                         #
###############################################################################
# check for essential binaries, and fail fast if missing
check_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: $1 not found in PATH." >&2
    exit 1
  }
}

for tool in fastqc trim_galore STAR samtools bedtools bedGraphToBigWig macs2 picard; do
  check_tool "$tool"
done

get_sample_basename() {
  # keep EVERYTHING up to the trailing _R1/_R2/_R1_001/_R2_001 token
  local b=${1##*/}               # strip directory
  b=${b%.fastq.gz}; b=${b%.fq.gz}
  b=${b%_R1_001}; b=${b%_R2_001}
  b=${b%_R1};     b=${b%_R2}
  b=${b%_1};      b=${b%_2}
  printf '%s\n' "$b"
}

trim_one_pair () {                # $1 R1  $2 R2  $3 SAMPLE multi-thread Cutadapt + pigz
  local R1="$1" R2="$2" BASE="$3"

  log Trim "$BASE" start

  trim_galore --paired --quality 20 --phred33 \
              --gzip \
              --cores "$NUM_PARALLEL_THREADS" \
              --length 25 \
              --output_dir "$ALIGNMENT_DIR" \
              "$R1" "$R2" \
              > "$LOG_DIR/trim_${BASE}.log" 2>&1

  # detect Trim Galore! outputs
  local VAL1 VAL2
  VAL1=$(find "$ALIGNMENT_DIR" -maxdepth 1 -name "${BASE}*val_1.fq.gz" | head -n1)
  VAL2=$(find "$ALIGNMENT_DIR" -maxdepth 1 -name "${BASE}*val_2.fq.gz" | head -n1)

  if [[ -s $VAL1 && -s $VAL2 ]]; then
      mv "$VAL1" "$ALIGNMENT_DIR/${BASE}_trimmed_R1.fq.gz"
      mv "$VAL2" "$ALIGNMENT_DIR/${BASE}_trimmed_R2.fq.gz"
      log Trim "$BASE" ok
  else
      log Trim "$BASE" FAIL "no trimmed FASTQs produced"
      return 1
  fi
}

run_star () {                     # R1  R2  OUTPREFIX  GENOMEDIR
  "$STAR_BIN" --runThreadN "$NUM_THREADS" \
              --genomeDir "$4" \
              --readFilesIn "$1" "$2" \
              --readFilesCommand zcat \
              --outSAMtype BAM SortedByCoordinate \
              --outFileNamePrefix "$3" \
              > "$LOG_DIR/$(basename "$3")STAR.log" \
              2> "$LOG_DIR/$(basename "$3")STAR_err.log"
}

get_unmapped_mates () {
  local base="$ALIGNMENT_DIR/$1"
  local m1="${base}.Unmapped.out.mate1"
  local m2="${base}.Unmapped.out.mate2"
  [[ -s $m1 && -s $m2 ]] && printf '%s %s\n' "$m1" "$m2"
}

run_fastqc () {                     # $1=R1  $2=R2  $3=SAMPLE
  local R1=$1 R2=$2 SAMPLE=$3
  log FastQC "$SAMPLE" start
  $FASTQC_BIN --quiet --extract -o "$FASTQC_DIR" "$R1" "$R2" \
      >  "$LOG_DIR/fastqc_${SAMPLE}.log" \
      2> "$LOG_DIR/fastqc_${SAMPLE}.err"
  if [[ $? -eq 0 ]]; then
      log FastQC "$SAMPLE" OK
  else
      log FastQC "$SAMPLE" FAIL  "see fastqc_${SAMPLE}.err"
  fi
}
export -f run_fastqc

# -----------------------------------------------------------------------------
# merge_bams  <out.bam>  <sampleName> ...                                     -
#   • 0 inputs  → skip with log                                               -
#   • 1 input   → ln -f (fast)                                                -
#   • ≥2 inputs → samtools merge + index                                      -
# -----------------------------------------------------------------------------
merge_bams () {
  local out=$1; shift
  local inputs=()

  # collect only existing filtered BAMs
  for s in "$@"; do
    bam="$ALIGNMENT_DIR/${s}.dedup.filtered.bam"
    [[ -s $bam ]] && inputs+=("$bam")
  done

  case ${#inputs[@]} in
    0)
        log MERGE "$(basename "$out")" "skip (no input BAMs)"
        return 1 ;;
    1)
        ln -f "${inputs[0]}" "$out"
        samtools index -@ "$NUM_THREADS" "$out" ;;
    *)
        samtools merge -@ "$NUM_THREADS" -f "$out" "${inputs[@]}" \
          && samtools index -@ "$NUM_THREADS" "$out" \
          && log MERGE "$(basename "$out")" "done (${#inputs[@]} inputs)" \
          || { log MERGE "$(basename "$out")" "FAIL (samtools merge)"; return 1; }
  esac
}

bam_to_bedgraph(){ bedtools genomecov -ibam "$1" -bg -pc | sort -k1,1 -k2,2n > "$2"; }
read_count(){ samtools view -c -F 2304 "$1"; }

###############################################################################
# 4  GET BASENAMES + GENOME SIZE                                              #
###############################################################################
TREAT_NAMES=(); for r1 in "${TREAT_R1[@]}"; do TREAT_NAMES+=( "$(get_sample_basename "$r1")" ); done
CTRL_NAMES=();  for r1 in "${CTRL_R1[@]}";  do CTRL_NAMES+=( "$(get_sample_basename "$r1")" ); done
SAMPLES=("${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}")

case $GENOME_SIZE_STRING in
  hs) GENOME_SIZE=2913022398 ;;  mm) GENOME_SIZE=2652783500 ;;
  dm) GENOME_SIZE=165000000  ;;  ce) GENOME_SIZE=1000000000 ;;
  sc) GENOME_SIZE=12000000   ;;  *)  GENOME_SIZE=$CUSTOM_GENOME_SIZE ;;
esac

# In the main script call it conditionally
if [[ $RUN_ONLY_BLOCK == "yes" ]]; then
    run_block
    exit 0
fi

###############################################################################
# 5  FASTQC  – per-sample logging (start / ok / FAIL)                         #
###############################################################################
log FastQC Conditions "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}"

ALL_FASTQS=( "${TREAT_R1[@]}" "${TREAT_R2[@]}"
             "${CTRL_R1[@]:-}" "${CTRL_R2[@]:-}" )

run_fastqc () {                       # $1 = FASTQ
  local fq="$1" base
  base=$(basename "$fq")
  log FastQC "$base" start
  if "$FASTQC_BIN" --extract -o "$FASTQC_DIR" "$fq" \
        >>"$LOG_DIR/fastqc_${base}.log" 2>&1 ; then
      log FastQC "$base" done
  else
      log FastQC "$base" FAIL
  fi
}
export -f run_fastqc log
export LOG_DIR FASTQC_DIR FASTQC_BIN

if command -v parallel >/dev/null 2>&1; then
  log FastQC ALL "running ${#ALL_FASTQS[@]} files with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           run_fastqc ::: "${ALL_FASTQS[@]}" \
    2>&1 | tee -a "$LOG_DIR/fastqc_parallel.log"
else
  log FastQC ALL "GNU parallel not found – running serially"
  for fq in "${ALL_FASTQS[@]}"; do run_fastqc "$fq"; done
fi

###############################################################################
# 6  TRIMMING  – per-sample logging (start / done)                            #
###############################################################################
log Trim ALL "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}  (Trim Galore! --cores $NUM_PARALLEL_THREADS)"

# ── treatment replicates ────────────────────────────────────────────────────
for i in "${!TREAT_R1[@]}"; do
  trim_one_pair "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}" || continue
done

# ── control replicates (if any) ─────────────────────────────────────────────
for i in "${!CTRL_R1[@]}"; do
  trim_one_pair "${CTRL_R1[$i]}" "${CTRL_R2[$i]}" "${CTRL_NAMES[$i]}" || continue
done

###############################################################################
# 7  HOST GENOME ALIGNMENT                                                    #
###############################################################################
for n in "${SAMPLES[@]}"; do
  log STARhost "$n" start
  run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" \
           "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           "$ALIGNMENT_DIR/${n}." "$REFERENCE_GENOME" \
           --outReadsUnmapped Fastx          # ← NEW OPTION
  log STARhost "$n" done
done

###############################################################################
# 8  SPIKE-IN ALIGNMENT (E. coli)                                             #
###############################################################################
for n in "${SAMPLES[@]}"; do
  # locate mate FASTQs produced in Step 8
  read -r U1 U2 < <(get_unmapped_mates "$n")

  if [[ -z $U1 ]]; then
    log SPIKE "$n" "skip (no unmapped mates)"
    continue
  fi

  log SPIKE "$n" start
  run_star "$U1" "$U2" \
           "$SPIKE_DIR/${n}_ecoli_" "$ECOLI_INDEX"
  mv "$SPIKE_DIR/${n}_ecoli_Aligned.sortedByCoord.out.bam" \
     "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  samtools index "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  log SPIKE "$n" done
done

###############################################################################
# 9  PICARD RG + DEDUP                                                        #
###############################################################################
picard_dedup () {                    # $1 = sample basename
  local n="$1"
  log Picard "$n" start

  local inbam="$ALIGNMENT_DIR/${n}.Aligned.sortedByCoord.out.bam"
  [[ -s $inbam ]] || { log Picard "$n" "skip (missing BAM)"; return; }

  "$PICARD_CMD" AddOrReplaceReadGroups \
       I="$inbam" \
       O="$ALIGNMENT_DIR/${n}.rg.bam" \
       RGID=1 RGLB=lib RGPL=ILM RGPU=unit RGSM="$n" \
       VALIDATION_STRINGENCY=LENIENT \
       >>"$LOG_DIR/picard_${n}.log" 2>&1

  "$PICARD_CMD" MarkDuplicates \
       I="$ALIGNMENT_DIR/${n}.rg.bam" \
       O="$ALIGNMENT_DIR/${n}.dedup.bam" \
       M="$LOG_DIR/${n}.metrics.txt" \
       REMOVE_DUPLICATES=true \
       VALIDATION_STRINGENCY=LENIENT \
       >>"$LOG_DIR/picard_${n}.log" 2>&1

  samtools index "$ALIGNMENT_DIR/${n}.dedup.bam" \
       >>"$LOG_DIR/picard_${n}.log" 2>&1

  log Picard "$n" done
}
export -f picard_dedup log

if command -v parallel >/dev/null 2>&1; then
  log Picard ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           picard_dedup ::: "${SAMPLES[@]}"
else
  log Picard ALL "GNU parallel not found – running serially"
  for n in "${SAMPLES[@]}"; do picard_dedup "$n"; done
fi

###############################################################################
# 10  FRAGMENT FILTER                                                         #
###############################################################################
# pick the awk condition once
case $FRAGMENT_SIZE_FILTER in
  histones)              AWK_CMD='{if($9>=130 && $9<=300 || $1~/^@/)print}';;
  transcription_factors) AWK_CMD='{if($9<130 || $1~/^@/)print}';;
  *)                     AWK_CMD='{if($9<1000||$1~/^@/)print}';;
esac

frag_filter () {                    # $1 = sample basename
  local n="$1"
  log FragFilt "$n" start

  samtools view -h "$ALIGNMENT_DIR/${n}.dedup.bam" \
    | awk "$AWK_CMD" \
    | samtools view -bS - \
        > "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" 2>>"$LOG_DIR/fragfilt_${n}.log"

  samtools index "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" \
        >>"$LOG_DIR/fragfilt_${n}.log" 2>&1

  log FragFilt "$n" done
}
export -f frag_filter log
export ALIGNMENT_DIR LOG_DIR NUM_THREADS AWK_CMD

if command -v parallel >/dev/null 2>&1; then
  log FragFilt ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           frag_filter ::: "${SAMPLES[@]}"
else
  log FragFilt ALL "GNU parallel not found – running serially"
  for n in "${SAMPLES[@]}"; do frag_filter "$n"; done
fi

###############################################################################
# 11  PRESEQ COMPLEXITY ESTIMATION                                            #
###############################################################################
preseq_one () {                      # $1 = sample name
  local n="$1"
  local bam="$ALIGNMENT_DIR/${n}.dedup.bam"
  local out="$OUTPUT_DIR/preseq/${n}_complexity.txt"

  log Preseq "$n" start

  if [[ ! -s $bam ]]; then
    log Preseq "$n" skip "missing BAM"
    return
  fi

  preseq lc_extrap -B -o "$out" -s 1000000 "$bam" \
    >>"$LOG_DIR/preseq_${n}.log" 2>&1 \
    && log Preseq "$n" done \
    || log Preseq "$n" FAIL
}

mkdir -p "$OUTPUT_DIR/preseq"

export -f preseq_one log
export ALIGNMENT_DIR OUTPUT_DIR LOG_DIR

if command -v parallel >/dev/null 2>&1; then
  log Preseq ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           preseq_one ::: "${SAMPLES[@]}"
else
  log Preseq ALL "GNU parallel not found – running serially"
  for n in "${SAMPLES[@]}"; do preseq_one "$n"; done
fi

###############################################################################
# 12  MERGE BAMs   (treatment & control groups)                               #
###############################################################################
# ── make merged BAMs *once*
T_MRG="$ALIGNMENT_DIR/treatment_merged.bam"
log MERGE Treatment "$T_MRG ${TREAT_NAMES[@]}"
merge_bams "$T_MRG" "${TREAT_NAMES[@]}"
log MERGE Treatment Done

if (( ${#CTRL_NAMES[@]} )); then
  CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"
  log MERGE Control "$CTRL_MRG ${CTRL_NAMES[@]}"
  merge_bams "$CTRL_MRG" "${CTRL_NAMES[@]}"
  log MERGE Control Done
fi

###############################################################################
# 13  MACS2 PEAKS: replicate, merged, pooled                                  #
###############################################################################
CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"
mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}
export PEAK_DIR ALIGNMENT_DIR GENOME_SIZE

# ── A  replicates (unchanged) ────────────────────────────────────────────────
for n in "${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}"; do
  BAM="$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
  macs2 callpeak -t "$BAM" \
       -f BAMPE -g "$GENOME_SIZE" \
       -n "$n" \
       --outdir "$PEAK_DIR/replicate" \
       2>&1 | tee -a "$LOG_DIR/macs2_replicate.log"
done

# ── B  treatment-merged vs control-merged ────────────────────────────────────
# Only if a control_merged.bam exists
if [[ -s $CTRL_MRG ]]; then
  MACS_CMDS+=(
    "macs2 callpeak \
        -t $CTRL_MRG \
        -f BAMPE -g $GENOME_SIZE \
        -n controlMerged \
        --outdir $PEAK_DIR/merged"
  )
fi

if [[ -s $T_MRG ]]; then
  if [[ -s $CTRL_MRG ]]; then
    MACS_CMDS+=("macs2 callpeak -t $T_MRG -c $CTRL_MRG -f BAMPE -g $GENOME_SIZE -n treatmentMerged_vs_controlMerged --outdir $PEAK_DIR/merged")
  else
    MACS_CMDS+=("macs2 callpeak -t $T_MRG -f BAMPE -g $GENOME_SIZE -n treatmentMerged --outdir $PEAK_DIR/merged")
  fi
else
  log MACS2merged ALL "skip (missing merged BAM)"
fi

# ── C  each replicate vs pooled control ──────────────────────────────────────
if [[ -s $CTRL_MRG ]]; then
  for n in "${TREAT_NAMES[@]}"; do
    MACS_CMDS+=("macs2 callpeak -t $ALIGNMENT_DIR/${n}.dedup.filtered.bam -c $CTRL_MRG -f BAMPE -g $GENOME_SIZE -n ${n}_vs_ctrlPooled --outdir $PEAK_DIR/pooled")
  done
fi

# ── run all commands with GNU parallel (or serial fallback) ──────────────────
if command -v parallel >/dev/null 2>&1; then
  log MACS2 ALL "running ${#MACS_CMDS[@]} jobs with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 ::: "${MACS_CMDS[@]}" \
    2>&1 | tee -a "$LOG_DIR/macs2_parallel.log"
else
  log MACS2 ALL "GNU parallel not found – running jobs serially"
  for cmd in "${MACS_CMDS[@]}"; do
    echo "[MACS2] $cmd" | tee -a "$LOG_DIR/macs2_serial.log"
    eval "$cmd" 2>&1 | tee -a "$LOG_DIR/macs2_serial.log"
  done
fi

###############################################################################
# 14  SPIKE SCALE FACTORS                                                     #
###############################################################################
declare -A SCALE                    # sample → host/spike factor
read_count () { samtools view -c -F 2304 "$1"; }

for n in "${SAMPLES[@]}"; do
  host_bam="$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
  spike_bam="$SPIKE_DIR/${n}.ecoli.sorted.bam"

  # ── skip if spike-in BAM is missing or empty ──────────────────────────────
  if [[ ! -s $spike_bam ]]; then
      log Scale "$n" "no spike-in BAM – scale=1"
      SCALE["$n"]=1
      continue
  fi
  # -------------------------------------------------------------------------

  h=$(read_count "$host_bam")
  s=$(read_count "$spike_bam")

  if (( s > 0 )); then
      SCALE["$n"]=$(awk -v h=$h -v s=$s 'BEGIN{printf "%.6f", h/s}')
  else
      SCALE["$n"]=1
  fi

  log Scale "$n" "host=$h  spike=$s  scale=${SCALE[$n]}"
done

###############################################################################
# 15  BIGWIG GENERATION                                                       #
###############################################################################
bigwig_one () {                       # $1 = sample basename
  local n="$1"
  local bg="$BW_DIR/${n}.bedgraph"
  local bam="$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
  local scalef="${SCALE[$n]:-1}"

  log BigWig "$n" "start  scale=$scalef"

  # BedGraph
  bedtools genomecov -ibam "$bam" -bg -pc > "$bg"

  # optional spike-in scaling
  if [[ $scalef != 1 ]]; then
      awk -v f="$scalef" '{$4=$4*f; print }' "$bg" > "${bg}.tmp" && mv "${bg}.tmp" "$bg"
  fi

  # BigWig
  bedGraphToBigWig "$bg" "$CHROM_SIZE" "$BW_DIR/${n}.bw"

  log BigWig "$n" done
}
export -f bigwig_one log
export BW_DIR ALIGNMENT_DIR CHROM_SIZE SCALE

if command -v parallel >/dev/null 2>&1; then
  log BigWig ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           bigwig_one ::: "${SAMPLES[@]}"
else
  log BigWig ALL "GNU parallel not found – running serially"
  for n in "${SAMPLES[@]}"; do bigwig_one "$n"; done
fi

###############################################################################
# 16  PEAK ANNOTATION                                                         #
###############################################################################
annotate_one () {                     # $1 = full path to narrowPeak
  local np="$1"
  [[ -f $np ]] || { log Annotating "$(basename "$np")" "skip (no file)"; return; }

  local base=${np##*/}; base=${base%.narrowPeak}
  log Annotate "$base" start

  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb \
      > "$ANN_DIR/${base}.annotated.bed"

  log Annotate "$base" done
}

export -f annotate_one log
export ANNOTATION_GENES ANN_DIR

# Collect all narrowPeak paths into an array
mapfile -t NP_FILES < <(find "$PEAK_DIR" -type f -name "*.narrowPeak" | sort)

if [[ ${#NP_FILES[@]} -eq 0 ]]; then
  log Annotate ALL "No peaks found – skipping annotation step"
else
  if command -v parallel >/dev/null 2>&1; then
    log Annotate ALL "running ${#NP_FILES[@]} peaks with GNU parallel (-j $NUM_PARALLEL_THREADS)"
    parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
             annotate_one ::: "${NP_FILES[@]}"
  else
    log Annotate ALL "GNU parallel not found – running serially"
    for np in "${NP_FILES[@]}"; do annotate_one "$np"; done
  fi
fi

###############################################################################
# 17  PRESEQ PLOTTING                                                         #
###############################################################################
plot_preseq_curves () {
  local preseq_dir="$OUTPUT_DIR/preseq"
  local plot_pdf="$OUTPUT_DIR/preseq/preseq_complexity_curves.pdf"
  local r_script="$OUTPUT_DIR/preseq/plot_preseq.R"

  log Preseq Plotting "start"

# Generate R script
export OUTPUT_DIR="$OUTPUT_DIR" # So the R script can access Sys.getenv("OUTPUT_DIR")
cat > "$r_script" <<'EOF'
# Load required libraries
library(ggplot2)
library(data.table)

# Function to read and validate a single preseq complexity file
plot_file <- function(file) {
  # Attempt to read the file, skipping the header line
  df <- tryCatch({
    fread(file, skip = 1)
  }, error = function(e) {
    stop(paste("ERROR reading", file, ":", e$message))
  })

  # Ensure the file has exactly 4 columns
  if (ncol(df) != 4) {
    stop(paste("ERROR:", file, "does not have exactly 4 columns (has", ncol(df), ")"))
  }

  # Assign standard column names
  setnames(df, c("total_reads", "expected_unique", "ci_lower", "ci_upper"))

  # Extract sample name from file name
  df[, sample := gsub("_complexity\\.txt$", "", basename(file))]

  return(df)
}

# Define the directory containing the preseq output files
preseq_dir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "."), "preseq")

# Find all complexity output files
files <- list.files(preseq_dir, pattern = "_complexity\\.txt$", full.names = TRUE)

# Stop if no files found
if (length(files) == 0) {
  stop("No _complexity.txt files found in 'preseq' directory: ", preseq_dir)
}

# Read and combine all complexity files into one data table
all_data <- rbindlist(lapply(files, plot_file), use.names = TRUE, fill = TRUE)

# Generate the plot
p <- ggplot(all_data, aes(x = total_reads, y = expected_unique, color = sample)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Preseq Library Complexity Projection",
    x = "Total Reads Sequenced",
    y = "Expected Unique Reads"
  ) +
  theme(legend.title = element_blank())

# Save the plot as PDF and PNG
ggsave(filename = file.path(preseq_dir, "preseq_complexity_curves.pdf"), plot = p, width = 8, height = 6)
ggsave(filename = file.path(preseq_dir, "preseq_complexity_curves.png"), plot = p, width = 8, height = 6, dpi = 300)
EOF

  # Run the R script
  Rscript "$r_script" >> "$LOG_DIR/preseq_plot.log" 2>&1 \
  && log Preseq Plotting "done: $plot_pdf" \
  || {
    log Preseq Plotting "FAIL – see preseq_plot.log"
    cat "$LOG_DIR/preseq_plot.log"
  }
}

plot_preseq_curves

# PIPELINE COMPLETED ##########################################################
log DONE ALL "Outputs in $OUTPUT_DIR"
