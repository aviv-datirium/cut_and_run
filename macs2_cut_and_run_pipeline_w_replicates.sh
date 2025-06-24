#!/bin/bash

set -o pipefail

cat <<'BANNER'

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ CUT&RUN PIPELINE (Paired-End) - Replicates · E. coli Spike-in · MACS2 Peaks ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
 This Bash workflow trims, aligns, deduplicates, filters, and peak-calls
 replicated CUT&RUN libraries.  It supports an optional IgG/empty-vector control
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
        replicate//*.narrowPeak
        pooled//*.narrowPeak
        merged/…/*.narrowPeak
     annotated_peaks/          ⇠   *.annotated.{bed,tsv}
     logs/                     ⇠   STAR, Picard, Trim, etc.

 MAJOR STEPS
 ───────────
  0  Read paths & parameters (config.json)
  1  FastQC on raw FASTQs
  2  Trimming (Trimmomatic) – all treatment & control replicates
  3  STAR   : host genome (hg38/mm10/…) → per-replicate BAM + index
  4  STAR   : E. coli (spike-in) → per-replicate BAM + index
  5  Picard : add-RG + duplicate removal
  6  Fragment-size filtering (histone/TF/≤1 kb)
  7  MACS2  : replicate, merged, pooled peak calling (BAMPE mode)
  8  Spike-in scale factors – host/spike read ratio per replicate
  9  BedGraph + BigWig generation (scaled if factors exist)
 10  Peak-to-gene annotation (bedtools intersect)
 11  DiffBind peak comparison

 USAGE
 ─────
   • Edit config_replicates_diffbind.json with absolute paths and replicate lists.
   • Run: bash macs2_cut_and_run_pipeline_w_replicates.sh (no CLI arguments)
   • Logs stream to stdout *and* to output/logs/pipeline.log

 REQUIREMENTS
 ────────────
   bash ≥4 · samtools ≥1.10 · bedtools ≥2.28 · STAR ≥2.7 · Java ≥17
   Trim Galore 0.6.10 · Picard ≥2.18 · MACS2 ≥2.2 · cutadapt ≥4.1
   R 4.x (for optional DiffBind) · GNU coreutils/awk · PIGZ · GNU Parallel

BANNER

###############################################################################
# 0  CONFIG + PATHS                                                           #
###############################################################################
CONFIG_FILE="/mnt/data/home/aviv/cut_and_run/config_replicates_diffbind.json"

ALIGNMENT_DIR=$(jq -r '.alignment_dir'  "$CONFIG_FILE")
OUTPUT_DIR=$(   jq -r '.output_dir'     "$CONFIG_FILE")
LOG_DIR=$(      jq -r '.log_dir'        "$CONFIG_FILE")

REFERENCE_GENOME=$(jq -r '.reference_genome'   "$CONFIG_FILE")
ECOLI_INDEX=$(     jq -r '.ecoli_index'        "$CONFIG_FILE")
ANNOTATION_GENES=$(jq -r '.annotation_genes'   "$CONFIG_FILE")
CHROM_SIZE=$(      jq -r '.chrom_sizes'        "$CONFIG_FILE")

GENOME_SIZE_STRING=$(jq -r '.genome_size'      "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' "$CONFIG_FILE")
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
NUM_THREADS=$(jq    -r '.num_threads'          "$CONFIG_FILE")
NUM_PARALLEL_THREADS=$(jq    -r '.num_parallel_threads'          "$CONFIG_FILE")

TREAT_R1=($(jq -r '.samples.treatment[]?.r1' "$CONFIG_FILE"))
TREAT_R2=($(jq -r '.samples.treatment[]?.r2' "$CONFIG_FILE"))
CTRL_R1=($( jq -r '.samples.control[]?.r1 // empty' "$CONFIG_FILE"))
CTRL_R2=($( jq -r '.samples.control[]?.r2 // empty' "$CONFIG_FILE"))

###############################################################################
# 1  TOOL LOCATIONS                                                           #
###############################################################################
PICARD_JAR="/mnt/data/home/aviv/tools/picard.jar"
FASTQC_BIN="/mnt/data/home/aviv/tools/FastQC/fastqc"
STAR_BIN="/mnt/data/home/aviv/tools/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR"
TRIMMOMATIC_JAR="/mnt/data/home/aviv/tools/Trimmomatic/trimmomatic-0.39.jar"
ADAPTER_FA="/mnt/data/home/aviv/tools/Trimmomatic/adapters/TruSeq3-PE.fa"
MACS2="macs2"

###############################################################################
# 2  FOLDERS + LOGGER                                                         #
###############################################################################
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports"
SPIKE_DIR="$ALIGNMENT_DIR/spikein"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks"
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

###############################################################################
# merge_bams  <out.bam>  <sampleName> ...                                     #
#   • 0 inputs  → skip with log                                               #
#   • 1 input   → ln -f (fast)                                                #
#   • ≥2 inputs → samtools merge + index                                      #
###############################################################################
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
# 4  BASENAMES + GENOME SIZE                                                  #
###############################################################################
TREAT_NAMES=(); for r1 in "${TREAT_R1[@]}"; do TREAT_NAMES+=( "$(get_sample_basename "$r1")" ); done
CTRL_NAMES=();  for r1 in "${CTRL_R1[@]}";  do CTRL_NAMES+=( "$(get_sample_basename "$r1")" ); done
SAMPLES=("${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}")

case $GENOME_SIZE_STRING in
  hs) GENOME_SIZE=2913022398 ;;  mm) GENOME_SIZE=2652783500 ;;
  dm) GENOME_SIZE=165000000  ;;  ce) GENOME_SIZE=1000000000 ;;
  sc) GENOME_SIZE=12000000   ;;  *)  GENOME_SIZE=$CUSTOM_GENOME_SIZE ;;
esac

###############################################################################
# 5  FASTQC  – per-sample logging (start / ok / FAIL)                         #
###############################################################################
log FastQC Conditions "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}"
log FastQC ALL "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}"

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
for n in "${SAMPLES[@]}"; do
  log Picard "$n"
  in="$ALIGNMENT_DIR/${n}.Aligned.sortedByCoord.out.bam"
  [[ -s $in ]] || continue
  java -jar "$PICARD_JAR" AddOrReplaceReadGroups I="$in" O="$ALIGNMENT_DIR/${n}.rg.bam" \
       RGID=1 RGLB=lib RGPL=ILM RGPU=unit RGSM="$n" > /dev/null
  java -jar "$PICARD_JAR" MarkDuplicates I="$ALIGNMENT_DIR/${n}.rg.bam" \
       O="$ALIGNMENT_DIR/${n}.dedup.bam" M="$LOG_DIR/${n}.metrics.txt" REMOVE_DUPLICATES=true > /dev/null
  samtools index "$ALIGNMENT_DIR/${n}.dedup.bam"
done

###############################################################################
# 10  FRAGMENT FILTER                                                         #
###############################################################################
case $FRAGMENT_SIZE_FILTER in
  histones)              CMD='{if($9>=130&&$9<=300||$1~/^@/)print}';;
  transcription_factors) CMD='{if($9<130||$1~/^@/)print}';;
  *)                     CMD='{if($9<1000||$1~/^@/)print}';;
esac
for n in "${SAMPLES[@]}"; do
  log FragFilt "$n"
  samtools view -h "$ALIGNMENT_DIR/${n}.dedup.bam" | awk "$CMD" | \
    samtools view -bS - > "$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
  samtools index "$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
done

###############################################################################
# 11  MERGE BAMs   (treatment & control groups)                               #
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
# 12  MACS2 PEAKS: replicate, merged, pooled                                  #
###############################################################################
mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}
CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"

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
# 13  SPIKE SCALE FACTORS                                                     #
###############################################################################
declare -A SCALE
for n in "${SAMPLES[@]}"; do
  h=$(read_count "$ALIGNMENT_DIR/${n}.dedup.filtered.bam")
  s=$(read_count "$SPIKE_DIR/${n}.ecoli.sorted.bam")
  ((s)) && SCALE["$n"]=$(awk -v h=$h -v s=$s 'BEGIN{printf "%.6f",h/s}')
  log Scale "$n" "host=$h spike=$s scale=${SCALE[$n]:-NA}"
done

###############################################################################
# 14  BIGWIG                                                                  #
###############################################################################
for n in "${SAMPLES[@]}"; do
  log BigWig "$n" "scale=${SCALE[$n]:-1}"
  bg="$BW_DIR/${n}.bedgraph"
  bam_to_bedgraph "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" "$bg"
  [[ -n ${SCALE[$n]} ]] && \
      awk -v f="${SCALE[$n]}" '{$4=$4*f;print}' "$bg" > "${bg}.tmp" && mv "${bg}.tmp" "$bg"
  bedGraphToBigWig "$bg" "$CHROM_SIZE" "$BW_DIR/${n}.bw"
done

###############################################################################
# 15  PEAK ANNOTATION                                                         #
###############################################################################
for np in "$PEAK_DIR"/{replicate,merged,pooled}/*.narrowPeak; do
  [[ -f $np ]] || continue
  b=${np##*/}; b=${b%.narrowPeak}
  log Annotate "$b"
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb \
    > "$ANN_DIR/${b}.annotated.bed"
done

###############################################################################
# 16  DIFFBIND  – merged treatment vs merged control peaks                    #
###############################################################################
DIFF_DIR="$OUTPUT_DIR/diffbind"
mkdir -p "$DIFF_DIR"

T_N=${#TREAT_NAMES[@]}
C_N=${#CTRL_NAMES[@]}
if (( T_N < 2 || C_N < 2 )); then
  log DiffBind ALL "skip (need ≥2 replicates per group; T=$T_N C=$C_N)"
  exit 0            # top-level; use 'return 0' if inside a function
fi

SAMPLE_SHEET="$DIFF_DIR/diffbind_samples.csv"
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,ScoreCol" > "$SAMPLE_SHEET"

# treatment replicates
rep=1
for s in "${TREAT_NAMES[@]}"; do
  peaks="$PEAK_DIR/replicate/${s}_peaks.narrowPeak"
  [[ -s $peaks ]] || { log DiffBind "$s" "skip (no peaks)"; continue; }
  echo "${s},NA,NA,treatment,${rep},$ALIGNMENT_DIR/${s}.dedup.filtered.bam,$peaks,7" \
       >> "$SAMPLE_SHEET"
  ((rep++))
done

# control replicates
rep=1
for s in "${CTRL_NAMES[@]}"; do
  peaks="$PEAK_DIR/replicate/${s}_peaks.narrowPeak"
  [[ -s $peaks ]] || { log DiffBind "$s" "skip (no peaks)"; continue; }
  echo "${s},NA,NA,control,${rep},$ALIGNMENT_DIR/${s}.dedup.filtered.bam,$peaks,7" \
       >> "$SAMPLE_SHEET"
  ((rep++))
done

# after filtering rows
t_rows=$(grep -c ',treatment,' "$SAMPLE_SHEET")
c_rows=$(grep -c ',control,'   "$SAMPLE_SHEET")
if (( t_rows < 2 || c_rows < 2 )); then
  log DiffBind ALL "skip (after filtering, T=$t_rows C=$c_rows)"
  exit 0
fi

log DiffBind ALL start
Rscript /mnt/data/home/aviv/cut_and_run/diffbind.R "$SAMPLE_SHEET" "$DIFF_DIR" "$NUM_THREADS"\
       > "$DIFF_DIR/diffbind.log" 2>&1 \
  && log DiffBind ALL ok \
  || log DiffBind ALL FAIL

# PIPELINE COMPLETED ##########################################################
log DONE ALL "Outputs in $OUTPUT_DIR"
