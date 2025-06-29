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
     seacr_peaks/
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
  7  SEACR  : replicate, merged, pooled peak calling (paired-end mode)
  8  Spike-in scale factors – host/spike read ratio per replicate
  9  BedGraph + BigWig generation (scaled if factors exist)
 10  Peak-to-gene annotation (bedtools intersect)

 USAGE
 ─────
   • Edit the config.json with absolute paths and replicate lists.
   • Run: bash seacr_cut_and_run_pipeline_w_replicates.sh (no CLI arguments)
   • Logs stream to stdout *and* to output/logs/pipeline.log

 REQUIREMENTS
 ────────────
   bash ≥4 · samtools ≥1.10 · bedtools ≥2.28 · STAR ≥2.7 · Java ≥17
   Trim Galore ≥0.6.10 · Picard ≥2.18 · SEACR ≥1.3 · cutadapt ≥4.1
   R 4.x (for optional DiffBind) · GNU coreutils/awk · PIGZ · GNU Parallel

BANNER

###############################################################################
# 0  CONFIG + PATHS                                                           #
###############################################################################
CONFIG_FILE="${1:-/mnt/data/home/aviv/cut_and_run/dockerfiles/config_for_docker.json}"
shift || true          # remove $1 so "$@" stays clean for future use

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
PICARD_CMD="$(command -v picard)"
export PICARD_CMD ALIGNMENT_DIR LOG_DIR NUM_THREADS
FASTQC_BIN="$(command -v fastqc)"
STAR_BIN="$(command -v STAR)"

###############################################################################
# 2  FOLDERS + LOGGER                                                         #
###############################################################################
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR"

# timestamp the old log (if any) before we start the new one
RUN_TS=$(date +'%Y%m%d_%H%M%S')                 # e.g. 20250624_104832
if [[ -f $LOG_DIR/pipeline.log ]]; then
    mv "$LOG_DIR/pipeline.log" \
       "$LOG_DIR/pipeline_${RUN_TS}.log"
fi
# create a fresh empty logfile for this run
: > "$LOG_DIR/pipeline.log"
                               
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports"
SPIKE_DIR="$ALIGNMENT_DIR/spikein"
PEAK_DIR="$OUTPUT_DIR/seacr_peaks"
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

# In the main script call it conditionally
if [[ $RUN_ONLY_BLOCK == "yes" ]]; then
    run_block
    exit 0
fi

#~ ###############################################################################
#~ # 5  FASTQC  – per-sample logging (start / ok / FAIL)                         #
#~ ###############################################################################
#~ log FastQC Conditions "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}"

#~ ALL_FASTQS=( "${TREAT_R1[@]}" "${TREAT_R2[@]}"
             #~ "${CTRL_R1[@]:-}" "${CTRL_R2[@]:-}" )

#~ run_fastqc () {                       # $1 = FASTQ
  #~ local fq="$1" base
  #~ base=$(basename "$fq")
  #~ log FastQC "$base" start
  #~ if "$FASTQC_BIN" --extract -o "$FASTQC_DIR" "$fq" \
        #~ >>"$LOG_DIR/fastqc_${base}.log" 2>&1 ; then
      #~ log FastQC "$base" done
  #~ else
      #~ log FastQC "$base" FAIL
  #~ fi
#~ }
#~ export -f run_fastqc log
#~ export LOG_DIR FASTQC_DIR FASTQC_BIN

#~ if command -v parallel >/dev/null 2>&1; then
  #~ log FastQC ALL "running ${#ALL_FASTQS[@]} files with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  #~ parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           #~ run_fastqc ::: "${ALL_FASTQS[@]}" \
    #~ 2>&1 | tee -a "$LOG_DIR/fastqc_parallel.log"
#~ else
  #~ log FastQC ALL "GNU parallel not found – running serially"
  #~ for fq in "${ALL_FASTQS[@]}"; do run_fastqc "$fq"; done
#~ fi

#~ ###############################################################################
#~ # 6  TRIMMING  – per-sample logging (start / done)                            #
#~ ###############################################################################
#~ log Trim ALL "Treatment=${#TREAT_R1[@]}  Control=${#CTRL_R1[@]}  (Trim Galore! --cores $NUM_PARALLEL_THREADS)"

#~ # ── treatment replicates ────────────────────────────────────────────────────
#~ for i in "${!TREAT_R1[@]}"; do
  #~ trim_one_pair "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}" || continue
#~ done

#~ # ── control replicates (if any) ─────────────────────────────────────────────
#~ for i in "${!CTRL_R1[@]}"; do
  #~ trim_one_pair "${CTRL_R1[$i]}" "${CTRL_R2[$i]}" "${CTRL_NAMES[$i]}" || continue
#~ done

#~ ###############################################################################
#~ # 7  HOST GENOME ALIGNMENT                                                    #
#~ ###############################################################################
#~ for n in "${SAMPLES[@]}"; do
  #~ log STARhost "$n" start
  #~ run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" \
           #~ "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           #~ "$ALIGNMENT_DIR/${n}." "$REFERENCE_GENOME" \
           #~ --outReadsUnmapped Fastx          # ← NEW OPTION
  #~ log STARhost "$n" done
#~ done

#~ ###############################################################################
#~ # 8  SPIKE-IN ALIGNMENT (E. coli)                                             #
#~ ###############################################################################
#~ for n in "${SAMPLES[@]}"; do
  #~ # locate mate FASTQs produced in Step 8
  #~ read -r U1 U2 < <(get_unmapped_mates "$n")

  #~ if [[ -z $U1 ]]; then
    #~ log SPIKE "$n" "skip (no unmapped mates)"
    #~ continue
  #~ fi

  #~ log SPIKE "$n" start
  #~ run_star "$U1" "$U2" \
           #~ "$SPIKE_DIR/${n}_ecoli_" "$ECOLI_INDEX"
  #~ mv "$SPIKE_DIR/${n}_ecoli_Aligned.sortedByCoord.out.bam" \
     #~ "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  #~ samtools index "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  #~ log SPIKE "$n" done
#~ done

#~ ###############################################################################
#~ # 9  PICARD RG + DEDUP                                                        #
#~ ###############################################################################
#~ picard_dedup () {                    # $1 = sample basename
  #~ local n="$1"
  #~ log Picard "$n" start

  #~ local inbam="$ALIGNMENT_DIR/${n}.Aligned.sortedByCoord.out.bam"
  #~ [[ -s $inbam ]] || { log Picard "$n" "skip (missing BAM)"; return; }

  #~ "$PICARD_CMD" AddOrReplaceReadGroups \
       #~ I="$inbam" \
       #~ O="$ALIGNMENT_DIR/${n}.rg.bam" \
       #~ RGID=1 RGLB=lib RGPL=ILM RGPU=unit RGSM="$n" \
       #~ VALIDATION_STRINGENCY=LENIENT \
       #~ >>"$LOG_DIR/picard_${n}.log" 2>&1

  #~ "$PICARD_CMD" MarkDuplicates \
       #~ I="$ALIGNMENT_DIR/${n}.rg.bam" \
       #~ O="$ALIGNMENT_DIR/${n}.dedup.bam" \
       #~ M="$LOG_DIR/${n}.metrics.txt" \
       #~ REMOVE_DUPLICATES=true \
       #~ VALIDATION_STRINGENCY=LENIENT \
       #~ >>"$LOG_DIR/picard_${n}.log" 2>&1

  #~ samtools index "$ALIGNMENT_DIR/${n}.dedup.bam" \
       #~ >>"$LOG_DIR/picard_${n}.log" 2>&1

  #~ log Picard "$n" done
#~ }
#~ export -f picard_dedup log

#~ if command -v parallel >/dev/null 2>&1; then
  #~ log Picard ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  #~ parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           #~ picard_dedup ::: "${SAMPLES[@]}"
#~ else
  #~ log Picard ALL "GNU parallel not found – running serially"
  #~ for n in "${SAMPLES[@]}"; do picard_dedup "$n"; done
#~ fi

###############################################################################
# 10  FRAGMENT FILTER                                                         #
###############################################################################
# pick the awk condition once
#~ case $FRAGMENT_SIZE_FILTER in
  #~ histones)              AWK_CMD='{if($9>=130 && $9<=300 || $1~/^@/)print}';;
  #~ transcription_factors) AWK_CMD='{if($9<130 || $1~/^@/)print}';;
  #~ *)                     AWK_CMD='{if($9<1000||$1~/^@/)print}';;
#~ esac

#~ frag_filter () {                    # $1 = sample basename
  #~ local n="$1"
  #~ log FragFilt "$n" start

  #~ samtools view -h "$ALIGNMENT_DIR/${n}.dedup.bam" \
    #~ | awk "$AWK_CMD" \
    #~ | samtools view -bS - \
        #~ > "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" 2>>"$LOG_DIR/fragfilt_${n}.log"

  #~ samtools index "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" \
        #~ >>"$LOG_DIR/fragfilt_${n}.log" 2>&1

  #~ log FragFilt "$n" done
#~ }
#~ export -f frag_filter log
#~ export ALIGNMENT_DIR LOG_DIR NUM_THREADS AWK_CMD

#~ if command -v parallel >/dev/null 2>&1; then
  #~ log FragFilt ALL "running ${#SAMPLES[@]} samples with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  #~ parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           #~ frag_filter ::: "${SAMPLES[@]}"
#~ else
  #~ log FragFilt ALL "GNU parallel not found – running serially"
  #~ for n in "${SAMPLES[@]}"; do frag_filter "$n"; done
#~ fi

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
  # ── NEW: make the bedGraph SEACR expects ────────────────────────────────
  mkdir -p "$BW_DIR"
  bedtools genomecov -bg -ibam "$CTRL_MRG" \
         > "$BW_DIR/control_merged.bedgraph"
fi

CTRL_BG="$BW_DIR/control_merged.bedgraph"
if [[ ! -s $CTRL_BG ]]; then
    mkdir -p "$BW_DIR"
    bedtools genomecov -bg -ibam "$CTRL_MRG" > "$CTRL_BG"
fi

###############################################################################
# 12  SEACR PEAKS: replicate, merged, pooled                                  #
###############################################################################
#  Make ONE self-contained, writable SEACR copy and shadow all others
export BW_DIR PEAK_DIR GENOME_SIZE LOG_DIR

mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}

# Writable home for wrapper + helper
SEACR_HOME=/tmp/seacr_bin
mkdir -p  "$SEACR_HOME"
chmod 1777 "$SEACR_HOME"               # world-writable, sticky

# Copy BOTH required files
cp  "$(command -v seacr)"                 "$SEACR_HOME/seacr"
cp  "$(dirname "$(command -v seacr)")/SEACR_1.3.R"  "$SEACR_HOME/"
chmod +x "$SEACR_HOME/seacr"

# Shadow every earlier copy *for this script and all subshells*
export PATH="$SEACR_HOME:$PATH"

# Helper: always run from SEACR_HOME so mktemp files land in a writable place
seacr_call () (
    cd "$SEACR_HOME"         # change cwd only for this subshell
    seacr "$@"               # all arguments stay absolute
)

echo "Using seacr at: $(command -v seacr)"
ls -l "$SEACR_HOME"

# ── A  replicate peaks ───────────────────────────────────────────────────────
for n in "${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}"; do
  IN_BG="$BW_DIR/${n}.bedgraph"                    # produced in step 14
  OUT_BED="$PEAK_DIR/replicate/${n}_seacr.bed"

  log SEACR "$n" start
  seacr_call "$IN_BG" 0.01 non stringent "$OUT_BED" \
      >>"$LOG_DIR/seacr_${n}.log" 2>&1
  if [[ $? -eq 0 && -s $OUT_BED ]]; then
      log SEACR "$n" done
  else
      log SEACR "$n" FAIL
  fi
done

# ── B  treatment-merged vs control-merged ────────────────────────────────────
if [[ -s $T_MRG && -s $CTRL_MRG ]]; then
  MERGED_T_BG="$BW_DIR/treatment_merged.bedgraph"
  MERGED_C_BG="$BW_DIR/control_merged.bedgraph"
  OUT_BED="$PEAK_DIR/merged/treatmentMerged_vs_controlMerged_seacr.bed"

  log SEACR treatmentMerged_vs_controlMerged start
  seacr_call "$MERGED_T_BG" "$MERGED_C_BG" non stringent "$OUT_BED" \
      >>"$LOG_DIR/seacr_merged.log" 2>&1 \
      && log SEACR treatmentMerged_vs_controlMerged done \
      || log SEACR treatmentMerged_vs_controlMerged FAIL
fi

# ── C  each replicate vs pooled control ──────────────────────────────────────
if [[ -s $CTRL_MRG ]]; then
  POOLED_C_BG="$BW_DIR/control_merged.bedgraph"
  for n in "${TREAT_NAMES[@]}"; do
    IN_BG="$BW_DIR/${n}.bedgraph"
    OUT_BED="$PEAK_DIR/pooled/${n}_vs_ctrlPooled_seacr.bed"

    log SEACR "${n}_vs_ctrlPooled" start
    seacr_call "$IN_BG" "$POOLED_C_BG" non stringent "$OUT_BED" \
        >>"$LOG_DIR/seacr_${n}.log" 2>&1 \
        && log SEACR "${n}_vs_ctrlPooled" done \
        || log SEACR "${n}_vs_ctrlPooled" FAIL
  done
fi

###############################################################################
# 13  SPIKE SCALE FACTORS                                                     #
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
# 14  BIGWIG GENERATION                                                       #
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
# 15  PEAK ANNOTATION                                                         #
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

if command -v parallel >/dev/null 2>&1; then
  log Annotate ALL "running ${#NP_FILES[@]} peaks with GNU parallel (-j $NUM_PARALLEL_THREADS)"
  parallel --line-buffer -j "$NUM_PARALLEL_THREADS" --halt now,fail=1 \
           annotate_one ::: "${NP_FILES[@]}"
else
  log Annotate ALL "GNU parallel not found – running serially"
  for np in "${NP_FILES[@]}"; do annotate_one "$np"; done
fi

# PIPELINE COMPLETED ##########################################################
log DONE ALL "Outputs in $OUTPUT_DIR"
