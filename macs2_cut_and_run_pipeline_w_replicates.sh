#!/bin/bash
set -o pipefail

cat <<'BANNER'

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ CUT&RUN PIPELINE (Paired-End) -- Replicates · E. coli Spike-in · MACS2 Peaks ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
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
  2  Trimming (Trimmomatic)   – all treatment & control replicates
  3  STAR   : E. coli (spike-in)    → per-replicate BAM + index
  4  STAR   : host genome (hg38/mm10/…) → per-replicate BAM + index
  5  Picard : add-RG + duplicate removal
  6  Fragment-size filtering (histone/TF/≤1 kb)
  7  MACS2  : replicate, merged, pooled peak calling (BAMPE mode)
  8  Spike-in scale factors   – host/spike read ratio per replicate
  9  BedGraph + BigWig generation (scaled if factors exist)
 10  Peak-to-gene annotation (bedtools intersect)
 11  DiffBind peak comparison

 USAGE
 ─────
   • Edit config.json with absolute paths and replicate lists.
   • Run:   bash cutrun_pipeline.sh   (no CLI arguments)
   • Logs stream to stdout *and* to output/logs/pipeline.log

 REQUIREMENTS
 ────────────
   bash ≥4   ·  samtools ≥1.10   ·  bedtools ≥2.28   ·  STAR ≥2.7
   Trim Galore 0.6.10 ·  Picard ≥2.18 ·  MACS2 ≥2.2 · cutadapt ≥4.1
   R 4.x (for optional DiffBind) · GNU coreutils/awk

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

log START ALL "Config=$CONFIG_FILE"

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
              --cores "$NUM_THREADS" \
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

run_fastqc () {                     # $1=R1  $2=R2  $3=SAMPLE
  local R1=$1 R2=$2 SAMPLE=$3
  log FastQC "$SAMPLE" start
  $FASTQC_BIN --quiet --extract -o "$FASTQC_DIR" "$R1" "$R2" \
      >  "$LOG_DIR/fastqc_${SAMPLE}.log" \
      2> "$LOG_DIR/fastqc_${SAMPLE}.err"
  if [[ $? -eq 0 ]]; then
      log FastQC "$SAMPLE" ok
  else
      log FastQC "$SAMPLE" FAIL  "see fastqc_${SAMPLE}.err"
  fi
}

bam_to_bedgraph(){ bedtools genomecov -ibam "$1" -bg -pc | sort -k1,1 -k2,2n > "$2"; }
merge_bams(){ [[ -f "$1" ]] || { samtools merge -f "$1" "${@:2}" && samtools index "$1"; }; }
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
#~ log FastQC ALL "T=${#TREAT_R1[@]}  C=${#CTRL_R1[@]}"

#~ for i in "${!TREAT_R1[@]}"; do
  #~ run_fastqc "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}"
#~ done
#~ for i in "${!CTRL_R1[@]}";  do
  #~ run_fastqc "${CTRL_R1[$i]}"  "${CTRL_R2[$i]}"  "${CTRL_NAMES[$i]}"
#~ done

###############################################################################
# 6  TRIMMING  – per-sample logging (start / done)                            #
###############################################################################
# Headline
log Trim ALL "T=${#TREAT_R1[@]}  C=${#CTRL_R1[@]}  (Trim Galore! --cores $NUM_THREADS)"

# ── treatment replicates ────────────────────────────────────────────────────
for i in "${!TREAT_R1[@]}"; do
  echo "DEBUG ${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}"
  trim_one_pair "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}" || continue
done

# ── control replicates (if any) ─────────────────────────────────────────────
for i in "${!CTRL_R1[@]}"; do
  echo "DEBUG ${CTRL_R1[$i]}" "${CTRL_R2[$i]}" "${CTRL_NAMES[$i]}"
  trim_one_pair "${CTRL_R1[$i]}" "${CTRL_R2[$i]}" "${CTRL_NAMES[$i]}" || continue
done

###############################################################################
# 7  SPIKE-IN ALIGNMENT (E. coli)                                             #
###############################################################################
for n in "${SAMPLES[@]}"; do
  log SPIKE "$n" start
  run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" \
           "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           "$SPIKE_DIR/${n}_ecoli_"  "$ECOLI_INDEX"

  mv "$SPIKE_DIR/${n}_ecoli_Aligned.sortedByCoord.out.bam" \
     "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  samtools index "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  log SPIKE    "$n" done
done

###############################################################################
# 8  HOST GENOME ALIGNMENT                                                    #
###############################################################################
for n in "${SAMPLES[@]}"; do
  log STARhost "$n" start
  run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" \
           "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           "$ALIGNMENT_DIR/${n}." "$REFERENCE_GENOME"
  log STARhost "$n" done

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
# 11  MACS2 PEAKS                                                             #
###############################################################################
mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}
merge_bams "$ALIGNMENT_DIR/control_merged.bam" "${CTRL_NAMES[@]/%/.dedup.filtered.bam}" 2>/dev/null
CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"

# replicate peaks
for n in "${TREAT_NAMES[@]}"; do
  log MACS2rep "$n"
  if [[ -f $CTRL_MRG ]]; then
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" -c "$CTRL_MRG" \
           -f BAMPE -g "$GENOME_SIZE" -n "$n" --outdir "$PEAK_DIR/replicate"
  else
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" \
           -f BAMPE -g "$GENOME_SIZE" -n "$n" --outdir "$PEAK_DIR/replicate"
  fi
done

# merged peaks
merge_bams "$ALIGNMENT_DIR/treat_merged.bam" "${TREAT_NAMES[@]/%/.dedup.filtered.bam}"
if [[ -f $CTRL_MRG ]]; then
  log MACS2mg Merged
  $MACS2 callpeak -t "$ALIGNMENT_DIR/treat_merged.bam" -c "$CTRL_MRG" \
         -f BAMPE -g "$GENOME_SIZE" -n treatMerged_vs_ctrlMerged --outdir "$PEAK_DIR/merged"
else
  log MACS2mg Merged
  $MACS2 callpeak -t "$ALIGNMENT_DIR/treat_merged.bam" \
         -f BAMPE -g "$GENOME_SIZE" -n treatMerged --outdir "$PEAK_DIR/merged"
fi

# pooled control vs each replicate
if [[ -f $CTRL_MRG ]]; then
  for n in "${TREAT_NAMES[@]}"; do
    log MACS2pl "$n"
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" -c "$CTRL_MRG" \
           -f BAMPE -g "$GENOME_SIZE" -n "${n}_vs_ctrlPooled" --outdir "$PEAK_DIR/pooled"
  done
fi

###############################################################################
# 12  SPIKE SCALE FACTORS                                                     #
###############################################################################
declare -A SCALE
for n in "${SAMPLES[@]}"; do
  h=$(read_count "$ALIGNMENT_DIR/${n}.dedup.filtered.bam")
  s=$(read_count "$SPIKE_DIR/${n}.ecoli.sorted.bam")
  ((s)) && SCALE["$n"]=$(awk -v h=$h -v s=$s 'BEGIN{printf "%.6f",h/s}')
  log Scale "$n" "host=$h spike=$s scale=${SCALE[$n]:-NA}"
done

###############################################################################
# 13  BIGWIG                                                                  #
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
# 14  PEAK ANNOTATION                                                         #
###############################################################################
for np in "$PEAK_DIR"/{replicate,merged,pooled}/*.narrowPeak; do
  [[ -f $np ]] || continue
  b=${np##*/}; b=${b%.narrowPeak}
  log Annotate "$b"
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb \
    > "$ANN_DIR/${b}.annotated.bed"
done

###############################################################################
# 15  DIFFBIND  – merged treatment vs merged control peaks                    #
###############################################################################
DIFF_DIR="$OUTPUT_DIR/diffbind"
mkdir -p "$DIFF_DIR"

# ── skip if no control or fewer than 2 treatment replicates ───────────
if [[ ${#CTRL_NAMES[@]} -eq 0 || ${#TREAT_NAMES[@]} -lt 2 ]]; then
  log DiffBind ALL "skipped (need ≥2 treatment reps + control)"
else
  log DiffBind ALL start

  # build sample sheet
  SAMPLE_SHEET="$DIFF_DIR/diffbind_samples.csv"
  {
    echo "SampleID,Condition,bamReads,Peaks,ScoreCol"
    echo "TreatmentMerged,treatment,$ALIGNMENT_DIR/treat_merged.bam,$PEAK_DIR/merged/treatMerged_vs_ctrlMerged_peaks.narrowPeak,7"
    echo "ControlMerged,control,$ALIGNMENT_DIR/control_merged.bam,$PEAK_DIR/merged/treatMerged_vs_ctrlMerged_control.narrowPeak,7"
  } > "$SAMPLE_SHEET"

  # run R script
  Rscript /mnt/data/home/aviv/cut_and_run/diffbind.R \
          "$SAMPLE_SHEET" "$DIFF_DIR" \
          > "$DIFF_DIR/diffbind.log" 2>&1

  if [[ $? -eq 0 ]]; then
    log DiffBind ALL ok
  else
    log DiffBind ALL FAIL "see $DIFF_DIR/diffbind.log"
  fi
fi

log DONE ALL "Outputs in $OUTPUT_DIR"
