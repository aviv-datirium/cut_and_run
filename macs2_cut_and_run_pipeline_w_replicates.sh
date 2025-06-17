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

 USAGE
 ─────
   • Edit config.json with absolute paths and replicate lists.
   • Run:   bash cutrun_pipeline.sh   (no CLI arguments)
   • Logs stream to stdout *and* to output/logs/pipeline.log

 REQUIREMENTS
 ────────────
   bash ≥4   ·  samtools ≥1.10   ·  bedtools ≥2.28   ·  STAR ≥2.7
   Trimmomatic 0.39   ·  Picard ≥2.18   ·  MACS2 ≥2.2
   R 4.x (for MultiQC/optional DiffBind)   ·  GNU coreutils / awk

BANNER

###############################################################################
# 0  LOAD FROM CONFIG                                                         #
###############################################################################
CONFIG_FILE="/mnt/data/home/aviv/cut_and_run/config_replicates.json"

RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir'  "$CONFIG_FILE")
ALIGNMENT_DIR=$(jq -r '.alignment_dir'  "$CONFIG_FILE")
OUTPUT_DIR=$(   jq -r '.output_dir'     "$CONFIG_FILE")
LOG_DIR=$(      jq -r '.log_dir'        "$CONFIG_FILE")

REFERENCE_GENOME=$(jq -r '.reference_genome' "$CONFIG_FILE")
ECOLI_INDEX=$(     jq -r '.ecoli_index'      "$CONFIG_FILE")
ANNOTATION_GENES=$(jq -r '.annotation_genes' "$CONFIG_FILE")
CHROM_SIZE=$(      jq -r '.chrom_sizes'      "$CONFIG_FILE")

GENOME_SIZE_STRING=$(jq -r '.genome_size'      "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' "$CONFIG_FILE")

FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
NUM_THREADS=$(jq    -r '.num_threads'          "$CONFIG_FILE")

TREAT_R1=($(jq -r '.samples.treatment[]?.r1' "$CONFIG_FILE"))
TREAT_R2=($(jq -r '.samples.treatment[]?.r2' "$CONFIG_FILE"))
CTRL_R1=($( jq -r '.samples.control[]?.r1 // empty' "$CONFIG_FILE"))
CTRL_R2=($( jq -r '.samples.control[]?.r2 // empty' "$CONFIG_FILE"))

###############################################################################
# 1  TOOL PATHS                                                               #
###############################################################################
PICARD_JAR="/mnt/data/home/aviv/tools/picard.jar"
FASTQC_BIN="/mnt/data/home/aviv/tools/FastQC/fastqc"
STAR_BIN="/mnt/data/home/aviv/tools/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR"
TRIMMOMATIC_JAR="/mnt/data/home/aviv/tools/Trimmomatic/trimmomatic-0.39.jar"
ADAPTER_FA="/mnt/data/home/aviv/tools/Trimmomatic/adapters/TruSeq3-PE.fa"
MACS2="macs2"

###############################################################################
# 2  DIRECTORIES + LOGGER                                                     #
###############################################################################
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports"
SPIKE_DIR="$ALIGNMENT_DIR/spikein"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks"
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs"
ANN_DIR="$OUTPUT_DIR/annotated_peaks"
for d in "$FASTQC_DIR" "$SPIKE_DIR" "$PEAK_DIR"/{replicate,merged,pooled} \
         "$BW_DIR" "$ANN_DIR"; do mkdir -p "$d"; done

log () {                                   # log "STEP" "detail ..."
  printf '[%(%F %T)T] %-10s %s\n' -1 "$1" "${*:2}" | tee -a "$LOG_DIR/pipeline.log"
}

log "START" "Config = $CONFIG_FILE"

###############################################################################
# 3  HELPER FUNCTIONS                                                         #
###############################################################################
get_sample_basename(){ local b=${1##*/}; echo "${b%%[_R12]*}"; }

trim_one_pair(){                        # $1 R1  $2 R2  $3 BASE
  log "Trim" "$3"
  java -jar "$TRIMMOMATIC_JAR" PE -threads "$NUM_THREADS" -phred33 \
       "$1" "$2" \
       "$ALIGNMENT_DIR/${3}_trimmed_R1.fq.gz" /dev/null \
       "$ALIGNMENT_DIR/${3}_trimmed_R2.fq.gz" /dev/null \
       ILLUMINACLIP:"$ADAPTER_FA":2:30:10 LEADING:5 TRAILING:5 \
       SLIDINGWINDOW:4:15 MINLEN:25 \
       > "$LOG_DIR/trim_${3}.log" 2>&1
}

run_star(){                              # $1 R1  $2 R2  $3 OUTPREFIX  $4 LOGBASE  $5 GENOME
  log "$4" "STAR align"
  "$STAR_BIN" --runThreadN "$NUM_THREADS" \
              --genomeDir "$5" \
              --readFilesIn "$1" "$2" \
              --readFilesCommand zcat \
              --outSAMtype BAM SortedByCoordinate \
              --outFileNamePrefix "$3" \
              > "$LOG_DIR/${4}.log" 2> "$LOG_DIR/${4}_err.log"
}

bam_to_bedgraph(){ bedtools genomecov -ibam "$1" -bg -pc | sort -k1,1 -k2,2n > "$2"; }

merge_bams(){ [[ -f "$1" ]] || { samtools merge -f "$1" "${@:2}" && samtools index "$1"; }; }

###############################################################################
# 4  BASENAMES & GENOME SIZE                                                  #
###############################################################################
TREAT_NAMES=(); for r1 in "${TREAT_R1[@]}"; do TREAT_NAMES+=( "$(get_sample_basename "$r1")" ); done
CTRL_NAMES=();  for r1 in "${CTRL_R1[@]}";  do CTRL_NAMES+=( "$(get_sample_basename "$r1")" ); done
SAMPLES=("${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}")

case $GENOME_SIZE_STRING in
  hs) GENOME_SIZE=2913022398 ;;
  mm) GENOME_SIZE=2652783500 ;;
  dm) GENOME_SIZE=165000000  ;;
  ce) GENOME_SIZE=1000000000 ;;
  sc) GENOME_SIZE=12000000   ;;
  *)  GENOME_SIZE=$CUSTOM_GENOME_SIZE ;;
esac

###############################################################################
# 5  FASTQC                                                                   #
###############################################################################
log "FastQC" "T=${#TREAT_R1[@]}  C=${#CTRL_R1[@]}"
for i in "${!TREAT_R1[@]}"; do $FASTQC_BIN --quiet -o "$FASTQC_DIR" "${TREAT_R1[$i]}" "${TREAT_R2[$i]}"  > /dev/null; done
for i in "${!CTRL_R1[@]}";  do $FASTQC_BIN --quiet -o "$FASTQC_DIR" "${CTRL_R1[$i]}"  "${CTRL_R2[$i]}"   > /dev/null; done

###############################################################################
# 6  TRIMMING                                                                 #
###############################################################################
for i in "${!TREAT_R1[@]}"; do trim_one_pair "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}"; done
for i in "${!CTRL_R1[@]}";  do trim_one_pair "${CTRL_R1[$i]}"  "${CTRL_R2[$i]}"  "${CTRL_NAMES[$i]}";  done

###############################################################################
# 7  SPIKE-IN ALIGNMENT                                                       #
###############################################################################
for n in "${SAMPLES[@]}"; do
  run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           "$SPIKE_DIR/${n}_ecoli_" "SPIKE_$n" "$ECOLI_INDEX"
  mv "$SPIKE_DIR/${n}_ecoli_Aligned.sortedByCoord.out.bam" "$SPIKE_DIR/${n}.ecoli.sorted.bam"
  samtools index "$SPIKE_DIR/${n}.ecoli.sorted.bam"
done

###############################################################################
# 8  HOST ALIGNMENT                                                           #
###############################################################################
for n in "${SAMPLES[@]}"; do
  run_star "$ALIGNMENT_DIR/${n}_trimmed_R1.fq.gz" "$ALIGNMENT_DIR/${n}_trimmed_R2.fq.gz" \
           "$ALIGNMENT_DIR/${n}." "STAR_$n" "$REFERENCE_GENOME"
done

###############################################################################
# 9  PICARD RG + DEDUP                                                        #
###############################################################################
for n in "${SAMPLES[@]}"; do
  in="$ALIGNMENT_DIR/${n}.Aligned.sortedByCoord.out.bam"
  [[ -s $in ]] || continue
  java -jar "$PICARD_JAR" AddOrReplaceReadGroups I="$in" O="$ALIGNMENT_DIR/${n}.rg.bam" \
       RGID=1 RGLB=lib RGPL=ILM RGPU=unit RGSM="$n"
  java -jar "$PICARD_JAR" MarkDuplicates I="$ALIGNMENT_DIR/${n}.rg.bam" \
       O="$ALIGNMENT_DIR/${n}.dedup.bam" M="$LOG_DIR/${n}.metrics.txt" REMOVE_DUPLICATES=true
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
  samtools view -h "$ALIGNMENT_DIR/${n}.dedup.bam" | awk "$CMD" | \
    samtools view -bS - > "$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
  samtools index "$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
done

###############################################################################
# 11  MACS2 PEAKS (replicate / merged / pooled)                               #
###############################################################################
mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}

# A  replicate peaks
if ((${#CTRL_NAMES[@]})); then
  CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"
  merge_bams "$CTRL_MRG" "${CTRL_NAMES[@]/%/.dedup.filtered.bam}"
fi
for n in "${TREAT_NAMES[@]}"; do
  if ((${#CTRL_NAMES[@]})); then
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" -c "$CTRL_MRG" \
           -f BAMPE -g $GENOME_SIZE -n "$n" --outdir "$PEAK_DIR/replicate"
  else
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" \
           -f BAMPE -g $GENOME_SIZE -n "$n" --outdir "$PEAK_DIR/replicate"
  fi
done

# B merged peaks
T_MRG="$ALIGNMENT_DIR/treat_merged.bam"
merge_bams "$T_MRG" "${TREAT_NAMES[@]/%/.dedup.filtered.bam}"
if ((${#CTRL_NAMES[@]})); then
  $MACS2 callpeak -t "$T_MRG" -c "$CTRL_MRG" -f BAMPE -g $GENOME_SIZE \
         -n treatMerged_vs_ctrlMerged --outdir "$PEAK_DIR/merged"
else
  $MACS2 callpeak -t "$T_MRG" -f BAMPE -g $GENOME_SIZE \
         -n treatMerged --outdir "$PEAK_DIR/merged"
fi

# C pooled-control vs each replicate
if ((${#CTRL_NAMES[@]})); then
  for n in "${TREAT_NAMES[@]}"; do
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" -c "$CTRL_MRG" \
           -f BAMPE -g $GENOME_SIZE -n "${n}_vs_ctrlPooled" --outdir "$PEAK_DIR/pooled"
  done
fi

###############################################################################
# 12  SPIKE SCALE FACTORS                                                     #
###############################################################################
declare -A SCALE
read_count(){ samtools view -c -F 2304 "$1"; }
for n in "${SAMPLES[@]}"; do
  h=$(read_count "$ALIGNMENT_DIR/${n}.dedup.filtered.bam")
  s=$(read_count "$SPIKE_DIR/${n}.ecoli.sorted.bam")
  ((s)) && SCALE["$n"]=$(awk -v h=$h -v s=$s 'BEGIN{printf "%.6f",h/s}')
done

###############################################################################
# 13  BIGWIG                                                                  #
###############################################################################
for n in "${SAMPLES[@]}"; do
  bg="$BW_DIR/${n}.bedgraph"
  bam_to_bedgraph "$ALIGNMENT_DIR/${n}.dedup.filtered.bam" "$bg"
  [[ -n ${SCALE[$n]} ]] && awk -v f="${SCALE[$n]}" '{$4*=$4*f;print}' "$bg" > "${bg}.tmp" && mv "${bg}.tmp" "$bg"
  bedGraphToBigWig "$bg" "$CHROM_SIZE" "$BW_DIR/${n}.bw"
done

###############################################################################
# 14  ANNOTATE PEAKS                                                          #
###############################################################################
for np in "$PEAK_DIR"/{replicate,merged,pooled}/*.narrowPeak; do
  [[ -f $np ]] || continue
  base=${np##*/}; base=${base%.narrowPeak}
  full="$ANN_DIR/${base}.annotated.bed"
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb > "$full"
done

log "DONE" "Outputs in $OUTPUT_DIR"
