#!/bin/bash

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

set -o pipefail

CONFIG_FILE="/mnt/data/home/aviv/cut_and_run/config_replicates.json"

###############################################################################
# 0  LOAD PATHS FROM CONFIG                                                   #
###############################################################################
RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir'  "$CONFIG_FILE")
ALIGNMENT_DIR=$(jq -r '.alignment_dir'  "$CONFIG_FILE")
OUTPUT_DIR=$( jq -r '.output_dir'       "$CONFIG_FILE")
LOG_DIR=$(    jq -r '.log_dir'          "$CONFIG_FILE")

REFERENCE_GENOME=$(jq -r '.reference_genome'   "$CONFIG_FILE")
ECOLI_INDEX=$(     jq -r '.ecoli_index'        "$CONFIG_FILE")
ANNOTATION_GENES=$(jq -r '.annotation_genes'   "$CONFIG_FILE")
CHROM_SIZE=$(      jq -r '.chrom_sizes'        "$CONFIG_FILE")

GENOME_SIZE_STRING=$(jq -r '.genome_size'      "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' "$CONFIG_FILE")
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
NUM_THREADS=$(jq    -r '.num_threads'          "$CONFIG_FILE")

# FASTQ replicate lists
TREAT_R1=($(jq -r '.samples.treatment[]?.r1' "$CONFIG_FILE"))
TREAT_R2=($(jq -r '.samples.treatment[]?.r2' "$CONFIG_FILE"))
CTRL_R1=($( jq -r '.samples.control[]?.r1 // empty' "$CONFIG_FILE"))
CTRL_R2=($( jq -r '.samples.control[]?.r2 // empty' "$CONFIG_FILE"))

NUM_TREAT=${#TREAT_R1[@]}
NUM_CTRL=${#CTRL_R1[@]}

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
# 2  CREATE DIRECTORIES                                                       #
###############################################################################
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports"
SPIKE_DIR="$ALIGNMENT_DIR/spikein"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks"
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs"
ANN_DIR="$OUTPUT_DIR/annotated_peaks"
for d in "$FASTQC_DIR" "$SPIKE_DIR" "$PEAK_DIR"/{replicate,merged,pooled} \
         "$BW_DIR" "$ANN_DIR" ; do mkdir -p "$d"; done

###############################################################################
# 3  UTILITY FUNCTIONS                                                        #
###############################################################################
get_sample_basename(){
  local f=${1##*/}; f=${f%.fastq.gz}; f=${f%.fq.gz}
  f=${f%_R1_001}; f=${f%_R2_001}; f=${f%_R1}; f=${f%_R2}; f=${f%_1}; f=${f%_2}
  printf '%s\n' "$f"
}

trim_one_pair(){        # $1 R1  $2 R2  $3 BASE
  java -jar "$TRIMMOMATIC_JAR" PE -threads "$NUM_THREADS" -phred33 \
       "$1" "$2" \
       "$ALIGNMENT_DIR/${3}_trimmed_R1.fq.gz" /dev/null \
       "$ALIGNMENT_DIR/${3}_trimmed_R2.fq.gz" /dev/null \
       ILLUMINACLIP:"$ADAPTER_FA":2:30:10 LEADING:5 TRAILING:5 \
       SLIDINGWINDOW:4:15 MINLEN:25 \
       > "$LOG_DIR/trim_${3}.log" 2>&1
}

run_star(){             # $1 R1  $2 R2  $3 OUTPREFIX  $4 LOGBASE  $5 GENOME
  "$STAR_BIN" --runThreadN "$NUM_THREADS" \
              --genomeDir "$5" \
              --readFilesIn "$1" "$2" \
              --readFilesCommand zcat \
              --outSAMtype BAM SortedByCoordinate \
              --outFileNamePrefix "$3" \
              > "$LOG_DIR/${4}.log" 2> "$LOG_DIR/${4}_err.log"
}

bam_to_bedgraph(){      # $1 in.bam  $2 out.bg
  bedtools genomecov -ibam "$1" -bg -pc | sort -k1,1 -k2,2n > "$2"
}

###############################################################################
# 4  BUILD BASENAME ARRAYS                                                    #
###############################################################################
TREAT_NAMES=(); for r1 in "${TREAT_R1[@]}"; do TREAT_NAMES+=( "$(get_sample_basename "$r1")" ); done
CTRL_NAMES=();  for r1 in "${CTRL_R1[@]}";  do CTRL_NAMES+=(  "$(get_sample_basename "$r1")" );  done

SAMPLES=("${TREAT_NAMES[@]}"); ((NUM_CTRL)) && SAMPLES+=("${CTRL_NAMES[@]}")

###############################################################################
# 5  FASTQC                                                                   #
###############################################################################
for ((i=0;i<NUM_TREAT;i++)); do $FASTQC_BIN -o "$FASTQC_DIR" --quiet "${TREAT_R1[$i]}" "${TREAT_R2[$i]}"; done
for ((i=0;i<NUM_CTRL;i++));  do $FASTQC_BIN -o "$FASTQC_DIR" --quiet "${CTRL_R1[$i]}"  "${CTRL_R2[$i]}";  done

###############################################################################
# 6  TRIM                                                                    #
###############################################################################
for ((i=0;i<NUM_TREAT;i++)); do trim_one_pair "${TREAT_R1[$i]}" "${TREAT_R2[$i]}" "${TREAT_NAMES[$i]}"; done
for ((i=0;i<NUM_CTRL;i++));  do trim_one_pair "${CTRL_R1[$i]}"  "${CTRL_R2[$i]}"  "${CTRL_NAMES[$i]}";  done

###############################################################################
# 7  SPIKE-IN ALIGNMENT (E. coli)                                             #
###############################################################################
for name in "${SAMPLES[@]}"; do
  run_star "$ALIGNMENT_DIR/${name}_trimmed_R1.fq.gz" \
           "$ALIGNMENT_DIR/${name}_trimmed_R2.fq.gz" \
           "$SPIKE_DIR/${name}_ecoli_"  "STAR_${name}_ecoli"  "$ECOLI_INDEX"
  mv "$SPIKE_DIR/${name}_ecoli_Aligned.sortedByCoord.out.bam" \
     "$SPIKE_DIR/${name}.ecoli.sorted.bam"
  samtools index "$SPIKE_DIR/${name}.ecoli.sorted.bam"
done

###############################################################################
# 8  HOST-GENOME ALIGNMENT (STAR)                                             #
###############################################################################
for name in "${SAMPLES[@]}"; do
  run_star "$ALIGNMENT_DIR/${name}_trimmed_R1.fq.gz" \
           "$ALIGNMENT_DIR/${name}_trimmed_R2.fq.gz" \
           "$ALIGNMENT_DIR/${name}." "STAR_${name}" "$REFERENCE_GENOME"
done

###############################################################################
# 9  PICARD RG + DEDUP                                                        #
###############################################################################
for name in "${SAMPLES[@]}"; do
  in_bam="$ALIGNMENT_DIR/${name}.Aligned.sortedByCoord.out.bam"
  [[ -s "$in_bam" ]] || continue
  java -jar "$PICARD_JAR" AddOrReplaceReadGroups I="$in_bam" \
       O="$ALIGNMENT_DIR/${name}.rg.bam" RGID=1 RGLB=lib RGPL=ILM RGPU=unit RGSM="$name"
  java -jar "$PICARD_JAR" MarkDuplicates I="$ALIGNMENT_DIR/${name}.rg.bam" \
       O="$ALIGNMENT_DIR/${name}.dedup.bam" M="$LOG_DIR/${name}.metrics.txt" REMOVE_DUPLICATES=true
  samtools index "$ALIGNMENT_DIR/${name}.dedup.bam"
done

###############################################################################
# 10  FRAGMENT FILTERING                                                      #
###############################################################################
case $FRAGMENT_SIZE_FILTER in
  histones)              FRAG_CMD='{if ($9>=130&&$9<=300 || $1~/^@/)print}';;
  transcription_factors) FRAG_CMD='{if ($9<130 || $1~/^@/)print}';;
  *)                     FRAG_CMD='{if ($9<1000||$1~/^@/)print}';;
esac
for name in "${SAMPLES[@]}"; do
  samtools view -h "$ALIGNMENT_DIR/${name}.dedup.bam" | \
    awk "$FRAG_CMD" | samtools view -bS - > "$ALIGNMENT_DIR/${name}.dedup.filtered.bam"
  samtools index "$ALIGNMENT_DIR/${name}.dedup.filtered.bam"
done

###############################################################################
# 11  PEAK CALLING – replicate / merged / pooled (MACS2)                      #
###############################################################################
mkdir -p "$PEAK_DIR"/{replicate,merged,pooled}
merge_bams(){ [[ -f "$1" ]] || { samtools merge -f "$1" "${@:2}" && samtools index "$1"; }; }

# A  replicate peaks
if ((NUM_CTRL)); then
  CTRL_MRG="$ALIGNMENT_DIR/control_merged.bam"
  merge_bams "$CTRL_MRG" "${CTRL_NAMES[@]/%/.dedup.filtered.bam}"
fi
for name in "${TREAT_NAMES[@]}"; do
  T_BAM="$ALIGNMENT_DIR/${name}.dedup.filtered.bam"
  if ((NUM_CTRL)); then
    $MACS2 callpeak -t "$T_BAM" -c "$CTRL_MRG" -f BAMPE -g $GENOME_SIZE \
           -n "$name" --outdir "$PEAK_DIR/replicate"
  else
    $MACS2 callpeak -t "$T_BAM" -f BAMPE -g $GENOME_SIZE \
           -n "$name" --outdir "$PEAK_DIR/replicate"
  fi
done

# B  merged peaks
T_MRG="$ALIGNMENT_DIR/treatment_merged.bam"
merge_bams "$T_MRG" "${TREAT_NAMES[@]/%/.dedup.filtered.bam}"
if ((NUM_CTRL)); then
  $MACS2 callpeak -t "$T_MRG" -c "$CTRL_MRG" -f BAMPE -g $GENOME_SIZE \
         -n treatmentMerged_vs_controlMerged --outdir "$PEAK_DIR/merged"
else
  $MACS2 callpeak -t "$T_MRG" -f BAMPE -g $GENOME_SIZE \
         -n treatmentMerged --outdir "$PEAK_DIR/merged"
fi

# C  pooled-control vs each replicate
if ((NUM_CTRL)); then
  for name in "${TREAT_NAMES[@]}"; do
    $MACS2 callpeak -t "$ALIGNMENT_DIR/${name}.dedup.filtered.bam" -c "$CTRL_MRG" \
           -f BAMPE -g $GENOME_SIZE -n "${name}_vs_ctrlPooled" --outdir "$PEAK_DIR/pooled"
  done
fi

###############################################################################
# 12  SPIKE-IN SCALE FACTORS                                                  #
###############################################################################
declare -A SCALE
count_reads(){ samtools view -c -F 2304 "$1"; }
for name in "${SAMPLES[@]}"; do
  host="$ALIGNMENT_DIR/${name}.dedup.filtered.bam"
  spike="$SPIKE_DIR/${name}.ecoli.sorted.bam"
  [[ -f $host && -f $spike ]] || continue
  h=$(count_reads "$host"); s=$(count_reads "$spike")
  ((s)) && SCALE["$name"]=$(awk -v h=$h -v s=$s 'BEGIN{printf "%.6f",h/s}')
done

###############################################################################
# 13  BIGWIG GENERATION                                                      #
###############################################################################
for name in "${SAMPLES[@]}"; do
  bg="$BW_DIR/${name}.bedgraph"
  bedtools genomecov -ibam "$ALIGNMENT_DIR/${name}.dedup.filtered.bam" -bg -pc > "$bg"
  [[ -n ${SCALE[$name]} ]] && \
      awk -v f="${SCALE[$name]}" '{ $4=$4*f; print }' "$bg" > "${bg}.tmp" && mv "${bg}.tmp" "$bg"
  bedGraphToBigWig "$bg" "$CHROM_SIZE" "$BW_DIR/${name}.bw"
done

###############################################################################
# 14  PEAK ANNOTATION                                                        #
###############################################################################
for np in "$PEAK_DIR"/{replicate,merged,pooled}/*.narrowPeak; do
  [[ -f $np ]] || continue
  base=$(basename "$np" .narrowPeak)
  full="$ANN_DIR/${base}.annotated.bed"
  tsv="$ANN_DIR/${base}.annotated.tsv"
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb > "$full"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10,$11,$12}' "$full" > "$tsv"
done

echo "🎉 Pipeline finished!  Outputs in $OUTPUT_DIR"
