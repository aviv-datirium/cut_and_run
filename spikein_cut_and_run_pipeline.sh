#!/bin/bash

cat <<'BANNER'

# ------------------------------------------------------------------------------
# CUT&RUN PIPELINE WITH E. coli SPIKE‑IN SCALING
# ------------------------------------------------------------------------------
# This script processes CUT&RUN paired‑end FASTQs through trimming, alignment,
# duplicate removal, fragment filtering, peak calling, BigWig generation,
# gene‑feature annotation, and graphical reporting. This version includes
# alignment of trimmed reads to an E. coli genome for normalization and scaling.
# ------------------------------------------------------------------------------

################################################################################
#              STEPS                                                           #
################################################################################
# ---- Preparatory steps                                                       # 
# 0  Load parameters from config.json                                          #
# 1  Paths to tools and software                                               #
# 2  Create required directories                                               #
# 3  Utility functions                                                         #
# 4  Derive filenames for downstream steps                                     #
# 5  Compute numeric genome size                                               #
# ---- Computational steps                                                     #
# 6  FASTQC (raw reads)                                                        #
# 7  Adapter trimming (Trim Galore!)  –  trim ALL declared FASTQ pairs         #
# 8  Include/exclude control sample logic (for downstream steps)               #
# 9  Spike-in alignment (E. coli) with STAR                                    #
# 10  Host‑genome alignment (STAR)                                             #
# 11  Picard: Add RG + MarkDuplicates                                          #
# 12  Fragment‑length filtering                                                #   
# 13  Peak calling (MACS2)                                                     #
# 14  Spike-in scaling factors (optional)                                      #
# 15  BigWig generation (with optional scaling)                                #
# 16  Peak annotation (optional)                                               #
# 17  MultiQC summary                                                          #
################################################################################

BANNER

set -euo pipefail

# ------------------------------------------------------------------------------
# 0  Load parameters from config.json
# ------------------------------------------------------------------------------
CONFIG_FILE="/mnt/data/home/aviv/cut_and_run/config.json"

RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir'              "$CONFIG_FILE")
ALIGNMENT_DIR=$(jq -r '.alignment_dir'              "$CONFIG_FILE")
OUTPUT_DIR=$(jq  -r '.output_dir'                   "$CONFIG_FILE")
LOG_DIR=$(jq     -r '.log_dir'                      "$CONFIG_FILE")

TREATMENT_R1=$(jq -r '.samples.treatment.r1'        "$CONFIG_FILE")
TREATMENT_R2=$(jq -r '.samples.treatment.r2'        "$CONFIG_FILE")
CONTROL_R1=$(  jq -r '.samples.control.r1 // empty' "$CONFIG_FILE")
CONTROL_R2=$(  jq -r '.samples.control.r2 // empty' "$CONFIG_FILE")

REFERENCE_GENOME=$(jq -r '.reference_genome'        "$CONFIG_FILE")
ANNOTATION_GENES=$(jq -r '.annotation_genes'        "$CONFIG_FILE")
CHROM_SIZE=$(      jq -r '.chrom_sizes'             "$CONFIG_FILE")
ECOLI_INDEX=$(     jq -r '.ecoli_index'             "$CONFIG_FILE")

# Optional explicit spike‑in BAM directory; default → <alignment>/spikein
SPIKE_BAM_DIR=$(jq -r '.spike_bam_dir // empty'     "$CONFIG_FILE")
[ -z "$SPIKE_BAM_DIR" ] || [ "$SPIKE_BAM_DIR" = "null" ] && \
  SPIKE_BAM_DIR="$ALIGNMENT_DIR/spikein"

GENOME_SIZE_STRING=$(jq -r '.genome_size'           "$CONFIG_FILE")
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$( jq -r '.custom_genome_size'    "$CONFIG_FILE")
NUM_THREADS=$(       jq -r '.num_threads'            "$CONFIG_FILE")

# In common in cases like CUT&RUN or when read count is low or fragment length is narrow (as expected in histone or TF targeting experiments), MACS2 simply falls back to non-model-based peak calling and recommends
BROAD_EXTSIZE=$( jq -r '.broad_peak_extsize'  "$CONFIG_FILE")
NARROW_EXTSIZE=$(jq -r '.narrow_peak_extsize' "$CONFIG_FILE")

#------------------------------------------------------------------------------
# 1 Paths to tools and software
#------------------------------------------------------------------------------
# Picard tools path
PICARD_PATH="/mnt/data/home/aviv/tools/picard.jar"  # Path to the Picard jar file (e.g., picard.jar)
# FastQC tools path
FASTQC_PATH="/mnt/data/home/aviv/tools/FastQC/fastqc"
# STAR path
STAR_PATH="/mnt/data/home/aviv/tools/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR"

# ------------------------------------------------------------------------------
# 2  Create required directories
# ------------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports";    mkdir -p "$FASTQC_DIR"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks";         mkdir -p "$PEAK_DIR"
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs";      mkdir -p "$BW_DIR"
ANN_DIR="$OUTPUT_DIR/annotated_peaks";      mkdir -p "$ANN_DIR"
MULTIQC_DIR="$OUTPUT_DIR/multiqc_reports";  mkdir -p "$MULTIQC_DIR"

# ------------------------------------------------------------------------------
# 3  Utility functions
# ------------------------------------------------------------------------------
get_sample_basename() {
  local r1=$1; local base=$(basename "$r1")
  base=${base%.fastq.gz}; base=${base%.fq.gz}; base=${base%_R1}; \
  base=${base%_1}; base=${base%.R1}; base=${base%.1}
  echo "$base" | sed 's/[^a-zA-Z0-9._-]//g'
}

# Usage: run_star <STAR_index> <R1> <R2> <outPrefix> <logBase>
# Runs STAR → sorted BAM + logs.
run_star() {
  local index=$1 r1=$2 r2=$3 prefix=$4 logbase=$5
  "$STAR_PATH" --runThreadN "$NUM_THREADS" \
               --genomeDir  "$REFERENCE_GENOME" \
               --readFilesIn "$r1" "$r2" \
               --readFilesCommand zcat \
               --outSAMtype BAM SortedByCoordinate \
               --outFileNamePrefix "$prefix" \
               > "$LOG_DIR/${logbase}.log"  \
               2> "$LOG_DIR/${logbase}_err.log"
}

run_spikein_align() {
  # Align trimmed reads to E. coli spike‑in genome with STAR.
  # Usage: run_spikein_align R1 R2 SAMPLE_BASENAME
  local r1=$1 r2=$2 sample=$3
  echo "[SPIKE‑IN] Aligning $sample to E. coli genome…" | tee -a "$LOG_DIR/pipeline.log"
  "$STAR_PATH" --runThreadN "$NUM_THREADS" \
               --genomeDir  "$ECOLI_INDEX" \
               --readFilesIn "$r1" "$r2" \
               --readFilesCommand zcat \
               --outSAMtype BAM SortedByCoordinate \
               --outFileNamePrefix "$SPIKE_BAM_DIR/${sample}_ecoli_" \
               > "$LOG_DIR/STAR_${sample}_ecoli.log"  \
               2> "$LOG_DIR/STAR_${sample}_ecoli_error.log"

  local tmp="$SPIKE_BAM_DIR/${sample}_ecoli_Aligned.sortedByCoord.out.bam"
  if [[ -f "$tmp" ]]; then
    mv "$tmp" "$SPIKE_BAM_DIR/${sample}.ecoli.sorted.bam"
    samtools index "$SPIKE_BAM_DIR/${sample}.ecoli.sorted.bam"
  else
    echo "❌ Spike‑in alignment failed for $sample — BAM not produced." | tee -a "$LOG_DIR/pipeline.log"
    exit 1
  fi
}

# ------------------------------------------------------------------------------
# 4  Derive filenames for downstream steps
# ------------------------------------------------------------------------------
TREATMENT_BASE=$(get_sample_basename "$TREATMENT_R1")
TREATMENT_TRIMMED_R1="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R1.fq.gz"
TREATMENT_TRIMMED_R2="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R2.fq.gz"

if [[ -n "$CONTROL_R1" ]]; then
  CONTROL_BASE=$(get_sample_basename "$CONTROL_R1")
  CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
  CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"
fi

# ------------------------------------------------------------------------------
# 5  Compute numeric genome size
# ------------------------------------------------------------------------------
GENOME_SIZE_HUMAN=2913022398; GENOME_SIZE_MOUSE=2652783500
GENOME_SIZE_DROSOPHILA=165000000; GENOME_SIZE_CELEGANS=1000000000
GENOME_SIZE_YEAST=12000000

case "$GENOME_SIZE_STRING" in
  hs) GENOME_SIZE=$GENOME_SIZE_HUMAN ;;
  mm) GENOME_SIZE=$GENOME_SIZE_MOUSE ;;
  dm) GENOME_SIZE=$GENOME_SIZE_DROSOPHILA ;;
  ce) GENOME_SIZE=$GENOME_SIZE_CELEGANS ;;
  sc) GENOME_SIZE=$GENOME_SIZE_YEAST ;;
  *)  if [[ "$CUSTOM_GENOME_SIZE" != "null" ]]; then
        GENOME_SIZE=$CUSTOM_GENOME_SIZE
      else
        echo "Error: Invalid genome size string $GENOME_SIZE_STRING" | tee -a "$LOG_DIR/pipeline.log"
        exit 1
      fi ;;
esac

echo "Using genome size $GENOME_SIZE for host genome: $GENOME_SIZE_STRING" | tee -a "$LOG_DIR/pipeline.log"

# ------------------------------------------------------------------------------
# 6  FASTQC (raw reads)
# ------------------------------------------------------------------------------
FASTQ_FILES=("$TREATMENT_R1" "$TREATMENT_R2")
if [[ -n "$CONTROL_R1" ]]; then FASTQ_FILES+=("$CONTROL_R1" "$CONTROL_R2"); fi

echo "Running FastQC…" | tee -a "$LOG_DIR/pipeline.log"
for fq in "${FASTQ_FILES[@]}"; do
  $FASTQC_PATH --extract -o "$FASTQC_DIR" "$fq" >> "$LOG_DIR/pipeline.log" 2>&1
done

# ------------------------------------------------------------------------------
# 7  Adapter trimming (Trim Galore!)  –  trim ALL declared FASTQ pairs
# ------------------------------------------------------------------------------
echo "[Trim Galore] starting…" | tee -a "$LOG_DIR/pipeline.log"
i=0
while [[ $i -lt ${#FASTQ_FILES[@]} ]]; do
  R1=${FASTQ_FILES[$i]}
  R2=${FASTQ_FILES[$((i+1))]}
  BASE=$(get_sample_basename "$R1")

  echo "  ↳ trimming $BASE" | tee -a "$LOG_DIR/pipeline.log"
  trim_galore --paired --quality 20 --phred33 \
              --output_dir "$ALIGNMENT_DIR" "$R1" "$R2" \
              > "$LOG_DIR/trim_${BASE}.log" 2>&1

  VAL1=$(find "$ALIGNMENT_DIR" -name "*_val_1.fq.gz" | grep "$BASE" | head -n1)
  VAL2=$(find "$ALIGNMENT_DIR" -name "*_val_2.fq.gz" | grep "$BASE" | head -n1)

  if [[ -f "$VAL1" && -f "$VAL2" ]]; then
    mv "$VAL1" "$ALIGNMENT_DIR/${BASE}_trimmed_R1.fq.gz"
    mv "$VAL2" "$ALIGNMENT_DIR/${BASE}_trimmed_R2.fq.gz"
  else
    echo "❌ Trim Galore did not produce trimmed files for $BASE — skipping." | tee -a "$LOG_DIR/pipeline.log"
  fi
  i=$((i+2))
done

# ------------------------------------------------------------------------------
# 8 Include/exclude control sample logic (for downstream steps)
# ------------------------------------------------------------------------------
USE_CONTROL=0
if [[ -n "${CONTROL_BASE:-}" ]] && \
   [[ -f "$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz" ]] && \
   [[ -f "$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz" ]]; then
  USE_CONTROL=1
  CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
  CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"
else
  echo "⚠️  No trimmed control FASTQs found — proceeding without control." | tee -a "$LOG_DIR/pipeline.log"
fi

# ------------------------------------------------------------------------------
# 9  Spike-in alignment (E. coli)
# ------------------------------------------------------------------------------
echo "Aligning to E. coli genome with STAR for subsequent spike-in scaling…" | tee -a "$LOG_DIR/pipeline.log"
run_spikein_align "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" "$TREATMENT_BASE"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_spikein_align "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" "$CONTROL_BASE"
fi

# ------------------------------------------------------------------------------
# 10  Host‑genome alignment (STAR)
# ------------------------------------------------------------------------------
echo "Aligning to host genome with STAR…" | tee -a "$LOG_DIR/pipeline.log"
run_star "$REFERENCE_GENOME" \
         "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" \
         "$ALIGNMENT_DIR/${TREATMENT_BASE}." "STAR_${TREATMENT_BASE}"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_star "$REFERENCE_GENOME" \
           "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" \
           "$ALIGNMENT_DIR/${CONTROL_BASE}." "STAR_${CONTROL_BASE}"
fi

# ------------------------------------------------------------------------------
# 11  Picard: Add RG + MarkDuplicates
# ------------------------------------------------------------------------------
echo "Running Picard (RG + dedup)…" | tee -a "$LOG_DIR/pipeline.log"
for bam in "$ALIGNMENT_DIR"/*.Aligned.sortedByCoord.out.bam; do
  [[ ! -s "$bam" ]] && { echo "⚠️ Empty BAM $bam — skipping" | tee -a "$LOG_DIR/pipeline.log"; continue; }

  base=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
  java -jar "$PICARD_PATH" AddOrReplaceReadGroups I="$bam" \
       O="$ALIGNMENT_DIR/${base}.rg.bam" \
       RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM="$base" \
       VALIDATION_STRINGENCY=LENIENT || { echo "AddRG failed for $base" | tee -a "$LOG_DIR/pipeline.log"; exit 1; }

  java -jar "$PICARD_PATH" MarkDuplicates I="$ALIGNMENT_DIR/${base}.rg.bam" \
       O="$ALIGNMENT_DIR/${base}.dedup.bam" \
       M="$LOG_DIR/${base}.metrics.txt" REMOVE_DUPLICATES=true \
       VALIDATION_STRINGENCY=LENIENT || { echo "MarkDuplicates failed for $base" | tee -a "$LOG_DIR/pipeline.log"; exit 1; }
done

# ------------------------------------------------------------------------------
# 12  Fragment‑length filtering
# ------------------------------------------------------------------------------
case "$FRAGMENT_SIZE_FILTER" in
  histones)              FRAG_CMD='{if ($9 >= 130 && $9 <= 300 || $1 ~ /^@/) print $0}' ;;
  transcription_factors) FRAG_CMD='{if ($9 < 130 || $1 ~ /^@/) print $0}' ;;
  *)                     FRAG_CMD='{if ($9 < 1000 || $1 ~ /^@/) print $0}' ;;
esac

echo "Filtering fragments by ($FRAGMENT_SIZE_FILTER)…" | tee -a "$LOG_DIR/pipeline.log"
for bam in "$ALIGNMENT_DIR"/*.dedup.bam; do
  base=$(basename "$bam" .dedup.bam)
  samtools view -h "$bam" | awk "$FRAG_CMD" | samtools view -bS - > "$ALIGNMENT_DIR/${base}.dedup.filtered.bam"
done

# ------------------------------------------------------------------------------
# 13  Peak calling (MACS2)
# ------------------------------------------------------------------------------
echo "Running MACS2 for peak calling (both broad and narrow peaks)..." | tee -a "$LOG_DIR/pipeline.log"
TREATMENT_FILTERED="$ALIGNMENT_DIR/${TREATMENT_BASE}.dedup.filtered.bam"
[[ -f "$TREATMENT_FILTERED" ]] || { echo "❌ Treatment BAM missing"; exit 1; }

if [[ $USE_CONTROL -eq 1 ]]; then
  CONTROL_FILTERED="$ALIGNMENT_DIR/${CONTROL_BASE}.dedup.filtered.bam"
  macs2 callpeak -t "$TREATMENT_FILTERED" -c "$CONTROL_FILTERED" \
        --broad         --outdir "$PEAK_DIR" -n "$TREATMENT_BASE"
  macs2 callpeak -t "$TREATMENT_FILTERED" -c "$CONTROL_FILTERED" \
        --call-summits  --outdir "$PEAK_DIR" -n "$TREATMENT_BASE"
else
  macs2 callpeak -t "$TREATMENT_FILTERED" \
        --broad --nomodel --extsize "$BROAD_EXTSIZE" \
        --outdir "$PEAK_DIR" -n "$TREATMENT_BASE"
  macs2 callpeak -t "$TREATMENT_FILTERED" \
        --call-summits --nomodel --extsize "$NARROW_EXTSIZE" \
        --outdir "$PEAK_DIR" -n "$TREATMENT_BASE"
fi

# ------------------------------------------------------------------------------
# 14  Spike-in scaling factors (optional)                                    
# ------------------------------------------------------------------------------
echo "[Spike-in] calculating scale factors" | tee -a "$LOG_DIR/pipeline.log"
declare -A SCALE  # associative array sample → factor

# function to count uniquely mapped reads
count_reads () { samtools view -c -F 2304 "$1"; }

for host_bam in "$ALIGNMENT_DIR"/*.dedup.filtered.bam; do
  samp=$(basename "$host_bam" .dedup.filtered.bam)

  ecoli_bam="$SPIKE_BAM_DIR/${samp}.ecoli.sorted.bam"
  if [[ -f "$ecoli_bam" ]]; then
    host_reads=$(count_reads "$host_bam")
    spike_reads=$(count_reads "$ecoli_bam")

    if (( spike_reads > 0 )); then
      factor=$(awk -v h=$host_reads -v s=$spike_reads 'BEGIN{printf "%.6f", h/s}')
      echo "  ↳ $samp : host=$host_reads  spike=$spike_reads  scale=$factor" \
           | tee -a "$LOG_DIR/pipeline.log"
      SCALE["$samp"]=$factor
    else
      echo "  ↳ $samp : spike-in reads = 0 — skipping scaling" \
           | tee -a "$LOG_DIR/pipeline.log"
    fi
  fi
done

# ------------------------------------------------------------------------------
# 15  BigWig generation (with optional scaling)                                
# ------------------------------------------------------------------------------
echo "Generating bigwig for track viewing in IGV..." | tee -a "$LOG_DIR/pipeline.log"
for host_bam in "$ALIGNMENT_DIR"/*.dedup.filtered.bam; do
  samp=$(basename "$host_bam" .dedup.filtered.bam)

  bedgraph="$BW_DIR/${samp}.bedgraph"
  bigwig="$BW_DIR/${samp}.bw"

  if [[ -n "${SCALE[$samp]:-}" ]]; then
    # scale coverage linearly by spike-in factor
    bedtools genomecov -ibam "$host_bam" -bg | \
      awk -v f="${SCALE[$samp]}" '{ $4=$4*f; print }' > "$bedgraph"
  else
    bedtools genomecov -ibam "$host_bam" -bg > "$bedgraph"
  fi

  bedGraphToBigWig "$bedgraph" "$CHROM_SIZE" "$bigwig"
done

# ------------------------------------------------------------------------------
# 16  Peak annotation
# ------------------------------------------------------------------------------
echo "Annotating peaks with bedtools intersect..." | tee -a "$LOG_DIR/pipeline.log"
for np in "$PEAK_DIR"/*.narrowPeak; do
  base=$(basename "$np" .narrowPeak)
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb > "$ANN_DIR/${base}.annotated.bed"

  # concise TSV: peak coords + gene name + strand
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$11,$12}' \
    "$ANN_DIR/${base}.annotated.bed" > "$ANN_DIR/${base}.annotated.tsv"
done

# ------------------------------------------------------------------------------
# 17 MultiQC summary                                                          
# ------------------------------------------------------------------------------
echo "[MultiQC] aggregating reports" | tee -a "$LOG_DIR/pipeline.log"
rm "$FASTQC_DIR"/*.zip 2>/dev/null # Clean up FastQC zip files to avoid MultiQC parsing issues
multiqc "$FASTQC_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR" \
       -o "$MULTIQC_DIR" \
       > "$LOG_DIR/multiqc.log" 2> "$LOG_DIR/multiqc_error.log"

echo "🎉  CUT&RUN pipeline complete!  Results are in: $OUTPUT_DIR" | tee -a "$LOG_DIR/pipeline.log"
