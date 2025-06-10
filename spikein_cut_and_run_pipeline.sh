#!/bin/bash

# -----------------------------------------------------------------------------
# CUT&RUN PIPELINE WITH E. coli SPIKE‑IN SUPPORT
# -----------------------------------------------------------------------------
# This script processes CUT&RUN paired‑end FASTQs through trimming, alignment,
# duplicate removal, fragment filtering, peak calling, BigWig generation, and
# optional gene‑feature annotation.  New in this version → automatic alignment
# of trimmed reads to an E. coli spike‑in genome, producing per‑sample BAMs that
# feed into MultiQC for QC/normalisation.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 0  Load parameters from config.json
# -----------------------------------------------------------------------------
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

BROAD_EXTSIZE=$( jq -r '.broad_peak_extsize'  "$CONFIG_FILE")
NARROW_EXTSIZE=$(jq -r '.narrow_peak_extsize' "$CONFIG_FILE")

# -----------------------------------------------------------------------------
# 1  Create required directories
# -----------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports";    mkdir -p "$FASTQC_DIR"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks";         mkdir -p "$PEAK_DIR"
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs";      mkdir -p "$BW_DIR"
ANN_DIR="$OUTPUT_DIR/annotated_peaks";      mkdir -p "$ANN_DIR"
MULTIQC_DIR="$OUTPUT_DIR/multiqc_reports";  mkdir -p "$MULTIQC_DIR"

# -----------------------------------------------------------------------------
# 2  Utility functions
# -----------------------------------------------------------------------------
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
               --genomeDir  "$index" \
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

  STAR --runThreadN "$NUM_THREADS" \
       --genomeDir "$ECOLI_INDEX" \
       --readFilesIn "$r1" "$r2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix "$SPIKE_BAM_DIR/${sample}_ecoli_" \
       > "$LOG_DIR/STAR_${sample}_ecoli.log" \
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

# -----------------------------------------------------------------------------
# 3  Derive filenames for downstream steps
# -----------------------------------------------------------------------------
TREATMENT_BASE=$(get_sample_basename "$TREATMENT_R1")
TREATMENT_TRIMMED_R1="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R1.fq.gz"
TREATMENT_TRIMMED_R2="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R2.fq.gz"

if [[ -n "$CONTROL_R1" ]]; then
  CONTROL_BASE=$(get_sample_basename "$CONTROL_R1")
  CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
  CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"
fi

# -----------------------------------------------------------------------------
# 4  Compute numeric genome size
# -----------------------------------------------------------------------------
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

echo "Using genome size $GENOME_SIZE for key $GENOME_SIZE_STRING" | tee -a "$LOG_DIR/pipeline.log"

# -----------------------------------------------------------------------------
# 5  FASTQC (raw reads)
# -----------------------------------------------------------------------------
FASTQ_FILES=("$TREATMENT_R1" "$TREATMENT_R2")
if [[ -n "$CONTROL_R1" ]]; then FASTQ_FILES+=("$CONTROL_R1" "$CONTROL_R2"); fi

echo "Running FastQC…" | tee -a "$LOG_DIR/pipeline.log"
for fq in "${FASTQ_FILES[@]}"; do
  fastqc --extract -o "$FASTQC_DIR" "$fq" >> "$LOG_DIR/pipeline.log" 2>&1
done

#~ # -----------------------------------------------------------------------------
#~ # 6  Adapter trimming (Trim Galore!)  –  trim ALL declared FASTQ pairs
#~ # -----------------------------------------------------------------------------
#~ echo "[Trim Galore] starting…" | tee -a "$LOG_DIR/pipeline.log"

#~ i=0
#~ while [[ $i -lt ${#FASTQ_FILES[@]} ]]; do
  #~ R1=${FASTQ_FILES[$i]}
  #~ R2=${FASTQ_FILES[$((i+1))]}
  #~ BASE=$(get_sample_basename "$R1")

  #~ echo "  ↳ trimming $BASE" | tee -a "$LOG_DIR/pipeline.log"
  #~ trim_galore --paired --quality 20 --phred33 \
              #~ --output_dir "$ALIGNMENT_DIR" "$R1" "$R2" \
              #~ > "$LOG_DIR/trim_${BASE}.log" 2>&1

  #~ VAL1=$(find "$ALIGNMENT_DIR" -name "*_val_1.fq.gz" | grep "$BASE" | head -n1)
  #~ VAL2=$(find "$ALIGNMENT_DIR" -name "*_val_2.fq.gz" | grep "$BASE" | head -n1)

  #~ if [[ -f "$VAL1" && -f "$VAL2" ]]; then
    #~ mv "$VAL1" "$ALIGNMENT_DIR/${BASE}_trimmed_R1.fq.gz"
    #~ mv "$VAL2" "$ALIGNMENT_DIR/${BASE}_trimmed_R2.fq.gz"
  #~ else
    #~ echo "❌ Trim Galore did not produce trimmed files for $BASE — skipping." | tee -a "$LOG_DIR/pipeline.log"
  #~ fi
  #~ i=$((i+2))
#~ done

#-----------------------------------------------------------------------
# Decide once if a usable control exists (for downstream steps)
#-----------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 7  Spike-in alignment (E. coli)
# -----------------------------------------------------------------------------
echo "Aligning to E. coli genome with STAR for subsequent spike-in scaling…" | tee -a "$LOG_DIR/pipeline.log"
run_spikein_align "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" "$TREATMENT_BASE"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_spikein_align "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" "$CONTROL_BASE"
fi

# -----------------------------------------------------------------------------
# 8  Host‑genome alignment (STAR)
# -----------------------------------------------------------------------------
echo "Aligning to host genome with STAR…" | tee -a "$LOG_DIR/pipeline.log"
run_star "$REFERENCE_GENOME" \
         "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" \
         "$ALIGNMENT_DIR/${TREATMENT_BASE}." "STAR_${TREATMENT_BASE}"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_star "$REFERENCE_GENOME" \
           "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" \
           "$ALIGNMENT_DIR/${CONTROL_BASE}." "STAR_${CONTROL_BASE}"
fi

# -----------------------------------------------------------------------------
# 9  Picard: Add RG + MarkDuplicates
# -----------------------------------------------------------------------------
PICARD_PATH="/mnt/data/home/aviv/tools/picard.jar"
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

# -----------------------------------------------------------------------------
# 10  Fragment‑length filtering
# -----------------------------------------------------------------------------
case "$FRAGMENT_SIZE_FILTER" in
  histones)              FRAG_CMD='{if ($9 >= 130 && $9 <= 300 || $1 ~ /^@/) print $0}' ;;
  transcription_factors) FRAG_CMD='{if ($9 < 130 || $1 ~ /^@/) print $0}' ;;
  *)                     FRAG_CMD='{if ($9 < 1000 || $1 ~ /^@/) print $0}' ;;
esac

echo "Filtering fragments ($FRAGMENT_SIZE_FILTER)…" | tee -a "$LOG_DIR/pipeline.log"
for bam in "$ALIGNMENT_DIR"/*.dedup.bam; do
  base=$(basename "$bam" .dedup.bam)
  samtools view -h "$bam" | awk "$FRAG_CMD" | samtools view -bS - > "$ALIGNMENT_DIR/${base}.dedup.filtered.bam"
done

# -----------------------------------------------------------------------------
# 11  Peak calling (MACS2)
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# 12  BigWig generation
# -----------------------------------------------------------------------------
for bam in "$ALIGNMENT_DIR"/*.dedup.filtered.bam; do
  base=$(basename "$bam" .dedup.filtered.bam)
  bedtools genomecov -ibam "$bam" -bg > "$BW_DIR/${base}.bedgraph"
  bedGraphToBigWig "$BW_DIR/${base}.bedgraph" "$CHROM_SIZE" "$BW_DIR/${base}.bw"
done

# -----------------------------------------------------------------------------
# 13  Peak annotation (optional)
# -----------------------------------------------------------------------------
for np in "$PEAK_DIR"/*.narrowPeak; do
  base=$(basename "$np" .narrowPeak)
  bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb > "$ANN_DIR/${base}.annotated.bed"

  # concise TSV: peak coords + gene name + strand
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$11,$12}' \
    "$ANN_DIR/${base}.annotated.bed" > "$ANN_DIR/${base}.annotated.tsv"
done

# -----------------------------------------------------------------------------
# 14 MultiQC summary                                                          
# -----------------------------------------------------------------------------
echo "[MultiQC] aggregating reports" | tee -a "$LOG_DIR/pipeline.log"
rm "$FASTQC_DIR"/*.zip 2>/dev/null # Clean up FastQC zip files to avoid MultiQC parsing issues
multiqc "$FASTQC_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR" \
       -o "$MULTIQC_DIR" \
       > "$LOG_DIR/multiqc.log" 2> "$LOG_DIR/multiqc_error.log"

echo "🎉  CUT&RUN pipeline complete!  Results are in: $OUTPUT_DIR" | tee -a "$LOG_DIR/pipeline.log"

#~ # -----------------------------------------------------------------------------
#~ # 11  MultiQC summary
#~ # -----------------------------------------------------------------------------
#~ echo "[MultiQC] aggregating QC reports" | tee -a "$LOG_DIR/pipeline.log"

#~ if [[ $USE_CONTROL -eq 1 ]]; then
  #~ multiqc "$FASTQC_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR" \
          #~ -o "$MULTIQC_DIR" \
          #~ > "$LOG_DIR/multiqc.log" 2> "$LOG_DIR/multiqc_error.log"
#~ else
  #~ multiqc "$FASTQC_DIR" "$ALIGNMENT_DIR" \
          #~ -o "$MULTIQC_DIR" \
          #~ > "$LOG_DIR/multiqc.log" 2> "$LOG_DIR/multiqc_error.log"
#~ fi
