#!/bin/bash

# PATHS: Read the config.json file to get input parameters

# Path to config
CONFIG_FILE="/mnt/data/home/aviv/cut_and_run/config.json"

# Input and output paths, incl. log and errors
RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir' $CONFIG_FILE)
ALIGNMENT_DIR=$(jq -r '.alignment_dir' $CONFIG_FILE)
OUTPUT_DIR=$(jq -r '.output_dir' $CONFIG_FILE)
LOG_DIR=$(jq -r '.log_dir' $CONFIG_FILE)

# Paths to treatment and optional control files
TREATMENT_R1=$(jq -r '.samples.treatment.r1' $CONFIG_FILE)
TREATMENT_R2=$(jq -r '.samples.treatment.r2' $CONFIG_FILE)
CONTROL_R1=$(jq -r '.samples.control.r1 // empty' $CONFIG_FILE)
CONTROL_R2=$(jq -r '.samples.control.r2 // empty' $CONFIG_FILE)

#~ # Defining base names for all samples
#~ TREATMENT_BASENAME=$(basename $TREATMENT_R1 | sed 's/_R1.*//;s/\.fastq.*//;s/[^a-zA-Z0-9_-]//g')
#~ CONTROL_BASENAME=""
#~ if [ -n "$CONTROL_R1" ]; then
    #~ CONTROL_BASENAME=$(basename $CONTROL_R1 | sed 's/_R1.*//;s/\.fastq.*//;s/[^a-zA-Z0-9_-]//g')
#~ fi

# Required data files
REFERENCE_GENOME=$(jq -r '.reference_genome' $CONFIG_FILE)
ANNOTATION_GENES=$(jq -r '.annotation_genes' $CONFIG_FILE)
CHROM_SIZE=$(jq -r '.chrom_sizes' $CONFIG_FILE)  # chromosome sizes for bedGraphToBigwig

# User run parameters
GENOME_SIZE_STRING=$(jq -r '.genome_size' $CONFIG_FILE)
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' $CONFIG_FILE)
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' $CONFIG_FILE)
NUM_THREADS=$(jq -r '.num_threads' $CONFIG_FILE)  # Get the number of threads from the config file

# Define genome sizes for various species (numeric values in base pairs)
GENOME_SIZE_HUMAN=2913022398  # Human genome size (hg38)
GENOME_SIZE_MOUSE=2652783500  # Mouse genome size (mm10)
GENOME_SIZE_DROSOPHILA=165000000  # Drosophila melanogaster (fruit fly) genome size
GENOME_SIZE_CELEGANS=1000000000  # Caenorhabditis elegans (nematode) genome size
GENOME_SIZE_YEAST=12000000  # Saccharomyces cerevisiae (yeast) genome size

# Paths to tools and software
# Picard tools path
PICARD_PATH="/mnt/data/home/aviv/tools/picard.jar"  # Path to the Picard jar file (e.g., picard.jar)
# FastQC tools path
FASTQC_PATH="/mnt/data/home/aviv/tools/FastQC/fastqc"
# STAR path
STAR_PATH="/mnt/data/home/aviv/tools/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR"

# Make output directories
mkdir -p $OUTPUT_DIR $ALIGNMENT_DIR
# Create FastQC output directory if it doesn't exist
mkdir -p $OUTPUT_DIR/fastqc_reports
# MACS2 peaks
mkdir -p $OUTPUT_DIR/macs2_peaks
# Create BigWig and BedGraph output directory if it doesn't exist
mkdir -p $OUTPUT_DIR/bigwig_bedgraphs
# Create annotated_peaks directory if it doesn't exist
mkdir -p $OUTPUT_DIR/annotated_peaks
# Create log and error files
mkdir -p $LOG_DIR

# Default fragment size threshold for filtering (below 1000 bp)
if [ "$FRAGMENT_SIZE_FILTER" == "histones" ]; then
    FRAGMENT_SIZE_FILTER="histones"
elif [ "$FRAGMENT_SIZE_FILTER" == "transcription_factors" ]; then
    FRAGMENT_SIZE_FILTER="transcription_factors"
else
    FRAGMENT_SIZE_FILTER="below_1000"
fi

# Map the GENOME_SIZE_STRING (user input from config.json) to the correct numeric genome size
case $GENOME_SIZE_STRING in
    "hs") GENOME_SIZE=$GENOME_SIZE_HUMAN ;;
    "mm") GENOME_SIZE=$GENOME_SIZE_MOUSE ;;
    "dm") GENOME_SIZE=$GENOME_SIZE_DROSOPHILA ;;
    "ce") GENOME_SIZE=$GENOME_SIZE_CELEGANS ;;
    "sc") GENOME_SIZE=$GENOME_SIZE_YEAST ;;
    *)
        if [ "$CUSTOM_GENOME_SIZE" != "null" ]; then
            # Use the custom genome size provided by the user
            GENOME_SIZE=$CUSTOM_GENOME_SIZE
        else
            echo "Error: Invalid genome size string provided in config.json. Please use one of: hs, mm, dm, ce, sc, or provide a custom size." | tee -a $LOG_DIR/pipeline.log
            exit 1
        fi
        ;;
esac

# Getting base names for all samples
# Helper to get clean base name
get_sample_basename() {
    local r1=$1
    local base=$(basename "$r1")
    base=${base%.fastq.gz}
    base=${base%.fq.gz}
    base=${base%_R1}
    base=${base%_1}
    base=${base%.R1}
    base=${base%.1}
    echo "$base" | sed 's/[^a-zA-Z0-9._-]//g'
}

get_trimmed_files() {
    local r1=$1
    local r2=$2
    local name1=$(basename "$r1")
    local name2=$(basename "$r2")
    name1=${name1%.fastq.gz}
    name1=${name1%.fq.gz}
    name2=${name2%.fastq.gz}
    name2=${name2%.fq.gz}
    echo "$ALIGNMENT_DIR/${name1}_val_1.fq.gz $ALIGNMENT_DIR/${name2}_val_2.fq.gz"
}

echo "Using genome size: $GENOME_SIZE for $GENOME_SIZE_STRING" | tee -a $LOG_DIR/pipeline.log

# Step 1: FastQC
echo "Running FastQC..." | tee -a $LOG_DIR/pipeline.log
for role in treatment control; do
    R1=$(jq -r ".samples.${role}.r1 // empty" $CONFIG_FILE)
    R2=$(jq -r ".samples.${role}.r2 // empty" $CONFIG_FILE)
    if [ -n "$R1" ] && [ -n "$R2" ]; then
        $FASTQC_PATH "$R1" "$R2" -o $OUTPUT_DIR/fastqc_reports \
            > "$LOG_DIR/fastqc_${role}.log" 2> "$LOG_DIR/fastqc_${role}_error.log"
    fi
done

# Step 2: Adapter trimming
echo "Running Trim Galore and renaming output..." | tee -a "$LOG_DIR/pipeline.log"
for role in treatment control; do
    R1=$(jq -r ".samples.${role}.r1 // empty" $CONFIG_FILE)
    R2=$(jq -r ".samples.${role}.r2 // empty" $CONFIG_FILE)

    if [ -n "$R1" ] && [ -n "$R2" ]; then
        BASE=$(get_sample_basename "$R1")

        trim_galore --paired --quality 20 --phred33 \
            --output_dir "$ALIGNMENT_DIR" \
            "$R1" "$R2" \
            > "$LOG_DIR/trim_galore_${BASE}.log" 2> "$LOG_DIR/trim_galore_${BASE}_error.log"

        # Detect output
        VAL1=$(find "$ALIGNMENT_DIR" -name "*_val_1.fq.gz" | grep "$BASE" | head -n1)
        VAL2=$(find "$ALIGNMENT_DIR" -name "*_val_2.fq.gz" | grep "$BASE" | head -n1)

        if [ "$role" == "treatment" ] && { [ ! -f "$VAL1" ] || [ ! -f "$VAL2" ]; }; then
            echo "❌ Trimmed FASTQ files for treatment not found. Aborting." | tee -a "$LOG_DIR/pipeline.log"
            exit 1
        fi

        if [ -f "$VAL1" ] && [ -f "$VAL2" ]; then
            mv "$VAL1" "$ALIGNMENT_DIR/${BASE}_trimmed_R1.fq.gz"
            mv "$VAL2" "$ALIGNMENT_DIR/${BASE}_trimmed_R2.fq.gz"
        else
            echo "⚠️  Trimmed files for optional sample '$role' not found. Skipping." | tee -a "$LOG_DIR/pipeline.log"
        fi
    fi
done


# Step 3: Align reads to the reference genome using STAR
echo "Aligning trimmed FASTQ paired-end reads to the reference genome using STAR..." | tee -a $LOG_DIR/pipeline.log

TREATMENT_BASE=$(get_sample_basename "$TREATMENT_R1")
TREATMENT_TRIMMED_R1="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R1.fq.gz"
TREATMENT_TRIMMED_R2="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R2.fq.gz"

$STAR_PATH --runThreadN "$NUM_THREADS" \
    --genomeDir "$REFERENCE_GENOME" \
    --readFilesIn "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$ALIGNMENT_DIR/${TREATMENT_BASE}." \
    --outSAMtype BAM SortedByCoordinate \
    > "$LOG_DIR/STAR_${TREATMENT_BASE}.log" 2> "$LOG_DIR/STAR_${TREATMENT_BASE}_error.log"

if [ -n "$CONTROL_R1" ] && [ -n "$CONTROL_R2" ]; then
    CONTROL_BASE=$(get_sample_basename "$CONTROL_R1")
    CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
    CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"

    $STAR_PATH --runThreadN "$NUM_THREADS" \
        --genomeDir "$REFERENCE_GENOME" \
        --readFilesIn "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$ALIGNMENT_DIR/${CONTROL_BASE}." \
        --outSAMtype BAM SortedByCoordinate \
        > "$LOG_DIR/STAR_${CONTROL_BASE}.log" 2> "$LOG_DIR/STAR_${CONTROL_BASE}_error.log"
fi

# Step 4.1: Remove duplicates using Picard's MarkDuplicates
echo "Removing duplicates using Picard MarkDuplicates..." | tee -a $LOG_DIR/pipeline.log
for bam in "$ALIGNMENT_DIR"/*.sortedByCoord.out.bam; do
    base=$(basename "$bam" .sortedByCoord.out.bam)
    java -jar "$PICARD_PATH" MarkDuplicates \
        I="$bam" \
        O="$ALIGNMENT_DIR/${base}.dedup.bam" \
        M="$LOG_DIR/${base}.metrics.txt" REMOVE_DUPLICATES=true
done

# Step 4.2: Filter by fragment size
echo "Filtering by fragment size: $FRAGMENT_SIZE_FILTER..." | tee -a $LOG_DIR/pipeline.log
for bam in "$ALIGNMENT_DIR"/*.dedup.bam; do
    base=$(basename "$bam" .dedup.bam)
    if [ "$FRAGMENT_SIZE_FILTER" == "histones" ]; then
        samtools view -h "$bam" | awk '{if ($9 >= 130 && $9 <= 300 || $1 ~ /^@/) print $0}' | samtools view -bS - > "$ALIGNMENT_DIR/${base}.filtered.bam"
    elif [ "$FRAGMENT_SIZE_FILTER" == "transcription_factors" ]; then
        samtools view -h "$bam" | awk '{if ($9 < 130 || $1 ~ /^@/) print $0}' | samtools view -bS - > "$ALIGNMENT_DIR/${base}.filtered.bam"
    else
        samtools view -h "$bam" | awk '{if ($9 < 1000 || $1 ~ /^@/) print $0}' | samtools view -bS - > "$ALIGNMENT_DIR/${base}.filtered.bam"
    fi
done

# Step 5: Peak calling with MACS2 (both broad and narrow peaks)
echo "Running MACS2 for peak calling (both broad and gapped peaks)..." | tee -a $LOG_DIR/pipeline.log
TREATMENT_FILTERED="$ALIGNMENT_DIR/${TREATMENT_BASE}.filtered.bam"
if [ -n "$CONTROL_R1" ] && [ -n "$CONTROL_R2" ]; then
    CONTROL_FILTERED="$ALIGNMENT_DIR/${CONTROL_BASE}.filtered.bam"
    macs2 callpeak -t "$TREATMENT_FILTERED" -c "$CONTROL_FILTERED" --broad --outdir "$OUTPUT_DIR/macs2_peaks" -n "$TREATMENT_BASE"
    macs2 callpeak -t "$TREATMENT_FILTERED" -c "$CONTROL_FILTERED" --call-summits --outdir "$OUTPUT_DIR/macs2_peaks" -n "$TREATMENT_BASE"
else
    macs2 callpeak -t "$TREATMENT_FILTERED" --broad --outdir "$OUTPUT_DIR/macs2_peaks" -n "$TREATMENT_BASE"
    macs2 callpeak -t "$TREATMENT_FILTERED" --call-summits --outdir "$OUTPUT_DIR/macs2_peaks" -n "$TREATMENT_BASE"
fi

# Step 6: Generate BigWig Files from BAM
echo "Generating BigWig files from BAM files..." | tee -a $LOG_DIR/pipeline.log
# BedGraph and BigWig
echo "Generating BigWig files..." | tee -a "$LOG_DIR/pipeline.log"
for bam in "$ALIGNMENT_DIR"/*.dedup.bam; do
    base=$(basename "$bam" .dedup.bam)
    bedtools genomecov -ibam "$bam" -bg > "$OUTPUT_DIR/bigwig_bedgraphs/${base}.bedgraph"
    bedGraphToBigWig "$OUTPUT_DIR/bigwig_bedgraphs/${base}.bedgraph" "$CHROM_SIZE" "$OUTPUT_DIR/bigwig_bedgraphs/${base}.bw"
done

# Step 7: Annotate Peaks (optional)
echo "Annotating peaks using bedtools..." | tee -a $LOG_DIR/pipeline.log
echo "Annotating peaks..." | tee -a "$LOG_DIR/pipeline.log"
for np in "$OUTPUT_DIR/macs2_peaks"/*.narrowPeak; do
    base=$(basename "$np" .narrowPeak)
    bedtools intersect -a "$np" -b "$ANNOTATION_GENES" -wa -wb > "$OUTPUT_DIR/annotated_peaks/${base}.annotated.bed"
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $10, $11, $12}' "$OUTPUT_DIR/annotated_peaks/${base}.annotated.bed" > "$OUTPUT_DIR/annotated_peaks/${base}.annotated.tsv"
done

# Step 8: MultiQC (aggregate QC results from both FastQC and Alignment)
echo "Generating MultiQC report for both FastQC and alignment results..." | tee -a $LOG_DIR/pipeline.log
multiqc "$OUTPUT_DIR/fastqc_reports" "$ALIGNMENT_DIR" -o "$OUTPUT_DIR/multiqc_report_combined" \
    > "$LOG_DIR/multiqc.log" 2> "$LOG_DIR/multiqc_error.log"

# Final report
echo "Analysis complete! Check the output directory for the results, including BigWig files, broad/narrow/gapped peaks, and annotated peaks in TSV format." | tee -a $LOG_DIR/pipeline.log
