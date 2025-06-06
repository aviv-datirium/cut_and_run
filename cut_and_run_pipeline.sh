#!/bin/bash

# Read the config.json file to get input parameters
CONFIG_FILE="/mnt/data/home/aviv/bash_scripts/config.json"

RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir' $CONFIG_FILE)
ALIGNMENT_DIR=$(jq -r '.alignment_dir' $CONFIG_FILE)
OUTPUT_DIR=$(jq -r '.output_dir' $CONFIG_FILE)
REFERENCE_GENOME=$(jq -r '.reference_genome' $CONFIG_FILE)
ANNOTATION_GENES=$(jq -r '.annotation_genes' $CONFIG_FILE)
GENOME_SIZE_STRING=$(jq -r '.genome_size' $CONFIG_FILE)
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' $CONFIG_FILE)
CUSTOM_GENOME_SIZE=$(jq -r '.custom_genome_size' $CONFIG_FILE)
LOG_DIR=$(jq -r '.log_dir' $CONFIG_FILE)

# Define genome sizes for various species (numeric values in base pairs)
GENOME_SIZE_HUMAN=2913022398  # Human genome size (hg38)
GENOME_SIZE_MOUSE=2652783500  # Mouse genome size (mm10)
GENOME_SIZE_DROSOPHILA=165000000  # Drosophila melanogaster (fruit fly) genome size
GENOME_SIZE_CELEGANS=1000000000  # Caenorhabditis elegans (nematode) genome size
GENOME_SIZE_YEAST=12000000  # Saccharomyces cerevisiae (yeast) genome size

# Picard tools path
PICARD_PATH="/mnt/data/home/aviv/tools/picard.jar"  # Path to the Picard jar file (e.g., picard.jar)

# FastQC tools path
FASTQC_PATH="/mnt/data/home/aviv/tools/FastQC/fastqc"

# Make output directories
mkdir -p $OUTPUT_DIR $ALIGNMENT_DIR

# Create FastQC output directory if it doesn't exist
mkdir -p $OUTPUT_DIR/fastqc_reports

# Define log and error files
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
    "hs")
        GENOME_SIZE=$GENOME_SIZE_HUMAN
        ;;
    "mm")
        GENOME_SIZE=$GENOME_SIZE_MOUSE
        ;;
    "dm")
        GENOME_SIZE=$GENOME_SIZE_DROSOPHILA
        ;;
    "ce")
        GENOME_SIZE=$GENOME_SIZE_CELEGANS
        ;;
    "sc")
        GENOME_SIZE=$GENOME_SIZE_YEAST
        ;;
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

echo "Using genome size: $GENOME_SIZE for $GENOME_SIZE_STRING" | tee -a $LOG_DIR/pipeline.log

# Step 1: Quality Control (FastQC)
#~ echo "Running FastQC on raw FASTQ files..." | tee -a $LOG_DIR/pipeline.log
#~ $FASTQC_PATH $RAW_FASTQ_DIR/*.{fastq,fq}.gz -o $OUTPUT_DIR/fastqc_reports > $LOG_DIR/fastqc_output.log 2> $LOG_DIR/fastqc_error.log
#~ if [ $? -ne 0 ]; then
    #~ echo "FastQC failed. Check $LOG_DIR/fastqc_error.log for details." | tee -a $LOG_DIR/pipeline.log
    #~ exit 1
#~ fi

# Step 2: Adapter trimming (optional, if necessary)
echo "Trimming adapters and low-quality reads..." | tee -a $LOG_DIR/pipeline.log

# Use find to list all .fq.gz and .fastq.gz files correctly
for fastq_file in $(find $RAW_FASTQ_DIR -type f \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \)); do
    # Extract the base name by removing the extensions (.fq.gz or .fastq.gz)
    base_name=$(basename "$fastq_file" .gz)  # Remove .gz first
    base_name=${base_name%.fastq}  # Remove .fastq extension
    base_name=${base_name%.fq}  # Remove .fq extension

    # Ensure base_name doesn't have issues with the file name
    base_name=$(echo "$base_name" | sed 's/[^a-zA-Z0-9_-]//g')  # Remove any special characters

    # Trim Galore with output redirection
    trim_galore --quality 20 --phred33 --output_dir "$ALIGNMENT_DIR" "$fastq_file" > "$LOG_DIR/trim_galore_${base_name}.log" 2> "$LOG_DIR/trim_galore_${base_name}_error.log"

    if [ $? -ne 0 ]; then
        echo "Trim Galore failed for $base_name. Check $LOG_DIR/trim_galore_${base_name}_error.log for details." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi
done

# Step 3: Align reads to the reference genome using STAR
echo "Aligning reads to the reference genome using STAR..." | tee -a $LOG_DIR/pipeline.log
for trimmed_file in $ALIGNMENT_DIR/*.{fastq,fq}.gz
do
    base_name=$(basename $trimmed_file .fastq)
    base_name=${base_name%.fq}  # Handle both fastq and fq extensions
    /mnt/data/home/aviv/tools/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 8 \
         --genomeDir $REFERENCE_GENOME \
         --readFilesIn $trimmed_file \
         --outFileNamePrefix $ALIGNMENT_DIR/$base_name. \
         --outSAMtype BAM SortedByCoordinate \
         > $LOG_DIR/STAR_${base_name}.log 2> $LOG_DIR/STAR_${base_name}_error.log
    if [ $? -ne 0 ]; then
        echo "STAR alignment failed for $base_name. Check $LOG_DIR/STAR_${base_name}_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi
done

# Step 4.1: Remove duplicates using Picard's MarkDuplicates
echo "Removing duplicates using Picard MarkDuplicates..." | tee -a $LOG_DIR/pipeline.log
for sorted_bam in $ALIGNMENT_DIR/*.sortedByCoordinate.bam
do
    base_name=$(basename $sorted_bam .sortedByCoordinate.bam)
    java -jar $PICARD_PATH MarkDuplicates I=$sorted_bam O=$ALIGNMENT_DIR/$base_name.dedup.bam M=$LOG_DIR/$base_name.metrics.txt REMOVE_DUPLICATES=true
    if [ $? -ne 0 ]; then
        echo "Picard MarkDuplicates failed for $base_name. Check $LOG_DIR/$base_name.metrics.txt for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi
done

# Step 4.2: Filter by fragment size
echo "Filtering by fragment size: $FRAGMENT_SIZE_FILTER..." | tee -a $LOG_DIR/pipeline.log
for dedup_bam in $ALIGNMENT_DIR/*.dedup.bam
do
    base_name=$(basename $dedup_bam .dedup.bam)
    if [ "$FRAGMENT_SIZE_FILTER" == "histones" ]; then
        samtools view -h $dedup_bam | awk 'BEGIN{OFS="\t"} {if ($9 >= 130 && $9 <= 300) print $0}' | samtools view -bS - > $ALIGNMENT_DIR/$base_name.filtered.bam
    elif [ "$FRAGMENT_SIZE_FILTER" == "transcription_factors" ]; then
        samtools view -h $dedup_bam | awk 'BEGIN{OFS="\t"} {if ($9 < 130) print $0}' | samtools view -bS - > $ALIGNMENT_DIR/$base_name.filtered.bam
    else
        # Default: Filter fragments below 1000 bp
        samtools view -h $dedup_bam | awk 'BEGIN{OFS="\t"} {if ($9 < 1000) print $0}' | samtools view -bS - > $ALIGNMENT_DIR/$base_name.filtered.bam
    fi
done

# Step 5: Peak calling with MACS2 (both broad, narrow, and gapped peaks)
echo "Running MACS2 for peak calling (both broad, narrow, and gapped peaks)..." | tee -a $LOG_DIR/pipeline.log
for filtered_bam in $ALIGNMENT_DIR/*.filtered.bam
do
    base_name=$(basename $filtered_bam .filtered.bam)

    # Narrow peaks (for TFs, etc.)
    macs2 callpeak -t $filtered_bam -f BAM -g $GENOME_SIZE -n $base_name --outdir $OUTPUT_DIR/macs2_peaks --call-summits > $LOG_DIR/macs2_${base_name}_narrow.log 2> $LOG_DIR/macs2_${base_name}_narrow_error.log
    if [ $? -ne 0 ]; then
        echo "MACS2 narrow peak calling failed for $base_name. Check $LOG_DIR/macs2_${base_name}_narrow_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi

    # Broad peaks (for histones, etc.)
    macs2 callpeak -t $filtered_bam -f BAM -g $GENOME_SIZE -n $base_name --outdir $OUTPUT_DIR/macs2_peaks --broad > $LOG_DIR/macs2_${base_name}_broad.log 2> $LOG_DIR/macs2_${base_name}_broad_error.log
    if [ $? -ne 0 ]; then
        echo "MACS2 broad peak calling failed for $base_name. Check $LOG_DIR/macs2_${base_name}_broad_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi

    # Gapped peaks
    macs2 callpeak -t $filtered_bam -f BAM -g $GENOME_SIZE -n $base_name --outdir $OUTPUT_DIR/macs2_peaks --gappedPeak > $LOG_DIR/macs2_${base_name}_gapped.log 2> $LOG_DIR/macs2_${base_name}_gapped_error.log
    if [ $? -ne 0 ]; then
        echo "MACS2 gapped peak calling failed for $base_name. Check $LOG_DIR/macs2_${base_name}_gapped_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi
done

# Step 6: Generate BigWig Files from BAM
echo "Generating BigWig files from BAM files..." | tee -a $LOG_DIR/pipeline.log
for dedup_bam in $ALIGNMENT_DIR/*.dedup.bam
do
    base_name=$(basename $dedup_bam .dedup.bam)

    # Step 5.1: Generate BedGraph from BAM using bedtools genomecov
    bedtools genomecov -ibam $dedup_bam -g $GENOME_SIZE -bg > $OUTPUT_DIR/bigwig/$base_name.bedgraph
    if [ $? -ne 0 ]; then
        echo "BedGraph generation failed for $base_name. Check $LOG_DIR/bedgraph_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi

    # Step 5.2: Convert BedGraph to BigWig using bedGraphToBigWig (UCSC tool)
    bedGraphToBigWig $OUTPUT_DIR/bigwig/$base_name.bedgraph $GENOME_SIZE $OUTPUT_DIR/bigwig/$base_name.bw
    if [ $? -ne 0 ]; then
        echo "BigWig conversion failed for $base_name. Check $LOG_DIR/bigwig_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi
done

# Step 7: Annotate Peaks (optional)
echo "Annotating peaks using bedtools..." | tee -a $LOG_DIR/pipeline.log
for narrowpeak_file in $OUTPUT_DIR/macs2_peaks/*.narrowPeak
do
    base_name=$(basename $narrowpeak_file .narrowPeak)
    bedtools intersect -a $narrowpeak_file -b $ANNOTATION_GENES -wa -wb > $OUTPUT_DIR/annotated_peaks/$base_name.annotated.bed
    if [ $? -ne 0 ]; then
        echo "Peak annotation failed for $base_name. Check $LOG_DIR/annotated_peaks_error.log for details." | tee -a $LOG_DIR/pipeline.log
        exit 1
    fi
    # Convert annotated peaks to TSV format
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $10, $11, $12}' $OUTPUT_DIR/annotated_peaks/$base_name.annotated.bed > $OUTPUT_DIR/annotated_peaks/$base_name.annotated.tsv
done

# Step 8: MultiQC (aggregate QC results from both FastQC and Alignment)
echo "Generating MultiQC report for both FastQC and alignment results..." | tee -a $LOG_DIR/pipeline.log
multiqc $OUTPUT_DIR/fastqc_reports $ALIGNMENT_DIR/ -o $OUTPUT_DIR/multiqc_report_combined > $LOG_DIR/multiqc_combined_output.log 2> $LOG_DIR/multiqc_combined_error.log
if [ $? -ne 0 ]; then
    echo "MultiQC failed for FastQC and alignment results. Check $LOG_DIR/multiqc_combined_error.log for details." | tee -a $LOG_DIR/pipeline.log
    exit 1
fi

# Final report
echo "Analysis complete! Check the output directory for the results, including BigWig files, broad/narrow/gapped peaks, and annotated peaks in TSV format." | tee -a $LOG_DIR/pipeline.log
