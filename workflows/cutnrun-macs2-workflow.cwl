cwlVersion: v1.2
class: Workflow

# -----------------------------------------------------------------------------
# WORKFLOW INPUTS
# -----------------------------------------------------------------------------
inputs:
  # Your JSON config (with relative paths inside the project_dir tree)
  config_json:
    type: File
    doc: |
      A JSON file listing sample R1/R2 pairs (with paths like
      "fastq/min_msto211h/…"), plus output_dir, reference paths, etc.

  # Directory of FASTQ files
  fastq_dir:
    type: Directory
    doc: "fastq/min_msto211h"

  # STAR index for the host genome (hg38)
  reference_genome_dir:
    type: Directory
    doc: "star_indices/hg38"

  # STAR index for the E. coli spike-in
  ecoli_index_dir:
    type: Directory
    doc: "star_indices/ecoli_canonical"

  # Chromosome sizes file
  chrom_sizes:
    type: File
    doc: "chrom/hg38.chrom.sizes"

  # Gene annotation GTF
  annotation_genes:
    type: File
    doc: "annotation/hg38.refGene.gtf"

# -----------------------------------------------------------------------------
# STEPS
# -----------------------------------------------------------------------------
steps:
  run_cutrun:
    run: ../tools/cutandrun-macs2.cwl
    in:
      config_json:         config_json
      fastq_dir:           fastq_dir
      reference_genome_dir: reference_genome_dir
      ecoli_index_dir:     ecoli_index_dir
      chrom_sizes:         chrom_sizes
      annotation_genes:    annotation_genes
    out: [ output_dir, log_stdout, log_stderr ]

# -----------------------------------------------------------------------------
# WORKFLOW OUTPUTS
# -----------------------------------------------------------------------------
outputs:
  # The directory with all of your pipeline outputs (bigwigs, peaks, etc.)
  cutrun_outputs:
    type: Directory
    outputSource: run_cutrun/output_dir

  # The stdout log from within the container
  log_stdout:
    type: File
    outputSource: run_cutrun/log_stdout

  # The stderr log from within the container
  log_stderr:
    type: File
    outputSource: run_cutrun/log_stderr
