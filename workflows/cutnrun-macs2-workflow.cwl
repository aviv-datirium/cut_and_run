cwlVersion: v1.2
class: Workflow

inputs:
  config_json:
    type: File
    doc: |
      Your JSON config (with relative paths inside the staged workspace,
      e.g. "fastq/min_msto211h/…", etc.)

  fastq_dir:
    type: Directory
    doc: "Directory of FASTQ files (e.g. fastq/min_msto211h)"

  reference_genome_dir:
    type: Directory
    doc: "STAR index for the host genome (e.g. star_indices/hg38)"

  ecoli_index_dir:
    type: Directory
    doc: "STAR index for the E. coli spike-in (e.g. star_indices/ecoli_canonical)"

  chrom_sizes:
    type: File
    doc: "Chromosome sizes file (e.g. chrom/hg38.chrom.sizes)"

  annotation_genes:
    type: File
    doc: "Gene annotation GTF (e.g. annotation/hg38.refGene.gtf)"

steps:
  run_cutrun:
    run: ../tools/cutandrun-macs2.cwl
    in:
      config_json:          config_json
      fastq_dir:            fastq_dir
      reference_genome_dir: reference_genome_dir
      ecoli_index_dir:      ecoli_index_dir
      chrom_sizes:          chrom_sizes
      annotation_genes:     annotation_genes
    out: [ output_dir, log_stdout, log_stderr ]

outputs:
  cutrun_outputs:
    type: Directory
    outputSource: run_cutrun/output_dir

  log_stdout:
    type: File
    outputSource: run_cutrun/log_stdout

  log_stderr:
    type: File
    outputSource: run_cutrun/log_stderr
