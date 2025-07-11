cwlVersion: v1.2
class: Workflow

inputs:
  config_json:
    type: File
  fastq_dir:
    type: Directory
  reference_genome_dir:
    type: Directory
  ecoli_index_dir:
    type: Directory
  chrom_sizes:
    type: File
  annotation_genes:
    type: File

steps:
  run_cutrun:
    run: ../tools/cutandrun-macs2.cwl
    in:
      config_json:        config_json
      fastq_dir:          fastq_dir
      reference_genome_dir: reference_genome_dir
      ecoli_index_dir:    ecoli_index_dir
      chrom_sizes:        chrom_sizes
      annotation_genes:   annotation_genes
    out:
      - alignment_replicates
      - output_replicates
      - cutrun_stdout
      - cutrun_stderr

outputs:
  alignment_replicates:
    type: Directory
    outputSource: run_cutrun/alignment_replicates

  output_replicates:
    type: Directory
    outputSource: run_cutrun/output_replicates

  log_stdout:
    type: File
    outputSource: run_cutrun/cutrun_stdout

  log_stderr:
    type: File
    outputSource: run_cutrun/cutrun_stderr
