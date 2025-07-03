# cutrun_workflow.cwl — wraps the single-tool CUT&RUN pipeline tool as a workflow

cwlVersion: v1.2
class: Workflow

inputs:
  config_json: File
  annotation_genes: File
  chrom_sizes: File
  ecoli_index_dir: Directory
  reference_genome_dir: Directory
  fastq_dir: Directory

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

steps:
  run_cutrun:
    run: ../tools/cutandrun-macs2.cwl
    in:
      config_json: config_json
      annotation_genes: annotation_genes
      chrom_sizes: chrom_sizes
      ecoli_index_dir: ecoli_index_dir
      reference_genome_dir: reference_genome_dir
      fastq_dir: fastq_dir
    out: [output_dir, log_stdout, log_stderr]
