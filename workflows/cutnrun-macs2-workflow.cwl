cwlVersion: v1.2
class: Workflow

inputs:
  project_dir: Directory
  config_json: File
  fastq_dir: Directory
  reference_genome_dir: Directory
  ecoli_index_dir: Directory
  chrom_sizes: File
  annotation_genes: File

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
      project_dir: project_dir
      config_json: config_json
      fastq_dir: fastq_dir
      reference_genome_dir: reference_genome_dir
      ecoli_index_dir: ecoli_index_dir
      chrom_sizes: chrom_sizes
      annotation_genes: annotation_genes
    out: [output_dir, log_stdout, log_stderr]
