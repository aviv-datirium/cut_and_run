# cutrun_workflow.cwl — wraps the single-tool CUT&RUN pipeline tool as a workflow

cwlVersion: v1.2
class: Workflow

inputs:
  config_json: File
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
      config_json: config_json
    out: [output_dir, log_stdout, log_stderr]
