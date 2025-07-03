cwlVersion: v1.2
class: Workflow

inputs:
  config_json: File

steps:
  run_cutrun:
    run: ../tools/cutandrun-macs2.cwl
    in:
      config_json: config_json
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
