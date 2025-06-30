# cutrun_workflow.cwl — wraps the single-tool CUT&RUN pipeline tool as a workflow

cwlVersion: v1.2
class: Workflow

inputs:
  config_json: File

outputs:
  cutrun_outputs:
    type: Directory
    outputSource: run_cutrun/output_dir

steps:
  run_cutrun:
    run: /mnt/data/home/aviv/cut_and_run/CommandLineTools/cutandrun-macs2.cwl
    in:
      config_json: config_json
    out: [output_dir]
