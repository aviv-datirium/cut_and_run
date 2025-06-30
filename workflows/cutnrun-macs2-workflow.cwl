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
    run: cutrun.cwl
    in:
      config_json: config_json
    out: [output_dir]
