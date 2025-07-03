cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest
  - class: InitialWorkDirRequirement
    listing:
      # mount your whole repo under `project/`
      - entry: $(inputs.project_dir)
        entryname: project
      # drop in the JSON as project/config_for_docker.json
      - entry: $(inputs.config_json)
        entryname: project/config_for_docker.json

baseCommand: ["bash","-c"]
arguments:
  - |
    cd project && \
    bash /usr/local/bin/cutrun.sh config_for_docker.json \
      > cutrun_stdout.log 2> cutrun_stderr.log

inputs:
  project_dir:
    type: Directory
    doc: your top-level cut_and_run/ working dir
  config_json:
    type: File
    doc: the JSON you’ve already been editing

outputs:
  # your main output tree
  output_dir:
    type: Directory
    outputBinding:
      glob: project/output_replicates

  # the two logs we explicitly redirected
  log_stdout:
    type: File
    outputBinding:
      glob: project/cutrun_stdout.log

  log_stderr:
    type: File
    outputBinding:
      glob: project/cutrun_stderr.log
