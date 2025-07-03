cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["bash","-c"]
requirements:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

inputs:
  # Mount your entire project here
  project_dir:
    type: Directory
    doc: "Your project root containing fastq/, star_indices/, dockerfiles/, etc."

  # Just the one config file
  config_json:
    type: File
    doc: "The config_for_macs2_cwl.json to drive cutrun.sh"
    # stage it into the workdir as config_for_docker.json
requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

outputs:
  # after the run, all of your outputs live under project/output_replicates
  output_dir:
    type: Directory
    outputBinding:
      glob: project/output_replicates
  log_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
  log_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

arguments:
  - |
    # cd into the mounted project, then run the pipeline
    cd "$(inputs.project_dir.path)" && \
    bash /usr/local/bin/cutrun.sh config_for_docker.json
