cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - -c
  - |
    cd project && bash /usr/local/bin/cutrun.sh config_for_docker.json

requirements:
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.project_dir)
        entryname: project
      - entry: $(inputs.config_json)
        entryname: project/config_for_docker.json

inputs:
  project_dir:
    type: Directory
    doc: |
      Your repo root (contains fastq/, annotation/, star_indices/, etc.)
  config_json:
    type: File
    doc: |
      JSON config file with relative paths — will be staged as
      `project/config_for_docker.json`.

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: project/output_replicates
    doc: Main output tree (fastqc_reports/, macs2_peaks/, …)

  log_stdout:
    type: File
    outputBinding:
      glob: project/cutrun_stdout.log
    doc: Pipeline stdout log

  log_stderr:
    type: File
    outputBinding:
      glob: project/cutrun_stderr.log
    doc: Pipeline stderr log

stdout: project/cutrun_stdout.log
stderr: project/cutrun_stderr.log
