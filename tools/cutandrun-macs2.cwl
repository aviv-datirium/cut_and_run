cwlVersion: v1.2
class: CommandLineTool

# We shell out via bash -c, pointing at the already‐inside‐container cutrun.sh.
baseCommand:
  - bash
  - -c
  - bash /usr/local/bin/cutrun.sh config_for_docker.json

requirements:
  # Use your Docker image
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  # Copy in the entire project_dir *and* the JSON as config_for_docker.json
  - class: InitialWorkDirRequirement
    listing:
      # this will replicate your cut_and_run/ tree as the CWL working dir
      - entry: $(inputs.project_dir)
        entryname: .
      # and slap your config into the same place
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

inputs:
  project_dir:
    type: Directory
    doc: |
      The root of your repository (contains fastq/, star_indices/, annotation/, etc.)
  config_json:
    type: File
    doc: |
      Your JSON‐pipeline config; will be staged as `config_for_docker.json`

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: output_replicates
    doc: The directory tree your pipeline writes (fastqc_reports/, macs2_peaks/, …)
  log_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
  log_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
