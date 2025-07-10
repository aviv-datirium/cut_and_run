cwlVersion: v1.2
class: CommandLineTool

# turn JS expressions back on
requirements:
  - class: InlineJavascriptRequirement

  # pull the Docker image
  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  # stage exactly what you need into the container's CWD,
  # and make the two output dirs up‐front so CWL can clean them
  - class: InitialWorkDirRequirement
    listing:
      # 1) launcher script: mkdir outputs, then call cutrun.sh
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) inline your JSON config file
      # this will copy the *contents* of the File into
      # config_for_docker.json in the container
      - entry: $(inputs.config_json.contents)
        entryname: config_for_docker.json

inputs:
  config_json:
    type: File

  fastq_dir:
    type: Directory

  reference_genome_dir:
    type: Directory

  ecoli_index_dir:
    type: Directory

  chrom_sizes:
    type: File

  annotation_genes:
    type: File

# we ignore all other inputs on the command line;
# run.sh will read config_for_docker.json directly.
baseCommand: [ bash, run.sh ]

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
