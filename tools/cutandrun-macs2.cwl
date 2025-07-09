cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ bash, run.sh ]

requirements:
  - class: DockerRequirement
    dockerPull: "biowardrobe2/cutrun-macs2-core:v1.1.1"
  - class: InitialWorkDirRequirement
    listing:
      # 1) the tiny wrapper script—no $(pwd) anywhere!
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd .                              # just stay in the staging dir
          exec bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) your config JSON
      - entry: $(inputs.config_json)
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

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: output_replicates

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
