cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  - class: InitialWorkDirRequirement
    listing:
      # 1) Launcher script: make the two output dirs, then call your pipeline.
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) Your config JSON, staged by copying it in.
      - class: File
        location: $(inputs.config_json.location)
        basename: config_for_docker.json

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
