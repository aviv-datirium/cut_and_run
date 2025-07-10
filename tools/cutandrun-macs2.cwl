cwlVersion: v1.2
class: CommandLineTool

baseCommand: bash
arguments:
  - run.sh

requirements:
  - class: InlineJavascriptRequirement

  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  - class: InitialWorkDirRequirement
    listing:
      # 1) our tiny launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) bring in your JSON as a File object
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
        writable: false

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

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
