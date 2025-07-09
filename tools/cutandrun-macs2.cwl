cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ bash, run.sh ]

requirements:
  - class: DockerRequirement
    dockerPull: "biowardrobe2/cutrun-macs2-core:V1.1.1"
  - class: InitialWorkDirRequirement
    listing:
      # 1) tiny wrapper script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd .                            # no more $(pwd)!
          exec bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) your config JSON
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

inputs:
  config_json:
    type: File
    doc: "Your config_for_macs2_cwl.json"
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

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
