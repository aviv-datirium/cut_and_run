cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}  
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # 1) launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) the JSON config only
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

baseCommand: [ bash, run.sh ]

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

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
