cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # 1) small launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          # run.sh and config_for_docker.json are in CWD
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) your JSON config
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

baseCommand: [ bash, run.sh ]

inputs:
  config_json:
    type: File
    inputBinding:                        # still needs to be on the command line
      position: 1

  # everything else will be bind-mounted by CWL under their staging paths:
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
