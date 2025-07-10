cwlVersion: v1.2
class: CommandLineTool

requirements:

  # Allow $(...) JS expressions
  InlineJavascriptRequirement: {}

  # Pull the image you built
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  # Stage just exactly what your script needs in its CWD
  InitialWorkDirRequirement:
    listing:

      # 1) The launcher that calls your pipeline with the right filename
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) Your config JSON, *renamed* to match the script’s argument
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

baseCommand: [ bash, run.sh ]

inputs:
  config_json:
    type: File
    doc: |
      A JSON file containing all of the parameters your cutnrun.sh needs.

  # All other inputs are bind-mounted read-only; no staging here:
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

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
