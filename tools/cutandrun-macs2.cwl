cwlVersion: v1.2
class: CommandLineTool

# 1) How we invoke your pipeline:
baseCommand:
  - bash
  - run.sh

# 2) We still need JS for CWL expressions, and we’ll only stage the wrapper + config.json
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # a tiny launcher script we write into the container
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # stage the config.json under that exact name:
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

inputs:
  # the only input we *write* into the workdir:
  config_json:
    type: File
    inputBinding:
      position: 1
      prefix: ""

  # the rest of your inputs are just mounted by CWL — no listing in IWDR:
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

# capture pipelines stdout/err
stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

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
