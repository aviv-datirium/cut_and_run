cwlVersion: v1.2
class: CommandLineTool

# We’ll need JS support for our $(...) expressions:
requirements:
  InlineJavascriptRequirement: {}

  # pull your freshly-built image
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  # only stage the things we actually write into the container’s CWD:
  InitialWorkDirRequirement:
    listing:
      # (1) a tiny launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # (2) your JSON config, under exactly this name:
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

baseCommand:
  - bash
  - run.sh

inputs:
  config_json:
    type: File
    inputBinding:
      position: 1
      prefix: ""    # no extra flag, just the path

  # everything else we mount read-only, no staging here:
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

# capture logs:
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
