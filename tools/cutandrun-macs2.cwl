cwlVersion: v1.2
class: CommandLineTool

# we invoke bash on our staged run.sh
baseCommand: bash
arguments:
  - run.sh

requirements:
  # allow JS expressions
  - class: InlineJavascriptRequirement

  # pull your image
  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  # stage exactly what you need into the container CWD
  - class: InitialWorkDirRequirement
    listing:
      # 1) launcher script: make the output dirs and call cutrun.sh
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) your JSON config, as a CWL File object
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
        writable: false

outputs:
  # these two directories must now exist in the workdir after run.sh
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  # capture container logs
  cutrun_stdout:
    type: File
    outputBinding:
      glob: $(self.basename + ".stdout.log")
  cutrun_stderr:
    type: File
    outputBinding:
      glob: $(self.basename + ".stderr.log")

# map your CWL inputs into File/Directory objects
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
