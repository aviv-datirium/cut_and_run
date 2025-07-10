cwlVersion: v1.2
class: CommandLineTool

requirements:
  # 1) run inside your freshly built image
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  # 2) drop into the container only the bits
  InitialWorkDirRequirement:
    listing:
      # (a) our tiny wrapper that just execs cutrun.sh
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          # invoke the "real" pipeline
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh

      # (b) stage in YOUR config JSON under the name it expects:
      - class: File
        path: $(inputs.config_json.path)
        basename: config_for_docker.json

# we DON'T stage the dirs manually here; CWL will bind-mount them
baseCommand: [ bash, run.sh ]

inputs:
  # the only thing we copy in is the config file
  config_json:
    type: File

  # these six inputs CWL will bind-mount read-only
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
