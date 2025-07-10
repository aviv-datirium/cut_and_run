cwlVersion: v1.2
class: CommandLineTool

requirements:
  # enable JS expressions
  InlineJavascriptRequirement: {}

  # pull the Docker image
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  # stage only what we need into the container’s CWD
  InitialWorkDirRequirement:
    listing:
      # 1) the wrapper script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config json
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) input dirs & files
      - entry: $(inputs.fastq_dir.path)
        entryname: fastq
      - entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38
      - entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical
      - entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

  ResourceRequirement:
    coresMin: 1
    ramMin: 4096

baseCommand: ["bash", "run.sh"]

inputs:
  config_json:
    type: File
    inputBinding: {}

  fastq_dir:
    type: Directory
    inputBinding: {}

  reference_genome_dir:
    type: Directory
    inputBinding: {}

  ecoli_index_dir:
    type: Directory
    inputBinding: {}

  chrom_sizes:
    type: File
    inputBinding: {}

  annotation_genes:
    type: File
    inputBinding: {}

outputs:
  cutrun_stdout:
    type: stdout

  cutrun_stderr:
    type: stderr

  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates
