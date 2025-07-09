cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - run.sh

requirements:
  DockerRequirement:
    dockerPull: "biowardrobe2/cutrun-macs2-core:v1.1.1"
    # overwrite the image’s ENTRYPOINT so it won’t try to do conda activate
    dockerRunOptions:
      - --entrypoint=/bin/bash

  InitialWorkDirRequirement:
    class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
      - entry: $(inputs.fastq_dir)
        entryname: fastq
      - entry: $(inputs.reference_genome_dir)
        entryname: star_indices/hg38
      - entry: $(inputs.ecoli_index_dir)
        entryname: star_indices/ecoli_canonical
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          exec bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

inputs:
  config_json:
    type: File
  fastq_dir: Directory
  reference_genome_dir: Directory
  ecoli_index_dir: Directory
  chrom_sizes: File
  annotation_genes: File

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
