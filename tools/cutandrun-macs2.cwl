cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest
  InitialWorkDirRequirement:
    listing:
      # 1) launcher
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) the JSON config
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json
        writable: false

      # 3) fastqs
      - entry: $(inputs.fastq_dir.path)
        entryname: fastq
        writable: false

      # 4) hg38 STAR index
      - entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38
        writable: false

      # 5) E. coli STAR index
      - entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical
        writable: false

      # 6) chromosome sizes
      - entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes
        writable: false

      # 7) gene annotation
      - entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf
        writable: false

baseCommand:
  - bash
  - run.sh

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
