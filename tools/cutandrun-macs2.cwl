# tools/cutandrun-seacr.cwl
cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: cutrun-seacr-core:latest

  InitialWorkDirRequirement:
    listing:
      # 1) your SEACR config JSON
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # 2) raw FASTQs
      - entry: $(inputs.fastq_dir)
        entryname: fastq

      # 3) genome indices (hg38 & ecoli)
      - entry: $(inputs.star_indices)
        entryname: star_indices

      # 4) chrom sizes
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes

      # 5) annotation GTF
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

baseCommand:
  # cutrun-seacr-core image ENTRYPOINT expects the JSON path
  - /usr/local/bin/cutrun.sh
  - config_for_docker.json

inputs:
  config_json:
    type: File

  fastq_dir:
    type: Directory

  star_indices:
    type: Directory

  chrom_sizes:
    type: File

  annotation_genes:
    type: File

outputs:
  alignment_seacr:
    type: Directory
    outputBinding:
      glob:
        - alignment_seacr
        - alignment_seacr/**

  output_seacr:
    type: Directory
    outputBinding:
      glob:
        - output_seacr
        - output_seacr/**
