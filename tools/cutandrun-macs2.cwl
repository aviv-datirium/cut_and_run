cwlVersion: v1.2
class: CommandLineTool

# Enable JS expressions and pull the Docker image
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      # 1) launcher script
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config JSON
      - class: File
        location: $(inputs.config_json.path)
        basename: config_for_docker.json

      # 3) fastq directory
      - class: Directory
        location: $(inputs.fastq_dir.path)
        basename: fastq

      # 4) reference genome (hg38)
      - class: Directory
        location: $(inputs.reference_genome_dir.path)
        basename: star_indices/hg38

      # 5) E. coli index
      - class: Directory
        location: $(inputs.ecoli_index_dir.path)
        basename: star_indices/ecoli_canonical

      # 6) chromosome sizes
      - class: File
        location: $(inputs.chrom_sizes.path)
        basename: chrom/hg38.chrom.sizes

      # 7) gene annotation
      - class: File
        location: $(inputs.annotation_genes.path)
        basename: annotation/hg38.refGene.gtf

baseCommand: [ bash, run.sh ]

inputs:
  config_json:
    type: File
    doc: "Your MACS2 Cut&Run config JSON"
  fastq_dir:
    type: Directory
    doc: "Directory containing your FASTQ files"
  reference_genome_dir:
    type: Directory
    doc: "STAR index dir for hg38"
  ecoli_index_dir:
    type: Directory
    doc: "STAR index dir for *E. coli*"
  chrom_sizes:
    type: File
    doc: "Chromosome sizes file"
  annotation_genes:
    type: File
    doc: "Gene annotation GTF"

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
