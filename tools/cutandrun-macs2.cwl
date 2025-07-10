cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - bash
  - run.sh
  - config_for_docker.json
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
    dockerOutputDirectory: /tmp
  ResourceRequirement:
    coresMin: 1
    ramMin: 256
  InitialWorkDirRequirement:
    listing:
      # 1) wrapper script
      - entryname: run.sh
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        writable: true
      # 2) our JSON config
      - entryname: config_for_docker.json
        entry: $(inputs.config_json)
      # 3) fastqs
      - entryname: fastq
        entry: $(inputs.fastq_dir)
      # 4) hg38 STAR index
      - entryname: star_indices/hg38
        entry: $(inputs.reference_genome_dir)
      # 5) E.coli STAR index
      - entryname: star_indices/ecoli_canonical
        entry: $(inputs.ecoli_index_dir)
      # 6) chrom sizes
      - entryname: chrom/hg38.chrom.sizes
        entry: $(inputs.chrom_sizes)
      # 7) gene annotation
      - entryname: annotation/hg38.refGene.gtf
        entry: $(inputs.annotation_genes)
stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
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
outputs:
  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates
  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates
