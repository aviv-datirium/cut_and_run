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
  ResourceRequirement:
    coresMin: 1
    ramMin: 256
  InitialWorkDirRequirement:
    listing:
      # 1) our tiny wrapper that makes the directories & symlinks
      - entryname: run.sh
        writable: true
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail

          # recreate the layout the pipeline expects
          mkdir -p star_indices/hg38 \
                   star_indices/ecoli_canonical \
                   fastq \
                   chrom \
                   annotation

          # point those dirs at the real mounts
          ln -sf "$(inputs.reference_genome_dir.path)" star_indices/hg38
          ln -sf "$(inputs.ecoli_index_dir.path)"    star_indices/ecoli_canonical
          ln -sf "$(inputs.fastq_dir.path)"          fastq
          ln -sf "$(inputs.chrom_sizes.path)"        chrom/hg38.chrom.sizes
          ln -sf "$(inputs.annotation_genes.path)"   annotation/hg38.refGene.gtf

          # now run the real entrypoint
          bash /usr/local/bin/cutrun.sh config_for_docker.json

          # relax perms so CWL can clean up
          chmod -R a+rwx . || true

      # 2) your JSON config
      - entryname: config_for_docker.json
        entry: $(inputs.config_json.path)


      # 2) the config JSON
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # 3) the fastq directory
      - entry: $(inputs.fastq_dir)
        entryname: fastq

      # 4) human genome index
      - entry: $(inputs.reference_genome_dir)
        entryname: star_indices/hg38

      # 5) E. coli index
      - entry: $(inputs.ecoli_index_dir)
        entryname: star_indices/ecoli_canonical

      # 6) chromosome sizes
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes

      # 7) annotation GTF
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

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
