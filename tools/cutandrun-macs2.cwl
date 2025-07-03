cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["bash", "-c"]

requirements:
  DockerRequirement:
    class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

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

inputs:
  config_json:
    type: File
    doc: "Your JSON, with only relative paths (e.g. annotation/hg38.refGene.gtf)"
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

arguments:
  - |
    # run from the CWL‑staged workspace
    bash /usr/local/bin/cutrun.sh config_for_docker.json

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

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
