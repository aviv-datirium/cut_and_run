cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - -c

requirements:
  DockerRequirement:
    class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
      - entry: fastq/min_msto211h
        entryname: fastq
      - entry: star_indices/hg38
        entryname: star_indices/hg38
      - entry: star_indices/ecoli_canonical
        entryname: star_indices/ecoli_canonical
      - entry: chrom/hg38.chrom.sizes
        entryname: chrom/hg38.chrom.sizes
      - entry: annotation/hg38.refGene.gtf
        entryname: annotation/hg38.refGene.gtf

arguments:
  # after staging everything exactly where your JSON expects it...
  - |
    bash /usr/local/bin/cutrun.sh config_for_docker.json

inputs:
  config_json:
    type: File
    doc: "Your JSON, with only relative paths (e.g. annotation/hg38.refGene.gtf)"

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
