cwlVersion: v1.2
class: CommandLineTool

# call bash with the pipeline script
baseCommand:
  - bash
  - /usr/local/bin/cutrun.sh

# mount everything into the container
requirements:
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # 1) your config JSON, renamed in‐place
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 2) data directories & files
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
    # no prefix, at argv position 2 (after bash, script)
    inputBinding:
      position: 2
      prefix: ""

  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

# capture stdout/stderr
stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

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
