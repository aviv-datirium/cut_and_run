cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:v1.1.1
    dockerImageId: biowardrobe2/cutrun-macs2-core:v1.1.1
  InitialWorkDirRequirement:
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

baseCommand:
  - bash
  - -lc
  - |
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate cutrun && \
    exec /usr/local/bin/cutrun.sh config_for_docker.json

outputs:
  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  cutrun_stdout:
    type: stdout

  cutrun_stderr:
    type: stderr
