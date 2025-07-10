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
          
          # Create directory structure
          mkdir -p star_indices chrom annotation
          
          # Create symlinks to the actual mounted paths in the container
          # CWL mounts directories with specific names, let's find them
          find /tmp -name "hg38" -type d | head -1 | xargs -I {} ln -sf {} star_indices/hg38
          find /tmp -name "ecoli_canonical" -type d | head -1 | xargs -I {} ln -sf {} star_indices/ecoli_canonical
          find /tmp -name "min_msto211h" -type d | head -1 | xargs -I {} ln -sf {} fastq
          find /tmp -name "hg38.chrom.sizes" -type f | head -1 | xargs -I {} ln -sf {} chrom/hg38.chrom.sizes
          find /tmp -name "hg38.refGene.gtf" -type f | head -1 | xargs -I {} ln -sf {} annotation/hg38.refGene.gtf
          
          # Run the actual pipeline
          bash /usr/local/bin/cutrun.sh config_for_docker.json
          
          # Fix permissions for cleanup
          chmod -R 755 . 2>/dev/null || true
        writable: true
      # 2) our JSON config
      - entryname: config_for_docker.json
        entry: $(inputs.config_json)
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
