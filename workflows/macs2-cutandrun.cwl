requirements:
  DockerRequirement:
    dockerPull: ghcr.io/<org>/cutrun-core:latest   # full registry/image:tag

baseCommand: [/usr/local/bin/cutrun.sh]

inputs:
  config:
    type: File
    inputBinding: {}
  # add other Files/Directories if needed

outputs: {}
