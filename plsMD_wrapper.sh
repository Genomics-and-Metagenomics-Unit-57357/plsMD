#!/bin/bash
# Local wrapper to run plsMD through Docker


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run through Docker with:
# - Current directory mounted to /data
# - All arguments passed through
# - Interactive TTY if available
docker run ${TTY:+-it} \
  -v "$DIR":/data \
  -v /tmp:/tmp \
  --user $(id -u):$(id -g) \
  plsmd "$@"
