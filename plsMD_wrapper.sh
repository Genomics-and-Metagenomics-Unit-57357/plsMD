#!/bin/bash
# Local wrapper to run plsMD through Docker

CURRENT_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=== plsMD Wrapper Debug ===" >&2
echo "Script directory: $SCRIPT_DIR" >&2
echo "Current directory: $CURRENT_DIR" >&2
echo "Arguments: $@" >&2


INPUT_DIR=""
ARGS=()
FOUND_DIR_FLAG=false


for arg in "$@"; do
    if [[ "$FOUND_DIR_FLAG" == true ]]; then
        INPUT_DIR="$arg"
        ARGS+=("$arg")
        FOUND_DIR_FLAG=false
    elif [[ "$arg" == "--dir" ]]; then
        FOUND_DIR_FLAG=true
        ARGS+=("$arg")
    elif [[ "$arg" == --dir=* ]]; then
        INPUT_DIR="${arg#*=}"
        ARGS+=("$arg")
    else
        ARGS+=("$arg")
    fi
done

echo "Detected input directory: $INPUT_DIR" >&2


if [[ -n "$INPUT_DIR" ]]; then
    if [[ ! "$INPUT_DIR" =~ ^/ ]]; then
        # Convert relative path to absolute
        INPUT_DIR="$(realpath "$INPUT_DIR" 2>/dev/null || echo "$CURRENT_DIR/$INPUT_DIR")"
    fi
    echo "Absolute input directory: $INPUT_DIR" >&2
    
    # Check if directory exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Directory $INPUT_DIR does not exist!" >&2
        exit 1
    fi
fi

VOLUMES="-v $CURRENT_DIR:/data"
VOLUMES="$VOLUMES -v $SCRIPT_DIR:/opt/plsMD/scripts"

VOLUMES="$VOLUMES -v /tmp:/tmp"
VOLUMES="$VOLUMES -v /media/db:/media/db"

if [[ -n "$INPUT_DIR" ]]; then
    VOLUMES="$VOLUMES -v $INPUT_DIR:$INPUT_DIR"
    echo "Mounting input directory: $INPUT_DIR" >&2
fi

echo "Volume mounts: $VOLUMES" >&2
echo "==========================" >&2

# Run through Docker
docker run ${TTY:+-it} \
  $VOLUMES \
  --user $(id -u):$(id -g) \
  -w /data \
  plsmd "${ARGS[@]}"
