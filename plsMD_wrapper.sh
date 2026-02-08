#!/bin/bash
# plsMD Docker Wrapper v1.0
# Handles Docker execution with proper volume mounting and path resolution

set -e

CURRENT_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Color output for better visibility
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Debug mode (set DEBUG=1 to enable)
DEBUG=${DEBUG:-0}

debug_log() {
    if [ "$DEBUG" -eq 1 ]; then
        echo -e "${YELLOW}[DEBUG]${NC} $1" >&2
    fi
}

error_log() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

info_log() {
    echo -e "${GREEN}[INFO]${NC} $1" >&2
}

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    error_log "Docker is not running or not installed"
    exit 1
fi

# Check if plsmd image exists
if ! docker image inspect plsmd >/dev/null 2>&1; then
    error_log "plsmd Docker image not found. Please build it first:"
    echo "  docker build -t plsmd ." >&2
    exit 1
fi

debug_log "Script directory: $SCRIPT_DIR"
debug_log "Current directory: $CURRENT_DIR"
debug_log "Arguments: $*"

# Parse arguments and collect directories to mount
VOLUMES=()
ARGS=()
DIRS_TO_MOUNT=()

# Always mount current directory
VOLUMES+=("-v" "$CURRENT_DIR:/data")
debug_log "Mounting current directory: $CURRENT_DIR -> /data"

# Parse arguments to find directories
i=0
for arg in "$@"; do
    # Check for --dir flag
    if [[ "$arg" == "--dir" ]]; then
        i=$((i + 1))
        next_arg="${@:$((i+2)):1}"
        if [[ -n "$next_arg" && ! "$next_arg" =~ ^-- ]]; then
            DIRS_TO_MOUNT+=("$next_arg")
        fi
    elif [[ "$arg" == --dir=* ]]; then
        dir="${arg#*=}"
        DIRS_TO_MOUNT+=("$dir")
    fi
    
    # Check for other path-containing arguments
    if [[ "$arg" == --input=* || "$arg" == --output=* || "$arg" == --path=* ]]; then
        path="${arg#*=}"
        if [[ -e "$path" ]]; then
            dir=$(dirname "$path")
            DIRS_TO_MOUNT+=("$dir")
        fi
    elif [[ "$arg" == "--input" || "$arg" == "--output" || "$arg" == "--path" ]]; then
        i=$((i + 1))
        next_arg="${@:$((i+2)):1}"
        if [[ -n "$next_arg" && -e "$next_arg" ]]; then
            dir=$(dirname "$next_arg")
            DIRS_TO_MOUNT+=("$dir")
        fi
    fi
    
    i=$((i + 1))
    ARGS+=("$arg")
done

# Mount all discovered directories
for dir in "${DIRS_TO_MOUNT[@]}"; do
    # Convert to absolute path
    if [[ ! "$dir" =~ ^/ ]]; then
        abs_dir="$(cd "$dir" 2>/dev/null && pwd)" || abs_dir="$CURRENT_DIR/$dir"
    else
        abs_dir="$dir"
    fi
    
    # Check if directory exists
    if [[ ! -d "$abs_dir" ]]; then
        error_log "Directory does not exist: $abs_dir"
        exit 1
    fi
    
    # Add to volumes if not already mounted
    if [[ ! " ${VOLUMES[@]} " =~ " -v $abs_dir:$abs_dir " ]]; then
        VOLUMES+=("-v" "$abs_dir:$abs_dir")
        debug_log "Mounting directory: $abs_dir"
    fi
done

# Mount tmp directory for temporary files
VOLUMES+=("-v" "/tmp:/tmp")

# Optional: Mount custom database directory if it exists
if [[ -d "/media/db" ]]; then
    VOLUMES+=("-v" "/media/db:/media/db")
    debug_log "Mounting custom database: /media/db"
fi

# Determine if we need interactive mode
TTY_FLAG=""
if [ -t 0 ]; then
    TTY_FLAG="-it"
fi

debug_log "Final volume mounts: ${VOLUMES[*]}"
debug_log "Final arguments: ${ARGS[*]}"

# Run Docker container
info_log "Running plsMD in Docker..."

docker run \
    --rm \
    $TTY_FLAG \
    "${VOLUMES[@]}" \
    --user "$(id -u):$(id -g)" \
    -w /data \
    plsmd "${ARGS[@]}"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    error_log "plsMD exited with error code $exit_code"
    echo ""
    echo "Troubleshooting tips:"
    echo "  1. Check that all input paths exist"
    echo "  2. Ensure you have read/write permissions"
    echo "  3. Run with DEBUG=1 for more information:"
    echo "     DEBUG=1 plsMD ${ARGS[*]}"
fi

exit $exit_code
