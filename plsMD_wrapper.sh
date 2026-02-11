#!/bin/bash
# plsMD Docker Wrapper

set -e

CURRENT_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' 

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

if ! docker info >/dev/null 2>&1; then
    error_log "Docker is not running or not installed"
    exit 1
fi


if ! docker image inspect plsmd >/dev/null 2>&1; then
    error_log "plsmd Docker image not found. Please build it first:"
    echo "  docker build -t plsmd ." >&2
    exit 1
fi

debug_log "Script directory: $SCRIPT_DIR"
debug_log "Current directory: $CURRENT_DIR"
debug_log "Arguments: $*"


declare -A DIRS_TO_MOUNT

DIRS_TO_MOUNT["$CURRENT_DIR"]=1
debug_log "Will mount current directory: $CURRENT_DIR"

add_dir_to_mount() {
    local path="$1"
    

    [ -z "$path" ] && return
    

    if [[ ! "$path" =~ ^/ ]]; then
        path="$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")" || path="$CURRENT_DIR/$path"
    fi
    

    if [ -f "$path" ]; then
        local dir=$(dirname "$path")
        DIRS_TO_MOUNT["$dir"]=1
        debug_log "Will mount directory (from file): $dir"
    elif [ -d "$path" ]; then
        DIRS_TO_MOUNT["$path"]=1
        debug_log "Will mount directory: $path"
    else

        local parent=$(dirname "$path")
        if [ -d "$parent" ]; then
            DIRS_TO_MOUNT["$parent"]=1
            debug_log "Will mount parent directory: $parent"
        fi
    fi
}

prev_arg=""
for arg in "$@"; do

    if [[ "$prev_arg" =~ ^--(dir|input|output|db|path|IS_db)$ ]]; then
        add_dir_to_mount "$arg"
    fi
    
    if [[ "$arg" =~ ^--(dir|input|output|db|path|IS_db)=(.+)$ ]]; then
        add_dir_to_mount "${BASH_REMATCH[2]}"
    fi

    if [[ "$arg" =~ ^/[^-] ]] && [ -e "$arg" ]; then
        add_dir_to_mount "$arg"
    fi
    
    prev_arg="$arg"
done

VOLUME_ARGS=()

VOLUME_ARGS+=("-v" "$CURRENT_DIR:/data")

for dir in "${!DIRS_TO_MOUNT[@]}"; do

    if [ "$dir" != "$CURRENT_DIR" ]; then
        VOLUME_ARGS+=("-v" "$dir:$dir")
    fi
done

VOLUME_ARGS+=("-v" "/tmp:/tmp")

debug_log "Volume mounts: ${VOLUME_ARGS[*]}"


TTY_FLAG=""
if [ -t 0 ]; then
    TTY_FLAG="-it"
fi

info_log "Running plsMD in Docker..."

debug_log "Docker command: docker run --rm $TTY_FLAG ${VOLUME_ARGS[*]} --user $(id -u):$(id -g) -w /data plsmd $*"

docker run \
    --rm \
    $TTY_FLAG \
    "${VOLUME_ARGS[@]}" \
    --user "$(id -u):$(id -g)" \
    -w /data \
    plsmd "$@"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    error_log "plsMD exited with error code $exit_code"
    echo ""
    echo "Troubleshooting tips:"
    echo "  1. Check that all input paths exist"
    echo "  2. Ensure you have read/write permissions"
    echo "  3. Run with DEBUG=1 for more information:"
    echo "     DEBUG=1 plsMD $*"
    echo "  4. Verify paths are accessible:"
    for arg in "$@"; do
        if [[ "$arg" =~ ^/ ]] && [ -e "$arg" ]; then
            echo "     ls -la $arg"
        fi
    done
fi

exit $exit_code
