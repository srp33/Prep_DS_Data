#!/usr/bin/env bash

set -euo pipefail

IMAGE_NAME="bioconductor-3-23"
CONTAINER_NAME="bioconductor-3-23-container"

# Stop and remove the container if it is already running or exists.
if docker ps -a --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"; then
  echo "Stopping existing container: $CONTAINER_NAME"
  docker stop "$CONTAINER_NAME" >/dev/null 2>&1 || true

  echo "Removing existing container: $CONTAINER_NAME"
  docker rm "$CONTAINER_NAME" >/dev/null 2>&1 || true
fi

# Build the Docker image from the Dockerfile in the current directory.
# Unchanged layers are reused from Docker's build cache.
echo "Building image: $IMAGE_NAME"
docker build -t "$IMAGE_NAME" .

# Remove abandoned containers/networks and dangling images.
# Avoids -a and builder prune so the tagged image and build cache stay intact.
echo "Pruning abandoned Docker resources..."
docker system prune -f

# Run an R script inside the container as the current host user.
# Extra arguments after the script name are forwarded to Rscript.
run_rscript() {
  local script_file="$1"
  shift

  if [[ -z "$script_file" ]]; then
    echo "Usage: run_rscript <script.R> [args...]" >&2
    return 1
  fi

  if [[ ! -f "$script_file" ]]; then
    echo "Error: script file not found: $script_file" >&2
    return 1
  fi

  echo "Running script in container: $script_file $*"
  docker run -it --rm \
    --name "$CONTAINER_NAME" \
    --user "$(id -u):$(id -g)" \
    -v "$PWD":/project \
    -w /project \
    "$IMAGE_NAME" \
    Rscript "$script_file" "$@"
}

# run_rscript "download_CEL_files.R"
# run_rscript "download_metadata.R"


#https://github.com/Pungetello/DownSyndromeCuration/blob/main/MetadataAttributes.R