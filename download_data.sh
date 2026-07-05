#!/bin/bash
# Exit on error
set -e

# Default settings (User can override these when running the script or editing this file)
USER="your-github-username"
REPO="your-repository-name"
TAG="v1.0.0-data"

show_help() {
    echo "Usage: bash download_data.sh [github_username] [repository_name] [release_tag]"
    echo "Example: bash download_data.sh Subaru-PFS netflow-pipeline v1.0.0-data"
}

# Allow override via arguments
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_help
    exit 0
fi

if [ -n "$1" ]; then
    USER="$1"
fi
if [ -n "$2" ]; then
    REPO="$2"
fi
if [ -n "$3" ]; then
    TAG="$3"
fi

echo "============================================================="
echo "PFS Netflow Dataset Downloader"
echo "GitHub Repository: ${USER}/${REPO}"
echo "Release Tag:       ${TAG}"
echo "============================================================="

# cosmos
if [ ! -d "cosmos" ]; then
    echo "Downloading cosmos dataset..."
    curl -L -o cosmos.tar.gz "https://github.com/${USER}/${REPO}/releases/download/${TAG}/cosmos.tar.gz"
    echo "Extracting cosmos dataset..."
    tar -zxvf cosmos.tar.gz
    rm cosmos.tar.gz
else
    echo "Directory 'cosmos' already exists. Skipping."
fi

# xmm_lss
if [ ! -d "xmm_lss" ]; then
    echo "Downloading xmm_lss dataset..."
    curl -L -o xmm_lss.tar.gz "https://github.com/${USER}/${REPO}/releases/download/${TAG}/xmm_lss.tar.gz"
    echo "Extracting xmm_lss dataset..."
    tar -zxvf xmm_lss.tar.gz
    rm xmm_lss.tar.gz
else
    echo "Directory 'xmm_lss' already exists. Skipping."
fi

echo "Dataset download and extraction completed successfully!"
