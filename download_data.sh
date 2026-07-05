#!/bin/bash
# Exit on error
set -e

# Settings
USER="yasudank"
REPO="NetflowByUser"
TAG="v1.0.0-data"

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
