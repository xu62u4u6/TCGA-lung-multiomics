#!/bin/bash
set -e

GDC_CLIENT=/home/blue0228/.local/bin/gdc-client
TOKEN=/home/blue0228/tcga-token.txt
MANIFEST_DIR=data/manifests
RAW_DIR=data/raw

for dtype in expression mutation cnv cnv_gene methylation clinical; do
    echo "=========================================="
    echo "  Downloading: $dtype"
    echo "=========================================="
    mkdir -p ${RAW_DIR}/${dtype}
    $GDC_CLIENT download \
        -m ${MANIFEST_DIR}/manifest_${dtype}.tsv \
        -t ${TOKEN} \
        -d ${RAW_DIR}/${dtype} \
        -n 8 \
        --retry-amount 5
    echo "  Done: $dtype"
    echo ""
done

echo "All downloads complete!"
du -sh ${RAW_DIR}/*/
