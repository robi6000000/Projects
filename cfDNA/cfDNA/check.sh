#!/bin/bash
# check_files.sh

MANIFEST=./data/manifest/Cristiano_manifest.csv
FRAG_DIR=./data/source/cristiano
WPS_DIR=./data/source/cristiano_features/cristiano_wps_download
HG19=./data/hg19/hg19.fa

echo "Checking file availability..."
echo ""

# Count manifest samples
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l)
echo "Total samples in manifest: $TOTAL"

# Check fragment files
FRAG_COUNT=$(ls -1 $FRAG_DIR/*.bgz 2>/dev/null | wc -l)
echo "Fragment files downloaded: $FRAG_COUNT / $TOTAL"

# Check WPS files
WPS_COUNT=$(ls -1 $WPS_DIR/*.bw 2>/dev/null | wc -l)
echo "WPS files downloaded: $WPS_COUNT / $TOTAL"

# Check hg19
if [ -f "$HG19" ]; then
    echo "✓ hg19.fa exists"
else
    echo "✗ hg19.fa MISSING!"
fi

echo ""
if [ $FRAG_COUNT -eq $TOTAL ] && [ $WPS_COUNT -eq $TOTAL ]; then
    echo "✓ All files ready!"
else
    echo "⚠ Missing files - download incomplete"
fi