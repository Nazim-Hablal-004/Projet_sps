#!/usr/bin/env bash
set -euo pipefail

URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MD5_URL="${URL}.md5"

RAW_DIR="data/raw"
FILE="$RAW_DIR/variant_summary.txt.gz"
mkdir -p "$RAW_DIR"

TMP_DIR=""
cleanup() {
    if [[ -n "$TMP_DIR" && -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

show_info() {
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "File: $FILE"
    echo "Size: $(stat -c %s "$FILE") bytes"
    echo "MD5:  $EXPECTED_MD5"
}

EXPECTED_MD5="$(curl -fsSL "$MD5_URL" | awk '{print $1}')" || {
    echo "ERROR: cannot fetch the checksum (network problem?)" >&2
    exit 1
}
echo "Expected MD5: $EXPECTED_MD5"

if [[ -f "$FILE" ]]; then
    ACTUAL_MD5="$(md5sum "$FILE" | awk '{print $1}')"
    if [[ "$ACTUAL_MD5" == "$EXPECTED_MD5" ]]; then
        chmod 444 "$FILE"
        echo "Already downloaded, checksum OK"
        show_info
        exit 0
    fi
fi

TMP_DIR="$(mktemp -d -p data)"
curl -fL --retry 20 --retry-delay 5 --retry-all-errors --speed-time 30 --speed-limit 1000 -C - -o "$TMP_DIR/variant_summary.txt.gz" "$URL"

ACTUAL_MD5="$(md5sum "$TMP_DIR/variant_summary.txt.gz" | awk '{print $1}')"
if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "ERROR: invalid checksum" >&2
    exit 4
fi

mv -f "$TMP_DIR/variant_summary.txt.gz" "$FILE"
chmod 444 "$FILE"
echo "Download OK, checksum OK"s
show_info
