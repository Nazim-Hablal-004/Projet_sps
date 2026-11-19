#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/lib.sh"

URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MD5_URL="${URL}.md5"

RAW_DIR="data/raw"
FILE="$RAW_DIR/variant_summary.txt.gz"
# copie locale de la somme du NCBI : permet de revérifier le fichier sans réseau
MD5_FILE="$FILE.md5"
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
    log "Fichier : $FILE"
    log "Taille : $(stat -c %s "$FILE") octets"
    log "MD5 : $EXPECTED_MD5"
}

# vrai si le fichier local correspond à la somme enregistrée lors du téléchargement
fichier_local_intact() {
    [[ -f "$FILE" && -f "$MD5_FILE" ]] \
        && (cd "$RAW_DIR" && md5sum -c --status "$(basename "$MD5_FILE")")
}

# --retry : une coupure réseau passagère ne doit pas arrêter le pipeline
EXPECTED_MD5="$(curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors "$MD5_URL" \
    | awk '{print $1}')" || EXPECTED_MD5=""

if [[ ! "$EXPECTED_MD5" =~ ^[0-9a-f]{32}$ ]]; then
    # pas de réseau : on garde le fichier déjà téléchargé s'il est intact
    if fichier_local_intact; then
        EXPECTED_MD5="$(awk '{print $1}' "$MD5_FILE")"
        log "Réseau indisponible : fichier local gardé (somme de contrôle locale OK)" "ATTENTION"
        show_info
        exit 0
    fi
    die "Impossible de récupérer la somme de contrôle (problème réseau ?)" 1
fi
log "MD5 attendu : $EXPECTED_MD5"

if [[ -f "$FILE" ]]; then
    ACTUAL_MD5="$(md5sum "$FILE" | awk '{print $1}')"
    if [[ "$ACTUAL_MD5" == "$EXPECTED_MD5" ]]; then
        chmod 444 "$FILE"
        echo "$EXPECTED_MD5  variant_summary.txt.gz" > "$MD5_FILE"
        log "Déjà téléchargé, somme de contrôle OK"
        show_info
        exit 0
    fi
fi

TMP_DIR="$(mktemp -d -p data)"
log "Téléchargement de $URL"
curl -fL --retry 20 --retry-delay 5 --retry-all-errors --speed-time 30 --speed-limit 1000 -C - -o "$TMP_DIR/variant_summary.txt.gz" "$URL"

ACTUAL_MD5="$(md5sum "$TMP_DIR/variant_summary.txt.gz" | awk '{print $1}')"
if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    die "Somme de contrôle invalide (attendu $EXPECTED_MD5, obtenu $ACTUAL_MD5)" 4
fi

mv -f "$TMP_DIR/variant_summary.txt.gz" "$FILE"
chmod 444 "$FILE"
echo "$EXPECTED_MD5  variant_summary.txt.gz" > "$MD5_FILE"
log "Téléchargement OK, somme de contrôle OK"
show_info
