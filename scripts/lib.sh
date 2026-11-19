#!/usr/bin/env bash
# lib.sh : fonctions communes du pipeline (log, die, check_file, make_tmpdir, cleanup)
# A charger avec : source "$(dirname "$0")/lib.sh"

if [[ -n "${LIB_SH_LOADED:-}" ]]; then
    return 0
fi
LIB_SH_LOADED=1

# Codes de retour
readonly E_OK=0
readonly E_GENERAL=1
readonly E_OPTION=2
readonly E_MISSING=3
readonly E_CHECKSUM=4

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs}"
mkdir -p "$LOG_DIR"
# exporté pour que tous les scripts écrivent dans le même journal
export LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_$(date '+%Y-%m-%d').log}"
VERBOSE="${VERBOSE:-0}"

# log "message" [NIVEAU]
log() {
    local message="$1"
    local niveau="${2:-INFO}"
    if [[ "$niveau" == "DEBUG" && "$VERBOSE" != "1" ]]; then
        return 0
    fi
    local ligne
    ligne="$(date '+%Y-%m-%d %H:%M:%S') [$niveau] $message"
    echo "$ligne" >> "$LOG_FILE"
    if [[ "$niveau" == "ERREUR" ]]; then
        echo "$ligne" >&2
    else
        echo "$ligne"
    fi
}

# die "message" [code]
die() {
    local message="$1"
    local code="${2:-$E_GENERAL}"
    log "$message (code $code)" "ERREUR"
    exit "$code"
}

# check_file chemin : le fichier doit exister et ne pas être vide
check_file() {
    local chemin="$1"
    if [[ ! -e "$chemin" ]]; then
        die "Fichier introuvable : $chemin" "$E_MISSING"
    fi
    if [[ ! -s "$chemin" ]]; then
        die "Fichier vide : $chemin" "$E_MISSING"
    fi
    log "Fichier OK : $chemin" "DEBUG"
}

make_tmpdir() {
    TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/pipeline_XXXXXX")"
    export TMP_DIR
    log "Dossier temporaire : $TMP_DIR" "DEBUG"
}

# appelée par trap à la fin du script
cleanup() {
    if [[ -n "${TMP_DIR:-}" && -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
        log "Nettoyage : dossier temporaire supprimé" "DEBUG"
    fi
}
