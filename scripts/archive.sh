#!/usr/bin/env bash
# archive.sh : crée une archive datée du dossier results/
# Usage : ./scripts/archive.sh [dossier_results]
set -euo pipefail

RESULTS_DIR="${1:-results}"
ARCHIVE_DIR="$RESULTS_DIR/archives"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { log "ERREUR : $*" >&2; exit 1; }

[ -d "$RESULTS_DIR" ] || die "dossier introuvable : $RESULTS_DIR"

mkdir -p "$ARCHIVE_DIR"
NOM="resultats_$(date '+%Y-%m-%d_%H%M%S').tar.gz"
ARCHIVE="$ARCHIVE_DIR/$NOM"

log "Création de l'archive $ARCHIVE"
tar -czf "$ARCHIVE" --exclude="$ARCHIVE_DIR" "$RESULTS_DIR"

# sha256sum sur Linux, shasum -a 256 sur macOS
if command -v sha256sum >/dev/null 2>&1; then
    HASH_CMD="sha256sum"
else
    HASH_CMD="shasum -a 256"
fi
( cd "$ARCHIVE_DIR" && $HASH_CMD "$NOM" > "$NOM.sha256" )

# vérifie que l'archive est lisible
tar -tzf "$ARCHIVE" > /dev/null || die "archive illisible : $ARCHIVE"

log "Archive OK : $ARCHIVE"
