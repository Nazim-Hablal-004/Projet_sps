#!/usr/bin/env bash
# archive.sh : crée une archive datée du dossier results/
# Usage : ./scripts/archive.sh [dossier_results]
set -euo pipefail

RESULTS_DIR="${1:-results}"
ARCHIVE_DIR="$RESULTS_DIR/archives"

source "$(dirname "$0")/lib.sh"

[ -d "$RESULTS_DIR" ] || die "dossier introuvable : $RESULTS_DIR" 3

mkdir -p "$ARCHIVE_DIR"
NOM="resultats_$(date '+%Y-%m-%d_%H%M%S').tar.gz"
ARCHIVE="$ARCHIVE_DIR/$NOM"

log "Création de l'archive $ARCHIVE"
# avec le journal du jour
JOURNAL=()
[ -f "$LOG_FILE" ] && JOURNAL=("$(realpath --relative-to=. "$LOG_FILE")")
tar -czf "$ARCHIVE" --exclude="$ARCHIVE_DIR" "$RESULTS_DIR" "${JOURNAL[@]}"

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
