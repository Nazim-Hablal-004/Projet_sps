#!/usr/bin/env bash
# =============================================================================
# lib.sh — fonctions communes du pipeline (Personne 1)
#
# Utilisation dans n'importe quel script :
#     source "$(dirname "$0")/lib.sh"        # depuis scripts/
#     source "scripts/lib.sh"                # depuis la racine
#
# Fonctions fournies :
#     log "message" [NIVEAU]   -> ligne horodatée à l'écran + dans le journal
#     die "message" [code]     -> message d'erreur puis sortie avec le code
#     check_file chemin        -> vérifie qu'un fichier existe et n'est pas vide
#     make_tmpdir              -> crée le dossier temporaire (mktemp -d)
#     cleanup                  -> supprime le dossier temporaire (pour trap)
# =============================================================================

# Évite de charger la bibliothèque deux fois
if [[ -n "${LIB_SH_LOADED:-}" ]]; then
    return 0
fi
LIB_SH_LOADED=1

# --- Codes de retour (documentés dans le README) -----------------------------
readonly E_OK=0          # succès
readonly E_GENERAL=1     # erreur générale
readonly E_OPTION=2      # mauvaise option
readonly E_MISSING=3     # fichier manquant ou vide
readonly E_CHECKSUM=4    # somme de contrôle invalide

# --- Chemins -----------------------------------------------------------------
# Racine du projet = dossier parent de scripts/
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs}"
mkdir -p "$LOG_DIR"
# Un journal par jour. "export" pour que les scripts lancés par run.sh
# écrivent dans le même fichier.
export LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_$(date '+%Y-%m-%d').log}"

# Mode verbeux (activé par run.sh -v)
VERBOSE="${VERBOSE:-0}"

# -----------------------------------------------------------------------------
# log "message" [NIVEAU]
# Écrit "2026-11-04 14:03:12 [INFO] message" à l'écran et dans le journal.
# Le niveau DEBUG n'est affiché qu'en mode verbeux.
# -----------------------------------------------------------------------------
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
        echo "$ligne" >&2          # les erreurs vont sur la sortie d'erreur
    else
        echo "$ligne"
    fi
}

# -----------------------------------------------------------------------------
# die "message" [code]
# Écrit l'erreur (journal + sortie d'erreur) puis quitte avec le code donné.
# -----------------------------------------------------------------------------
die() {
    local message="$1"
    local code="${2:-$E_GENERAL}"
    log "$message (code $code)" "ERREUR"
    exit "$code"
}

# -----------------------------------------------------------------------------
# check_file chemin
# Arrête le script avec le code 3 si le fichier n'existe pas ou est vide.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# make_tmpdir
# Crée un dossier temporaire au nom unique et l'exporte dans TMP_DIR,
# pour que les autres scripts puissent y écrire leurs fichiers intermédiaires.
# -----------------------------------------------------------------------------
make_tmpdir() {
    TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/pipeline_XXXXXX")"
    export TMP_DIR
    log "Dossier temporaire : $TMP_DIR" "DEBUG"
}

# -----------------------------------------------------------------------------
# cleanup
# Supprime le dossier temporaire. Appelée automatiquement par :
#     trap cleanup EXIT
# -----------------------------------------------------------------------------
cleanup() {
    if [[ -n "${TMP_DIR:-}" && -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
        log "Nettoyage : dossier temporaire supprimé" "DEBUG"
    fi
}
