#!/usr/bin/env bash
# run.sh : lance tout le pipeline ClinVar (sein / ovaire)
# Exemple : ./run.sh -g BRCA1,BRCA2 -o sein
set -euo pipefail

# racine du projet (utile pour cron)
cd "$(dirname "$0")"

source scripts/lib.sh

# valeurs par défaut
GENES=""
ORGANE="les_deux"
CLASSIFICATION="Pathogenic"
RUN_TESTS=0

RAW_DIR="data/raw"
CLEAN_DIR="data/clean"
RESULTS_DIR="results"
RAW_FILE="$RAW_DIR/variant_summary.txt.gz"
CLEAN_FILE="$CLEAN_DIR/variants_filtres.tsv"
STATS_FILE="$RESULTS_DIR/stats.csv"
COUNTS_FILE="$RESULTS_DIR/comptages_shell.txt"
FIG_DIR="$RESULTS_DIR/figures"
RESUME_FILE="$RESULTS_DIR/resume.json"

if [[ -x ".venv/bin/python" ]]; then
    PYTHON=".venv/bin/python"
else
    PYTHON="python3"
fi

usage() {
    cat <<EOF
Utilisation : ./run.sh [-h] [-g GENES] [-o ORGANE] [-c CLASSIF] [-v] [-t]

Pipeline ClinVar : cancer héréditaire du sein et de l'ovaire.
Télécharge les données, les filtre, calcule les statistiques, produit
les figures, un résumé JSON et une archive datée des résultats.

Options :
  -h           affiche cette aide et quitte
  -g GENES     gènes à garder, séparés par des virgules (ex. BRCA1,BRCA2)
               par défaut : tous les gènes
  -o ORGANE    sein, ovaire ou les_deux (par défaut : les_deux)
  -c CLASSIF   classification étudiée (par défaut : Pathogenic)
  -v           mode verbeux (affiche les messages DEBUG)
  -t           lance aussi les tests automatiques (pytest tests/)

Codes de retour :
  0 succès | 1 erreur générale | 2 mauvaise option
  3 fichier manquant | 4 somme de contrôle invalide

Exemple :
  ./run.sh -g BRCA1,BRCA2 -o sein
EOF
}

# run_step "nom" commande... : lance une étape et s'arrête si elle échoue
run_step() {
    local nom="$1"
    shift
    log "=== Étape : $nom ==="
    local debut=$SECONDS
    local code=0
    "$@" || code=$?
    if (( code != 0 )); then
        die "L'étape « $nom » a échoué" "$code"
    fi
    log "Étape « $nom » terminée en $(( SECONDS - debut )) s"
}

# select_organe fichier organe : garde l'organe demandé + les_deux
select_organe() {
    local fichier="$1"
    local organe="$2"
    local tmp="$TMP_DIR/organe.tsv"
    awk -F'\t' -v o="$organe" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == "organe") c = i; print; next }
        $c == o || $c == "les_deux"
    ' "$fichier" > "$tmp"
    log "Organe $organe : $(( $(wc -l < "$tmp") - 1 )) variants gardés"
    mv "$tmp" "$fichier"
}

while getopts ":hg:o:c:vt" opt; do
    case "$opt" in
        h) usage; exit "$E_OK" ;;
        g) GENES="$OPTARG" ;;
        o) ORGANE="$OPTARG" ;;
        c) CLASSIFICATION="$OPTARG" ;;
        v) VERBOSE=1 ;;
        t) RUN_TESTS=1 ;;
        :) echo "Erreur : l'option -$OPTARG attend une valeur." >&2
           usage >&2; exit "$E_OPTION" ;;
        \?) echo "Erreur : option inconnue -$OPTARG" >&2
           usage >&2; exit "$E_OPTION" ;;
    esac
done
shift $((OPTIND - 1))
export VERBOSE

case "$ORGANE" in
    sein|ovaire|les_deux) ;;
    *) echo "Erreur : organe « $ORGANE » invalide (sein, ovaire ou les_deux)." >&2
       usage >&2; exit "$E_OPTION" ;;
esac
if [[ -n "$GENES" && ! "$GENES" =~ ^[A-Za-z0-9_-]+(,[A-Za-z0-9_-]+)*$ ]]; then
    echo "Erreur : liste de gènes « $GENES » invalide (ex. BRCA1,BRCA2)." >&2
    exit "$E_OPTION"
fi

trap cleanup EXIT
trap 'die "Interruption par l utilisateur (Ctrl+C)" 130' INT
trap 'die "Arrêt demandé (signal TERM)" 143' TERM
trap 'die "Erreur ligne $LINENO : la commande « $BASH_COMMAND » a échoué" "$E_GENERAL"' ERR

mkdir -p "$RAW_DIR" "$CLEAN_DIR" "$RESULTS_DIR" "$FIG_DIR" \
         "$RESULTS_DIR/archives" logs
make_tmpdir

log "Démarrage du pipeline (gènes=${GENES:-tous}, organe=$ORGANE, classification=$CLASSIFICATION)"
log "Journal : $LOG_FILE"

for s in scripts/download.sh scripts/check_data.py scripts/filter.sh \
         scripts/quick_counts.sh scripts/analyse.py scripts/figures.py scripts/resume.py scripts/archive.sh; do
    check_file "$s"
done

run_step "1. Téléchargement" \
    bash scripts/download.sh
check_file "$RAW_FILE"

run_step "1b. Contrôle des colonnes" \
    "$PYTHON" scripts/check_data.py --input "$RAW_FILE"

run_step "2. Filtrage" \
    bash scripts/filter.sh -i "$RAW_FILE" -o "$CLEAN_FILE" -g "$GENES"
check_file "$CLEAN_FILE"

if [[ "$ORGANE" != "les_deux" ]]; then
    run_step "2b. Sélection de l'organe : $ORGANE" \
        select_organe "$CLEAN_FILE" "$ORGANE"
    check_file "$CLEAN_FILE"
fi

run_step "3a. Comptages shell" \
    bash scripts/quick_counts.sh "$CLEAN_FILE" "$COUNTS_FILE"

run_step "3b. Analyse statistique" \
    "$PYTHON" scripts/analyse.py --input "$CLEAN_FILE" --output "$STATS_FILE" \
        --classification "$CLASSIFICATION" --comptages-shell "$COUNTS_FILE"
check_file "$STATS_FILE"

run_step "4. Figures" \
    "$PYTHON" scripts/figures.py --input "$STATS_FILE" --outdir "$FIG_DIR"

run_step "5. Résumé JSON" \
    "$PYTHON" scripts/resume.py --stats "$STATS_FILE" --log "$LOG_FILE" \
        --output "$RESUME_FILE" --genes "$GENES" --organe "$ORGANE" \
        --classification "$CLASSIFICATION"
check_file "$RESUME_FILE"

run_step "6. Archivage" \
    bash scripts/archive.sh "$RESULTS_DIR"

if (( RUN_TESTS == 1 )); then
    run_step "Tests automatiques" "$PYTHON" -m pytest -q tests/
fi

log "Pipeline terminé avec succès. Résumé : $RESUME_FILE"
exit "$E_OK"
