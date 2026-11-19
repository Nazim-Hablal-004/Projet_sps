#!/usr/bin/env bash
#
# quick_counts.sh - Comptages rapides en shell (Personne 4)
#
# Résume le fichier filtré avec des tubes, avant l'analyse Python. Les nombres
# servent de contrôle : analyse.py doit retrouver exactement les mêmes.
#
# Usage : scripts/quick_counts.sh [entree.tsv] [sortie.txt]
#   entree.tsv  fichier filtré (défaut : data/clean/variants_filtres.tsv)
#   sortie.txt  fichier de résultats (défaut : results/comptages_shell.txt)
#
# Codes de retour : 0 = succès, 1 = colonne absente, 3 = fichier manquant ou vide.

set -euo pipefail

# Le dossier du script permet de charger lib.sh quel que soit le dossier courant.
DOSSIER_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib.sh
source "$DOSSIER_SCRIPT/lib.sh"

ENTREE="${1:-data/clean/variants_filtres.tsv}"
SORTIE="${2:-results/comptages_shell.txt}"
NB_GENES=10

check_file "$ENTREE"

# Donne le numéro d'une colonne à partir de son nom dans l'en-tête.
# On ne code pas les numéros en dur : l'ordre des colonnes peut changer.
numero_colonne() {
    local nom="$1"
    local numero
    numero="$(awk -F'\t' -v nom="$nom" 'NR == 1 {
        for (i = 1; i <= NF; i++) if ($i == nom) print i
        exit
    }' "$ENTREE")"
    [[ -n "$numero" ]] || die "colonne absente de $ENTREE : $nom" 1
    echo "$numero"
}

COL_TYPE="$(numero_colonne type)"
COL_GENE="$(numero_colonne gene)"
COL_CLASSIF="$(numero_colonne classification)"
COL_ORGANE="$(numero_colonne organe)"

# Saute l'en-tête et retire les variants sans classification.
# Même règle que retirer_non_classes() dans analyse.py, sinon les nombres diffèrent.
lignes_classees() {
    tail -n +2 "$ENTREE" | awk -F'\t' -v c="$COL_CLASSIF" '
        { valeur = tolower($c) }
        valeur != "" && valeur != "-" && valeur != "not provided" && valeur !~ /^no classification/
    '
}

# Compte les valeurs d'une colonne, de la plus fréquente à la moins fréquente.
# uniq -c ne compte que les lignes identiques qui se suivent : il faut trier avant.
# LC_ALL=C donne le même ordre de tri sur toutes les machines.
compter() {
    local colonne="$1"
    lignes_classees | cut -f"$colonne" | LC_ALL=C sort | uniq -c | LC_ALL=C sort -nr
}

mkdir -p "$(dirname "$SORTIE")"
TEMPORAIRE="$SORTIE.tmp"

# awk 'NR <= N' plutôt que head -n N : head se ferme dès qu'il a ses lignes, sort
# reçoit alors SIGPIPE (code 141) et, avec pipefail, tout le script s'arrêterait.
{
    echo "## type"
    compter "$COL_TYPE"
    echo "## gene"
    compter "$COL_GENE" | awk -v n="$NB_GENES" 'NR <= n'
    echo "## classification"
    compter "$COL_CLASSIF"
    echo "## organe"
    compter "$COL_ORGANE"
} > "$TEMPORAIRE"

# Le fichier final n'apparaît qu'une fois complet.
mv "$TEMPORAIRE" "$SORTIE"

log "Comptages shell : $(lignes_classees | wc -l) variants classés, écrits dans $SORTIE"
