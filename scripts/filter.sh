#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib.sh"

input="data/raw/variant_summary.txt.gz"
output="data/clean/variants_filtres.tsv"
genes=""
origin=""   # "" = toutes origines ; "germline" = germline seulement

usage() {
  cat <<EOF
Usage: $0 [-i input.gz] [-o output.tsv] [-g BRCA1,BRCA2] [-r germline|all]

Filtre ClinVar (corrections P4 du 6/11) :
  - garde Assembly = GRCh38
  - découpe PhenotypeList sur | et ; ; garde si UNE maladie contient
    un mot d organe (breast/ovarian/ovary) ET un mot de cancer
  - retire les classifications -, not provided, no classification...
  - garde TOUTES les autres classifications (Benign, VUS, etc.)
  - conserve VariationID pour ne pas confondre deux mutations
  - sortie : variation_id gene type classification organe

Options :
  -i  entrée (.gz ou .tsv)    défaut : data/raw/variant_summary.txt.gz
  -o  sortie                   défaut : data/clean/variants_filtres.tsv
  -g  gènes (BRCA1,BRCA2)      défaut : aucun filtre
  -r  germline  ou  all        défaut : toutes origines
  -h  aide
EOF
}

while getopts "hi:o:g:r:" opt; do
  case "$opt" in
    h) usage; exit 0 ;;
    i) input="$OPTARG" ;;
    o) output="$OPTARG" ;;
    g) genes="$OPTARG" ;;
    r) origin="$OPTARG" ;;
    \?) usage >&2; exit 2 ;;
    :)  usage >&2; exit 2 ;;
  esac
done

# Valider l'option -r
if [[ -n "$origin" && "$origin" != "germline" && "$origin" != "all" ]]; then
  die "Option -r invalide : '$origin' (attendu : germline ou all)" 2
fi
[[ "$origin" == "all" ]] && origin=""

check_file "$input"
mkdir -p "$(dirname "$output")" logs

if [[ "$input" == *.gz ]]; then
  reader=zcat
else
  reader=cat
fi

before=$($reader "$input" | wc -l)
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT INT TERM

# En-tête écrit à part : sort -u ne doit pas le trier
printf 'variation_id\tgene\ttype\tclassification\torgane\n' > "$tmp"

$reader "$input" | awk -F'\t' -v OFS='\t' -v genes="$genes" -v origin="$origin" '
# --- Fonctions de matching (correction 3) ---
function is_cancer(d,   x) {
  x = tolower(d)
  return (x ~ /cancer|carcinoma|neoplasm|tumor|tumour|malignan/)
}
function has_breast(d,   x) {
  x = tolower(d)
  return (x ~ /breast/)
}
function has_ovarian(d,   x) {
  x = tolower(d)
  return (x ~ /ovarian|ovary/)
}

BEGIN {
  gene_filter = (genes != "")
  n = split(genes, g, ",")
  for (i = 1; i <= n; i++) {
    gsub(/^[ \t]+|[ \t]+$/, "", g[i])
    if (g[i] != "") wanted[toupper(g[i])] = 1
  }
}

NR == 1 {
  for (i = 1; i <= NF; i++) idx[$i] = i
  # VariationID est obligatoire (correction 1)
  split("VariationID GeneSymbol Type ClinicalSignificance PhenotypeList Assembly", req, " ")
  for (i = 1; i <= 6; i++) {
    if (!(req[i] in idx)) {
      print "Colonne manquante : " req[i] > "/dev/stderr"
      exit 1
    }
  }
  next
}

{
  # Filtre assembly
  if ($(idx["Assembly"]) != "GRCh38") next

  # Filtre origine (option -r)
  if (origin != "") {
    if (!("OriginSimple" in idx)) {
      print "Option -r utilisee mais colonne OriginSimple absente" > "/dev/stderr"
      exit 1
    }
    if ($(idx["OriginSimple"]) != origin) next
  }

  # Correction 3 : split et matching par maladie
  pheno = $(idx["PhenotypeList"])
  ndis = split(pheno, dis, /[|;]/)
  hb = 0; ho = 0
  for (i = 1; i <= ndis; i++) {
    if (!is_cancer(dis[i])) continue
    if (has_breast(dis[i]))  hb = 1
    if (has_ovarian(dis[i])) ho = 1
  }
  if (!hb && !ho) next

  # Filtre gènes
  gene = $(idx["GeneSymbol"])
  if (gene_filter && !(toupper(gene) in wanted)) next

  # Correction 2 : nettoyage classification
  cls = tolower($(idx["ClinicalSignificance"]))
  gsub(/^[ \t]+|[ \t]+$/, "", cls)
  if (cls == "" \
      || cls == "-" \
      || cls == "not provided" \
      || cls == "not_provided" \
      || cls ~ /^no classification/) next

  # Organe
  if (hb && ho) organe = "les_deux"
  else if (hb)  organe = "sein"
  else          organe = "ovaire"

  # Correction 1 : VariationID en sortie
  vid = $(idx["VariationID"])
  print vid, gene, $(idx["Type"]), cls, organe
}' | cut -f1-5 | LC_ALL=C sort -u >> "$tmp"

after=$(($(wc -l < "$tmp") - 1))
mv "$tmp" "$output"
trap - EXIT INT TERM

log "Filtrage termine : $before lignes lues, $after variants conserves dans $output"
