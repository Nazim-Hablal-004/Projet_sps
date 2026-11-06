#!/usr/bin/env python3
"""Analyse statistique des variants du sein et de l'ovaire (Personne 4).

Répond à la question du projet : quels gènes et quels types de mutations sont
le plus souvent pathogènes, et le profil diffère-t-il entre les deux organes ?

Entrée  : le fichier filtré (colonnes gene, type, classification, organe).
Sorties : results/stats.csv (comptages et proportions)
          results/stats.json (tests du khi-deux et métadonnées)

Codes de retour : 0 = succès, 1 = autre erreur, 3 = fichier manquant ou vide.
"""

import argparse
import json
import math
import os
import re
import sys

import pandas as pd
from scipy.stats import chi2_contingency

CODE_ERREUR = 1
CODE_FICHIER_MANQUANT = 3

COLONNES_REQUISES = ["gene", "type", "classification", "organe"]
COLONNES_STATS = [
    "section",
    "categorie",
    "organe",
    "n_variants",
    "n_pathogenes",
    "proportion_pathogene",
    "part_parmi_pathogenes",
]

# Chaque variant appartient à un seul de ces trois groupes.
ORGANES = ["sein", "ovaire", "les_deux"]
# Vue complète d'un organe : les variants qui lui sont propres et ceux partagés.
# Ces groupes se recoupent, on ne les utilise donc pas pour le khi-deux.
GROUPES_ORGANES = {
    "sein": ["sein"],
    "ovaire": ["ovaire"],
    "les_deux": ["les_deux"],
    "sein_total": ["sein", "les_deux"],
    "ovaire_total": ["ovaire", "les_deux"],
}

AUTRES = "autres"
NB_GENES_STATS = 15
# Les tests regroupent les catégories rares : le khi-deux demande des effectifs
# attendus d'au moins 5 dans chaque case.
NB_TYPES_TEST = 3
NB_GENES_TEST = 8
EFFECTIF_ATTENDU_MIN = 5


def erreur(message, code=CODE_ERREUR):
    """Affiche un message sur la sortie d'erreur et arrête le programme."""
    print(f"ERREUR : {message}", file=sys.stderr)
    sys.exit(code)


def normaliser_organe(valeur):
    """Accepte « les deux » comme « les_deux », pour une seule écriture ensuite."""
    return valeur.strip().lower().replace(" ", "_")


def charger_variants(chemin):
    """Lit le fichier filtré et vérifie qu'il respecte le contrat entre les étapes."""
    if not os.path.isfile(chemin):
        erreur(f"fichier introuvable : {chemin}", CODE_FICHIER_MANQUANT)

    # dtype=str et keep_default_na=False : on garde le texte tel quel. Sinon
    # pandas transformerait « NA » ou un champ vide en valeur manquante, groupby
    # les ignorerait, et les nombres ne seraient plus ceux du shell.
    try:
        df = pd.read_csv(chemin, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        erreur(f"fichier vide : {chemin}", CODE_FICHIER_MANQUANT)

    manquantes = [c for c in COLONNES_REQUISES if c not in df.columns]
    if manquantes:
        erreur(f"colonnes absentes de {chemin} : {', '.join(manquantes)}")

    # Sans identifiant, le sort -u du filtrage fusionne des variants différents.
    if "variation_id" not in df.columns:
        print(
            "ATTENTION : pas de colonne variation_id, les comptages sont faux "
            "si le filtrage a dédoublonné sur les seules colonnes restantes.",
            file=sys.stderr,
        )

    df["organe"] = df["organe"].map(normaliser_organe)
    inconnus = sorted(set(df["organe"]) - set(ORGANES))
    if inconnus:
        erreur(f"valeurs inattendues dans la colonne organe : {', '.join(inconnus)}")
    return df


def retirer_non_classes(df):
    """Retire les variants sans classification (« - » dans ClinVar).

    Ils ne sont ni pathogènes ni bénins : les compter ferait baisser les
    proportions à tort. Même règle que lignes_classees() dans quick_counts.sh.
    """
    valeur = df["classification"].str.lower()
    non_classe = valeur.isin(["", "-", "not provided"]) | valeur.str.startswith(
        "no classification"
    )
    return df[~non_classe].copy()


def marquer_pathogenes(df, motif="Pathogenic"):
    r"""Ajoute la colonne est_pathogene (vrai ou faux).

    \b impose un mot entier. « Pathogenic » trouve donc « Likely pathogenic » et
    « Pathogenic/Likely pathogenic », mais pas « Conflicting classifications of
    pathogenicity », qui n'est pas une classification pathogène.
    """
    expression = rf"\b{re.escape(motif)}\b"
    resultat = df.copy()
    resultat["est_pathogene"] = df["classification"].str.contains(
        expression, case=False, regex=True
    )
    return resultat


def proportion_pathogene(df, colonne):
    """Calcule, pour chaque valeur d'une colonne, la proportion de pathogènes.

    est_pathogene vaut vrai (1) ou faux (0) : sa moyenne est donc une proportion.
    Tri par nombre de pathogènes décroissant, puis par nom pour un ordre stable.
    """
    table = (
        df.groupby(colonne)["est_pathogene"]
        .agg(n_variants="size", n_pathogenes="sum", proportion_pathogene="mean")
        .reset_index()
    )
    table["n_pathogenes"] = table["n_pathogenes"].astype(int)
    return table.sort_values(
        ["n_pathogenes", "n_variants", colonne],
        ascending=[False, False, True],
        ignore_index=True,
    )


def regrouper_rares(df, colonne, n_gardes):
    """Garde les n_gardes valeurs les plus pathogènes et remplace le reste par « autres »."""
    gardes = proportion_pathogene(df, colonne)[colonne].head(n_gardes)
    return df[colonne].where(df[colonne].isin(gardes), AUTRES)


def resumer_par_categorie(df, colonne, section, organe="tous"):
    """Construit les lignes de stats.csv pour une colonne, avec « autres » en dernier."""
    table = proportion_pathogene(df, colonne).rename(columns={colonne: "categorie"})
    table["section"] = section
    table["organe"] = organe

    # La part décrit le profil : sur 100 pathogènes de ce groupe, combien par catégorie.
    total = table["n_pathogenes"].sum()
    table["part_parmi_pathogenes"] = table["n_pathogenes"] / total if total else 0.0

    est_autres = table["categorie"] == AUTRES
    return pd.concat([table[~est_autres], table[est_autres]], ignore_index=True)


def profil_par_organe(df, colonne, section):
    """Répète le résumé dans chacun des trois groupes d'organes."""
    morceaux = [
        resumer_par_categorie(df[df["organe"] == organe], colonne, section, organe)
        for organe in ORGANES
        if (df["organe"] == organe).any()
    ]
    return pd.concat(morceaux, ignore_index=True)


def resumer_groupes(df):
    """Construit les sections global et organe de stats.csv."""
    total_pathogenes = df["est_pathogene"].sum()
    groupes = [("global", "tous", df)]
    for nom, organes in GROUPES_ORGANES.items():
        sous_ensemble = df[df["organe"].isin(organes)]
        if not sous_ensemble.empty:
            groupes.append(("organe", nom, sous_ensemble))

    lignes = []
    for section, organe, groupe in groupes:
        n_pathogenes = int(groupe["est_pathogene"].sum())
        lignes.append(
            {
                "section": section,
                "categorie": "tous",
                "organe": organe,
                "n_variants": len(groupe),
                "n_pathogenes": n_pathogenes,
                "proportion_pathogene": groupe["est_pathogene"].mean(),
                "part_parmi_pathogenes": (
                    n_pathogenes / total_pathogenes if total_pathogenes else 0.0
                ),
            }
        )
    return pd.DataFrame(lignes)


def construire_stats(df):
    """Assemble toutes les sections de stats.csv en un seul tableau."""
    avec_genes_regroupes = df.assign(gene=regrouper_rares(df, "gene", NB_GENES_STATS))
    sections = [
        resumer_groupes(df),
        resumer_par_categorie(avec_genes_regroupes, "gene", "gene"),
        resumer_par_categorie(df, "type", "type"),
        resumer_par_categorie(df, "classification", "classification"),
        profil_par_organe(avec_genes_regroupes, "gene", "gene_organe"),
        profil_par_organe(df, "type", "type_organe"),
    ]
    return pd.concat(sections, ignore_index=True)[COLONNES_STATS]


def khi_deux(table):
    """Applique le test du khi-deux à un tableau de contingence.

    Le test mesure si l'écart à la situation « aucun lien » peut venir du hasard
    (p-value). Avec de grands effectifs, un écart minime devient significatif :
    le V de Cramér (de 0 à 1) indique en plus si le lien est fort.
    """
    # Une ligne ou une colonne vide rendrait le test impossible.
    table = table.loc[table.sum(axis=1) > 0, table.sum(axis=0) > 0]
    if table.shape[0] < 2 or table.shape[1] < 2:
        return {
            "statut": "non calculable",
            "raison": "il faut au moins deux groupes et deux catégories non vides",
        }

    # correction=False : la formule classique, sans la correction de Yates
    # (appliquée sinon aux tableaux 2x2), inutile avec de tels effectifs.
    chi2, p_value, ddl, attendus = chi2_contingency(table, correction=False)
    effectif = int(table.to_numpy().sum())
    v_cramer = math.sqrt(chi2 / (effectif * (min(table.shape) - 1)))
    return {
        "statut": "calculé",
        "chi2": round(float(chi2), 3),
        "ddl": int(ddl),
        "p_value": float(p_value),
        "v_cramer": round(v_cramer, 4),
        "effectif": effectif,
        "effectif_attendu_min": round(float(attendus.min()), 2),
        "fiable": bool(attendus.min() >= EFFECTIF_ATTENDU_MIN),
    }


def comparer_organes(df):
    """Compare le sein seul et l'ovaire seul avec trois tests du khi-deux.

    Le groupe les_deux est exclu : le test exige que chaque variant ne compte
    que dans un seul groupe.
    """
    exclusifs = df[df["organe"].isin(["sein", "ovaire"])]
    pathogenes = exclusifs[exclusifs["est_pathogene"]]
    return {
        "proportion_pathogene": khi_deux(
            pd.crosstab(exclusifs["organe"], exclusifs["est_pathogene"])
        ),
        "profil_types": khi_deux(
            pd.crosstab(
                pathogenes["organe"], regrouper_rares(pathogenes, "type", NB_TYPES_TEST)
            )
        ),
        "profil_genes": khi_deux(
            pd.crosstab(
                pathogenes["organe"], regrouper_rares(pathogenes, "gene", NB_GENES_TEST)
            )
        ),
    }


def lire_comptages_shell(chemin):
    """Lit comptages_shell.txt et renvoie {section: {valeur: nombre}}."""
    if not os.path.isfile(chemin):
        erreur(f"fichier introuvable : {chemin}", CODE_FICHIER_MANQUANT)

    comptages = {}
    section = None
    with open(chemin, encoding="utf-8") as fichier:
        for ligne in fichier:
            ligne = ligne.rstrip("\n")
            if ligne.startswith("## "):
                section = ligne[3:]
                comptages[section] = {}
            else:
                # Format de uniq -c : le nombre, une espace, puis la valeur.
                nombre, _, valeur = ligne.lstrip().partition(" ")
                comptages[section][valeur] = int(nombre)
    return comptages


def verifier_comptages(df, chemin):
    """Compare les nombres du shell à ceux de Python et renvoie la liste des écarts."""
    ecarts = []
    for section, valeurs in lire_comptages_shell(chemin).items():
        # Une section qui ne correspond pas à une colonne n'est pas comparable.
        if section not in df.columns:
            continue
        comptes_python = df[section].value_counts()
        for valeur, nombre_shell in valeurs.items():
            cle = normaliser_organe(valeur) if section == "organe" else valeur
            nombre_python = int(comptes_python.get(cle, 0))
            if nombre_python != nombre_shell:
                ecarts.append(
                    f"{section} « {valeur} » : shell = {nombre_shell}, "
                    f"python = {nombre_python}"
                )
        # Le shell ne liste que les premiers gènes ; le reste doit être complet.
        if section != "gene" and sum(valeurs.values()) != len(df):
            ecarts.append(
                f"{section} : total shell = {sum(valeurs.values())}, "
                f"total python = {len(df)}"
            )
    return ecarts


def lire_arguments():
    """Définit et lit les arguments de la ligne de commande."""
    parser = argparse.ArgumentParser(
        description="Calcule les statistiques des variants du sein et de l'ovaire."
    )
    parser.add_argument("--input", required=True, help="fichier filtré (TSV)")
    parser.add_argument("--output", required=True, help="fichier de sortie stats.csv")
    parser.add_argument(
        "--classification",
        default="Pathogenic",
        help="mot recherché dans la classification, sans tenir compte de la casse "
        "(défaut : Pathogenic, qui inclut Likely pathogenic)",
    )
    parser.add_argument(
        "--comptages-shell",
        help="comptages_shell.txt à comparer aux résultats ; un écart arrête le script",
    )
    return parser.parse_args()


def main():
    args = lire_arguments()

    lus = charger_variants(args.input)
    classes = retirer_non_classes(lus)
    if classes.empty:
        erreur(f"aucun variant classé dans {args.input}", CODE_FICHIER_MANQUANT)
    df = marquer_pathogenes(classes, args.classification)

    if args.comptages_shell:
        ecarts = verifier_comptages(df, args.comptages_shell)
        if ecarts:
            erreur("le shell et Python ne donnent pas les mêmes nombres :\n" + "\n".join(ecarts))
        print("Contrôle : le shell et Python donnent les mêmes nombres.")

    stats = construire_stats(df)
    tests = comparer_organes(df)
    n_pathogenes = int(df["est_pathogene"].sum())

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    stats.to_csv(args.output, index=False, float_format="%.4f")

    chemin_json = os.path.splitext(args.output)[0] + ".json"
    with open(chemin_json, "w", encoding="utf-8") as fichier:
        json.dump(
            {
                "fichier_entree": args.input,
                "motif_classification": args.classification,
                "lignes_lues": len(lus),
                "lignes_retirees_non_classees": len(lus) - len(classes),
                "variants_analyses": len(df),
                "variants_pathogenes": n_pathogenes,
                "tests_sein_seul_vs_ovaire_seul": tests,
            },
            fichier,
            indent=2,
            ensure_ascii=False,
        )

    for nom, test in tests.items():
        if test["statut"] == "calculé" and not test["fiable"]:
            print(f"ATTENTION : test {nom} peu fiable (effectifs attendus < 5).", file=sys.stderr)

    print(
        f"OK : {len(df)} variants analysés, {n_pathogenes} pathogènes "
        f"({n_pathogenes / len(df):.1%}). Résultats : {args.output}, {chemin_json}"
    )


if __name__ == "__main__":
    main()
