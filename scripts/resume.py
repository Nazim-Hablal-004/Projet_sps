#!/usr/bin/env python3
"""resume.py : résumé du pipeline (chiffres clés de stats.csv) dans un fichier JSON."""

import argparse
import csv
import json
import os
import sys
from datetime import datetime

COLONNES = ["section", "categorie", "organe", "n_variants", "n_pathogenes",
            "proportion_pathogene"]
AUTRES = "autres"


def lire_arguments():
    """Lit les arguments de la ligne de commande."""
    parser = argparse.ArgumentParser(description="Résumé JSON du pipeline ClinVar")
    parser.add_argument("--stats", required=True, help="chemin de stats.csv")
    parser.add_argument("--log", required=True, help="chemin du journal")
    parser.add_argument("--output", required=True, help="chemin du JSON à écrire")
    parser.add_argument("--genes", default="", help="gènes demandés (option -g)")
    parser.add_argument("--organe", default="les_deux", help="organe (option -o)")
    parser.add_argument("--classification", default="Pathogenic",
                        help="classification (option -c)")
    return parser.parse_args()


def lire_stats(chemin):
    """Renvoie la liste des lignes de stats.csv sous forme de dictionnaires."""
    with open(chemin, newline="", encoding="utf-8") as f:
        lecteur = csv.DictReader(f)
        manquantes = [c for c in COLONNES if c not in (lecteur.fieldnames or [])]
        if manquantes:
            print(f"Erreur : colonnes absentes de {chemin} : {manquantes}",
                  file=sys.stderr)
            sys.exit(1)
        lignes = []
        for ligne in lecteur:
            ligne["n_variants"] = int(ligne["n_variants"])
            ligne["n_pathogenes"] = int(ligne["n_pathogenes"])
            ligne["proportion_pathogene"] = float(ligne["proportion_pathogene"])
            lignes.append(ligne)
    if not lignes:
        print(f"Erreur : {chemin} ne contient aucune donnée", file=sys.stderr)
        sys.exit(1)
    return lignes


def filtrer(lignes, section):
    """Garde les lignes d'une section (global, organe, gene, type...)."""
    return [l for l in lignes if l["section"] == section]


def chiffres(ligne):
    """Les trois chiffres utiles d'une ligne de stats.csv."""
    return {"n_variants": ligne["n_variants"],
            "n_pathogenes": ligne["n_pathogenes"],
            "proportion_pathogene": round(ligne["proportion_pathogene"], 4)}


def plus_pathogene(lignes):
    """Catégorie avec le plus de pathogènes, sans compter « autres »."""
    lignes = [l for l in lignes if l["categorie"] != AUTRES]
    if not lignes:
        return None
    meilleure = max(lignes, key=lambda l: l["n_pathogenes"])
    return {"nom": meilleure["categorie"], **chiffres(meilleure)}


def lire_tests(chemin_stats):
    """Reprend les tests du khi-deux de stats.json s'il existe (sinon None)."""
    chemin = os.path.splitext(chemin_stats)[0] + ".json"
    if not os.path.isfile(chemin):
        return None
    with open(chemin, encoding="utf-8") as f:
        return json.load(f).get("tests_sein_seul_vs_ovaire_seul")


def compter_log(chemin):
    """Compte les avertissements et erreurs du journal (0 si absent)."""
    compte = {"lignes": 0, "avertissements": 0, "erreurs": 0}
    if not os.path.exists(chemin):
        return compte
    with open(chemin, encoding="utf-8") as f:
        for ligne in f:
            compte["lignes"] += 1
            if "[ATTENTION]" in ligne:
                compte["avertissements"] += 1
            if "[ERREUR]" in ligne:
                compte["erreurs"] += 1
    return compte


def construire_resume(lignes, args):
    """Rassemble toutes les informations dans un dictionnaire."""
    total = filtrer(lignes, "global")
    return {
        "date_execution": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "options": {
            "genes": args.genes.split(",") if args.genes else "tous",
            "organe": args.organe,
            "classification": args.classification,
        },
        "total": chiffres(total[0]) if total else None,
        "pathogenes_par_organe": {
            l["organe"]: chiffres(l) for l in filtrer(lignes, "organe")
        },
        "gene_le_plus_touche": plus_pathogene(filtrer(lignes, "gene")),
        "type_le_plus_frequent": plus_pathogene(filtrer(lignes, "type")),
        "tests_sein_vs_ovaire": lire_tests(args.stats),
        "journal": compter_log(args.log),
        "fichiers": {"stats": args.stats, "log": args.log},
    }


def main():
    args = lire_arguments()
    if not os.path.isfile(args.stats):
        print(f"Erreur : fichier introuvable : {args.stats}", file=sys.stderr)
        sys.exit(1)

    resume = construire_resume(lire_stats(args.stats), args)

    dossier = os.path.dirname(args.output)
    if dossier:
        os.makedirs(dossier, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(resume, f, indent=2, ensure_ascii=False)
    print(f"Résumé écrit dans {args.output}")


if __name__ == "__main__":
    main()
