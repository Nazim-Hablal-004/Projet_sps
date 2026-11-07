#!/usr/bin/env python3
"""figures.py : dessine les figures à partir de results/stats.csv"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SEUIL_MIN_VARIANTS = 10


def lire_arguments():
    parser = argparse.ArgumentParser(description="Génère les figures à partir de stats.csv")
    parser.add_argument("--input", required=True, help="chemin de stats.csv")
    parser.add_argument("--outdir", required=True, help="dossier des figures")
    return parser.parse_args()


def figure_genes(df, outdir):
    """Figure 1 : proportion de pathogènes par gène."""
    genes = df[(df["section"] == "gene") & (df["categorie"] != "autres")]
    genes = genes.sort_values("proportion_pathogene", ascending=False)

    plt.figure(figsize=(10, 5))
    plt.bar(genes["categorie"], genes["proportion_pathogene"])
    plt.xlabel("Gène")
    plt.ylabel("Proportion de variants pathogènes")
    plt.title("Proportion de variants pathogènes par gène")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    chemin = os.path.join(outdir, "figure1_genes.png")
    plt.savefig(chemin, dpi=150)
    plt.close()
    print(f"Figure enregistrée : {chemin}")


def figure_types(df, outdir):
    """Figure 2 : proportion de pathogènes par type de mutation, sein vs ovaire."""
    t = df[(df["section"] == "type_organe") & (df["organe"].isin(["sein", "ovaire"]))]
    prop = t.pivot(index="categorie", columns="organe", values="proportion_pathogene")
    effectifs = t.pivot(index="categorie", columns="organe", values="n_variants")

    # on garde les types assez fréquents dans les deux organes
    garder = effectifs[(effectifs["sein"] >= SEUIL_MIN_VARIANTS)
                       & (effectifs["ovaire"] >= SEUIL_MIN_VARIANTS)].index
    prop = prop.loc[garder].sort_values("sein", ascending=False)

    prop[["sein", "ovaire"]].plot(kind="bar", figsize=(10, 5))
    plt.xlabel("Type de mutation")
    plt.ylabel("Proportion de variants pathogènes")
    plt.title("Proportion de variants pathogènes par type de mutation : sein vs ovaire")
    plt.xticks(rotation=30, ha="right")
    plt.legend(title="Organe")
    plt.tight_layout()

    chemin = os.path.join(outdir, "figure2_types.png")
    plt.savefig(chemin, dpi=150)
    plt.close()
    print(f"Figure enregistrée : {chemin}")


def main():
    args = lire_arguments()

    if not os.path.isfile(args.input):
        print(f"Erreur : fichier introuvable : {args.input}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)
    figure_genes(df, args.outdir)
    figure_types(df, args.outdir)


if __name__ == "__main__":
    main()
