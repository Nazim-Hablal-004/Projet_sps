#!/usr/bin/env python3
"""figures.py : dessine les figures à partir de results/stats.csv"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


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


def main():
    args = lire_arguments()

    if not os.path.isfile(args.input):
        print(f"Erreur : fichier introuvable : {args.input}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)
    figure_genes(df, args.outdir)


if __name__ == "__main__":
    main()
