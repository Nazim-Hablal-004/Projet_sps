"""Tests de scripts/analyse.py et de scripts/quick_counts.sh (Personne 4).

Lancement depuis la racine du projet : pytest tests/test_analyse.py
"""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

RACINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RACINE / "scripts"))

import analyse  # noqa: E402

ANALYSE = RACINE / "scripts" / "analyse.py"
QUICK_COUNTS = RACINE / "scripts" / "quick_counts.sh"

# Petit fichier dont on connaît les résultats à la main :
# - sein : 3 classés dont 2 pathogènes
# - ovaire : 2 classés dont 1 pathogène
# - les_deux : 2 classés dont 2 pathogènes
# - 2 lignes sans classification, à retirer
FICHIER_FILTRE = """\
variation_id\tgene\ttype\tclassification\torgane
1\tBRCA1\tDeletion\tPathogenic\tles_deux
2\tBRCA2\tsingle nucleotide variant\tLikely pathogenic\tles_deux
3\tATM\tDeletion\tPathogenic/Likely pathogenic\tsein
4\tATM\tsingle nucleotide variant\tPathogenic\tsein
5\tATM\tsingle nucleotide variant\tConflicting classifications of pathogenicity\tsein
6\tTP53\tsingle nucleotide variant\tPathogenic\tovaire
7\tTP53\tsingle nucleotide variant\tBenign\tovaire
8\tTP53\tsingle nucleotide variant\t-\tovaire
9\tATM\tDeletion\tnot provided\tsein
"""


def variants(lignes):
    """Construit un tableau à partir de tuples (gene, type, classification, organe)."""
    return pd.DataFrame(lignes, columns=["gene", "type", "classification", "organe"])


def lancer_analyse(*arguments):
    return subprocess.run(
        [sys.executable, str(ANALYSE), *map(str, arguments)],
        capture_output=True,
        text=True,
    )


@pytest.fixture
def fichier_filtre(tmp_path):
    chemin = tmp_path / "variants_filtres.tsv"
    chemin.write_text(FICHIER_FILTRE, encoding="utf-8")
    return chemin


@pytest.mark.parametrize(
    "classification, attendu",
    [
        ("Pathogenic", True),
        ("Likely pathogenic", True),
        ("Pathogenic/Likely pathogenic", True),
        ("Pathogenic/Likely pathogenic; association", True),
        ("Likely pathogenic, low penetrance", True),
        ("pathogenic", True),
        ("Conflicting classifications of pathogenicity", False),
        ("conflicting classifications of pathogenicity", False),
        ("Conflicting classifications of pathogenicity; risk factor", False),
        ("Benign", False),
        ("Uncertain significance", False),
    ],
)
def test_marquer_pathogenes(classification, attendu):
    df = variants([("BRCA1", "Deletion", classification, "sein")])
    assert bool(analyse.marquer_pathogenes(df)["est_pathogene"].iloc[0]) is attendu


@pytest.mark.parametrize(
    "classification, garde",
    [
        ("-", False),
        ("", False),
        ("not provided", False),
        ("Not Provided", False),
        ("no classifications from unflagged records", False),
        ("Benign", True),
        ("Uncertain significance", True),
    ],
)
def test_retirer_non_classes(classification, garde):
    df = variants([("BRCA1", "Deletion", classification, "sein")])
    assert len(analyse.retirer_non_classes(df)) == (1 if garde else 0)


def test_proportion_pathogene_calculee_a_la_main():
    df = analyse.marquer_pathogenes(
        variants(
            [
                ("BRCA1", "Deletion", "Pathogenic", "sein"),
                ("BRCA1", "Deletion", "Likely pathogenic", "sein"),
                ("BRCA1", "Deletion", "Benign", "sein"),
                ("BRCA2", "Deletion", "Pathogenic", "sein"),
                ("BRCA2", "Deletion", "Benign", "sein"),
                ("BRCA2", "Deletion", "Benign", "sein"),
                ("BRCA2", "Deletion", "Benign", "sein"),
            ]
        )
    )
    table = analyse.proportion_pathogene(df, "gene")
    # Le gène qui a le plus de pathogènes arrive en premier.
    assert list(table["gene"]) == ["BRCA1", "BRCA2"]
    assert list(table["n_variants"]) == [3, 4]
    assert list(table["n_pathogenes"]) == [2, 1]
    assert table["proportion_pathogene"].tolist() == pytest.approx([2 / 3, 1 / 4])


def test_khi_deux_sur_une_table_connue():
    # Effectifs attendus de 15 partout : chi2 = 4 x (5 x 5) / 15 = 6,667.
    resultat = analyse.khi_deux(pd.DataFrame([[10, 20], [20, 10]]))
    assert resultat["statut"] == "calculé"
    assert resultat["chi2"] == pytest.approx(6.667, abs=0.001)
    assert resultat["ddl"] == 1
    assert resultat["v_cramer"] == pytest.approx(1 / 3, abs=0.001)
    assert resultat["effectif_attendu_min"] == 15
    assert resultat["fiable"] is True


def test_khi_deux_signale_un_test_peu_fiable():
    assert analyse.khi_deux(pd.DataFrame([[1, 2], [3, 4]]))["fiable"] is False


def test_tests_non_calculables_avec_un_seul_organe():
    # L'exemple du document : ./run.sh -g BRCA1,BRCA2 -o sein
    df = analyse.marquer_pathogenes(
        variants(
            [
                ("BRCA1", "Deletion", "Pathogenic", "sein"),
                ("BRCA1", "Deletion", "Benign", "sein"),
            ]
        )
    )
    for test in analyse.comparer_organes(df).values():
        assert test["statut"] == "non calculable"


def test_tests_non_calculables_sans_pathogene():
    df = analyse.marquer_pathogenes(
        variants(
            [
                ("BRCA1", "Deletion", "Benign", "sein"),
                ("BRCA1", "Deletion", "Benign", "ovaire"),
            ]
        )
    )
    for test in analyse.comparer_organes(df).values():
        assert test["statut"] == "non calculable"


def test_les_deux_ne_change_pas_les_tests():
    sein_ovaire = [
        ("BRCA1", "Deletion", "Pathogenic", "sein"),
        ("BRCA1", "Deletion", "Benign", "sein"),
        ("BRCA1", "Deletion", "Pathogenic", "ovaire"),
        ("BRCA1", "Deletion", "Benign", "ovaire"),
        ("BRCA1", "Deletion", "Benign", "ovaire"),
    ]
    sans = analyse.marquer_pathogenes(variants(sein_ovaire))
    avec = analyse.marquer_pathogenes(
        variants(sein_ovaire + [("BRCA2", "Deletion", "Pathogenic", "les_deux")] * 50)
    )
    assert analyse.comparer_organes(sans) == analyse.comparer_organes(avec)


def test_autres_arrive_en_dernier():
    # 3 gènes rares totalisent plus de pathogènes que le gène gardé, mais le
    # premier de la section doit rester un vrai gène (utilisé par resume.py).
    lignes = [("BRCA2", "Deletion", "Pathogenic", "sein")] * 3
    for gene in ["G1", "G2", "G3"]:
        lignes += [(gene, "Deletion", "Pathogenic", "sein")] * 2
    df = analyse.marquer_pathogenes(variants(lignes))
    avec_autres = df.assign(gene=analyse.regrouper_rares(df, "gene", 1))
    table = analyse.resumer_par_categorie(avec_autres, "gene", "gene")
    assert list(table["categorie"]) == ["BRCA2", "autres"]
    assert list(table["n_pathogenes"]) == [3, 6]


def test_charger_variants_accepte_les_deux_avec_espace(tmp_path):
    chemin = tmp_path / "f.tsv"
    chemin.write_text(
        "gene\ttype\tclassification\torgane\nBRCA1\tDeletion\tPathogenic\tles deux\n",
        encoding="utf-8",
    )
    assert list(analyse.charger_variants(str(chemin))["organe"]) == ["les_deux"]


def test_script_complet(fichier_filtre, tmp_path):
    sortie = tmp_path / "resultats" / "stats.csv"
    resultat = lancer_analyse("--input", fichier_filtre, "--output", sortie)
    assert resultat.returncode == 0, resultat.stderr

    stats = pd.read_csv(sortie, keep_default_na=False)
    assert list(stats.columns) == analyse.COLONNES_STATS

    def ligne(section, categorie, organe):
        choix = stats[
            (stats["section"] == section)
            & (stats["categorie"] == categorie)
            & (stats["organe"] == organe)
        ]
        assert len(choix) == 1
        return choix.iloc[0]

    # 7 variants classés, 5 pathogènes : la ligne « pathogenicity » n'est pas comptée.
    assert ligne("global", "tous", "tous")["n_variants"] == 7
    assert ligne("global", "tous", "tous")["n_pathogenes"] == 5
    assert ligne("organe", "tous", "sein")["n_variants"] == 3
    assert ligne("organe", "tous", "sein")["n_pathogenes"] == 2
    assert ligne("organe", "tous", "ovaire")["n_variants"] == 2
    assert ligne("organe", "tous", "ovaire")["n_pathogenes"] == 1
    assert ligne("organe", "tous", "les_deux")["n_pathogenes"] == 2
    # Le total d'un organe inclut les variants partagés.
    assert ligne("organe", "tous", "sein_total")["n_variants"] == 5
    assert ligne("organe", "tous", "sein_total")["n_pathogenes"] == 4
    assert ligne("organe", "tous", "ovaire_total")["n_variants"] == 4
    assert ligne("organe", "tous", "ovaire_total")["n_pathogenes"] == 3
    assert ligne("gene", "ATM", "tous")["n_pathogenes"] == 2
    assert ligne("type", "Deletion", "tous")["proportion_pathogene"] == 1
    # Profil : 1 des 2 pathogènes du sein est une délétion.
    assert ligne("type_organe", "Deletion", "sein")["part_parmi_pathogenes"] == 0.5
    assert ligne("gene_organe", "TP53", "ovaire")["n_variants"] == 2

    metadonnees = json.loads(sortie.with_suffix(".json").read_text(encoding="utf-8"))
    assert metadonnees["lignes_lues"] == 9
    assert metadonnees["lignes_retirees_non_classees"] == 2
    assert metadonnees["variants_pathogenes"] == 5
    assert set(metadonnees["tests_sein_seul_vs_ovaire_seul"]) == {
        "proportion_pathogene",
        "profil_types",
        "profil_genes",
    }


def test_script_fichier_absent(tmp_path):
    resultat = lancer_analyse(
        "--input", tmp_path / "absent.tsv", "--output", tmp_path / "stats.csv"
    )
    assert resultat.returncode == analyse.CODE_FICHIER_MANQUANT
    assert "introuvable" in resultat.stderr


@pytest.mark.parametrize(
    "contenu",
    ["", "\n", "variation_id\tgene\ttype\tclassification\torgane\n"],
    ids=["vide", "ligne_vide", "entete_seule"],
)
def test_script_fichier_sans_donnees(tmp_path, contenu):
    chemin = tmp_path / "f.tsv"
    chemin.write_text(contenu, encoding="utf-8")
    resultat = lancer_analyse("--input", chemin, "--output", tmp_path / "stats.csv")
    assert resultat.returncode == analyse.CODE_FICHIER_MANQUANT
    assert "ERREUR" in resultat.stderr


def test_script_colonne_absente(tmp_path):
    chemin = tmp_path / "f.tsv"
    chemin.write_text("gene\ttype\nBRCA1\tDeletion\n", encoding="utf-8")
    resultat = lancer_analyse("--input", chemin, "--output", tmp_path / "stats.csv")
    assert resultat.returncode == analyse.CODE_ERREUR
    assert "classification" in resultat.stderr


def test_script_organe_inconnu(tmp_path):
    chemin = tmp_path / "f.tsv"
    chemin.write_text(
        "gene\ttype\tclassification\torgane\nBRCA1\tDeletion\tPathogenic\tpoumon\n",
        encoding="utf-8",
    )
    resultat = lancer_analyse("--input", chemin, "--output", tmp_path / "stats.csv")
    assert resultat.returncode == analyse.CODE_ERREUR
    assert "poumon" in resultat.stderr


def test_avertissement_sans_variation_id(tmp_path):
    chemin = tmp_path / "f.tsv"
    chemin.write_text(
        "gene\ttype\tclassification\torgane\nBRCA1\tDeletion\tPathogenic\tsein\n",
        encoding="utf-8",
    )
    resultat = lancer_analyse("--input", chemin, "--output", tmp_path / "stats.csv")
    assert resultat.returncode == 0
    assert "variation_id" in resultat.stderr


def test_ecart_entre_shell_et_python(fichier_filtre, tmp_path):
    comptages = tmp_path / "comptages.txt"
    comptages.write_text(
        "## type\n     99 Deletion\n      5 single nucleotide variant\n", encoding="utf-8"
    )
    resultat = lancer_analyse(
        "--input", fichier_filtre,
        "--output", tmp_path / "stats.csv",
        "--comptages-shell", comptages,
    )
    assert resultat.returncode == analyse.CODE_ERREUR
    assert "shell = 99" in resultat.stderr


def test_valeur_oubliee_par_le_shell(fichier_filtre, tmp_path):
    # La valeur listée est juste, mais le total ne correspond pas.
    comptages = tmp_path / "comptages.txt"
    comptages.write_text("## type\n      2 Deletion\n", encoding="utf-8")
    resultat = lancer_analyse(
        "--input", fichier_filtre,
        "--output", tmp_path / "stats.csv",
        "--comptages-shell", comptages,
    )
    assert resultat.returncode == analyse.CODE_ERREUR
    assert "total shell" in resultat.stderr


def test_section_shell_sans_colonne_ignoree(fichier_filtre, tmp_path):
    comptages = tmp_path / "comptages.txt"
    comptages.write_text("## type du sein\n      2 Deletion\n", encoding="utf-8")
    resultat = lancer_analyse(
        "--input", fichier_filtre,
        "--output", tmp_path / "stats.csv",
        "--comptages-shell", comptages,
    )
    assert resultat.returncode == 0, resultat.stderr


@pytest.mark.skipif(
    not (RACINE / "scripts" / "lib.sh").exists(),
    reason="scripts/lib.sh (Personne 1) n'est pas encore disponible",
)
def test_shell_et_python_donnent_les_memes_nombres(fichier_filtre, tmp_path):
    comptages = tmp_path / "comptages.txt"
    shell = subprocess.run(
        ["bash", str(QUICK_COUNTS), str(fichier_filtre), str(comptages)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert shell.returncode == 0, shell.stderr
    assert "      2 Deletion" in comptages.read_text(encoding="utf-8")

    resultat = lancer_analyse(
        "--input", fichier_filtre,
        "--output", tmp_path / "stats.csv",
        "--comptages-shell", comptages,
    )
    assert resultat.returncode == 0, resultat.stderr
    assert "mêmes nombres" in resultat.stdout
