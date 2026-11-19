from collections import Counter
from pathlib import Path
import csv
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MINI = ROOT / "tests" / "mini_clinvar.tsv"
FILTER = ROOT / "scripts" / "filter.sh"


def run_filter(tmp_path, genes=None, origin=None):
    out = tmp_path / "variants_filtres.tsv"
    cmd = ["bash", str(FILTER), "-i", str(MINI), "-o", str(out)]
    if genes:
        cmd += ["-g", genes]
    if origin:
        cmd += ["-r", origin]
    subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    return out


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# --- Tests basiques ---

def test_nombre_de_lignes_filtrees(tmp_path):
    rows = read_rows(run_filter(tmp_path))
    assert len(rows) == 11


def test_aucune_ligne_grch37(tmp_path):
    rows = read_rows(run_filter(tmp_path))
    assert all(r["gene"] != "GRCH37ONLY" for r in rows)


def test_colonnes_attendues(tmp_path):
    rows = read_rows(run_filter(tmp_path))
    assert rows, "aucune ligne en sortie"
    attendu = {"variation_id", "gene", "type", "classification", "organe"}
    assert set(rows[0].keys()) == attendu


def test_variation_id_conserve(tmp_path):
    """Correction 1 : 1015 et 1016 ont tout identique sauf VariationID."""
    rows = read_rows(run_filter(tmp_path))
    ids = {r["variation_id"] for r in rows}
    assert "1015" in ids
    assert "1016" in ids


def test_classification_tiret_retiree(tmp_path):
    """Correction 2 : la classification '-' est retirée."""
    rows = read_rows(run_filter(tmp_path))
    for r in rows:
        assert r["classification"] != "-"
        assert not r["classification"].startswith("no classification")


def test_maladies_non_cancereuses_rejetees(tmp_path):
    """Correction 3 : Premature ovarian failure et Ovarian dysgenesis rejetés."""
    rows = read_rows(run_filter(tmp_path))
    ids = {r["variation_id"] for r in rows}
    assert "1009" not in ids   # Ovarian dysgenesis
    assert "1017" not in ids   # Premature ovarian failure


def test_deux_maladies_separees(tmp_path):
    """Correction 3 : 1019 a 'Breast cancer; Ovarian cancer' -> les_deux."""
    rows = read_rows(run_filter(tmp_path))
    row = next(r for r in rows if r["variation_id"] == "1019")
    assert row["organe"] == "les_deux"


def test_organes(tmp_path):
    rows = read_rows(run_filter(tmp_path))
    assert Counter(r["organe"] for r in rows) == {
        "sein": 6,
        "ovaire": 3,
        "les_deux": 2,
    }

def test_classifications_non_pathogenes_gardees(tmp_path):
    """Correction 4 : Benign et VUS doivent apparaître en sortie."""
    rows = read_rows(run_filter(tmp_path))
    classes = {r["classification"] for r in rows}
    assert "benign" in classes
    assert "uncertain significance" in classes


def test_filtre_genes(tmp_path):
    rows = read_rows(run_filter(tmp_path, genes="BRCA1,BRCA2"))
    assert {r["gene"] for r in rows} <= {"BRCA1", "BRCA2"}


# --- Tests de l'option -r germline ---

def test_germline_seulement(tmp_path):
    rows = read_rows(run_filter(tmp_path, origin="germline"))
    # 1018 est somatic -> rejeté
    ids = {r["variation_id"] for r in rows}
    assert "1018" not in ids
    assert len(rows) == 10


def test_toutes_origines_par_defaut(tmp_path):
    rows = read_rows(run_filter(tmp_path))
    ids = {r["variation_id"] for r in rows}
    assert "1018" in ids   # somatic conservé par défaut
