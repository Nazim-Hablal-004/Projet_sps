# Projet « Enquête sur une maladie » — ClinVar, cancer du sein et de l'ovaire

**Question :** parmi les variants ClinVar associés au cancer du sein et de l'ovaire,
quels gènes et quels types de mutations sont le plus souvent pathogènes, et le
profil diffère-t-il entre les deux organes ?

## 1. Installation

Outils système : `bash`, `wget` ou `curl`, `gzip`, `awk`, `md5sum`, `tar`, `python3`.

```bash
sudo apt install python3-venv    # Ubuntu/Debian : venv est un paquet à part
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
chmod +x run.sh scripts/*.sh scripts/*.py
```

## 2. Mode d'emploi

```bash
./run.sh                         # tout le pipeline, options par défaut
./run.sh -g BRCA1,BRCA2 -o sein  # seulement BRCA1/BRCA2, cancer du sein
./run.sh -h                      # aide
```

| Option | Rôle | Défaut |
|---|---|---|
| `-h` | affiche l'aide | |
| `-g GENES` | gènes séparés par des virgules | tous |
| `-o ORGANE` | `sein` ou `ovaire` : garde cet organe et les variants communs aux deux ; `les_deux` : tout | `les_deux` |
| `-c CLASSIF` | classification étudiée | `Pathogenic` |
| `-v` | mode verbeux (messages DEBUG) | |
| `-t` | lance aussi `pytest tests/` | |

**Codes de retour :** 0 succès · 1 erreur générale · 2 mauvaise option ·
3 fichier manquant · 4 somme de contrôle invalide · 130 Ctrl+C · 143 arrêt (TERM).

Le script peut être relancé sans rien nettoyer : les dossiers sont créés avec
`mkdir -p`, les données déjà téléchargées et valides ne sont pas re-téléchargées,
et le dossier temporaire est supprimé automatiquement (`trap`). Sans réseau, le
fichier déjà téléchargé est réutilisé s'il correspond toujours à la somme de
contrôle enregistrée dans `data/raw/variant_summary.txt.gz.md5`.

**Version des données.** Le NCBI met à jour `variant_summary.txt.gz` régulièrement
(environ chaque semaine). Les résultats fournis dans `results/` ont été obtenus avec
la version du 13/11/2026 (MD5 `d2bb73808fcfbcc1f17ac243f5ef994d`, voir `SOURCES.md`).
Un lancement plus tardif télécharge la version du moment : les chiffres peuvent
alors différer légèrement.

## 3. Arborescence

```
Projet_Maladie/
├── run.sh                  chef d'orchestre (P1)
├── scripts/
│   ├── lib.sh              fonctions communes : log, die, check_file, cleanup (P1)
│   ├── download.sh         téléchargement + md5 + lecture seule (P2)
│   ├── check_data.py       vérification des colonnes (P2)
│   ├── filter.sh           filtrage sein/ovaire, GRCh38 (P3)
│   ├── quick_counts.sh     comptages rapides par tubes (P4)
│   ├── analyse.py          statistiques -> stats.csv + stats.json (P4)
│   ├── figures.py          graphiques PNG (P5)
│   ├── archive.sh          archive datée (P5)
│   └── resume.py           résumé JSON (P1)
├── tests/                  mini_clinvar.tsv, test_filter.py (P3), test_analyse.py (P4)
├── data/raw/               données brutes, lecture seule (444)
├── data/clean/             variants_filtres.tsv
├── results/                stats.csv, resume.json, figures/, archives/
├── logs/                   pipeline_AAAA-MM-JJ.log (journal horodaté)
├── cron/                   pipeline.cron (P5)
├── rapport.pdf             rapport du projet (4 pages)
├── requirements.txt        (P5)
├── SOURCES.md              (P2)
└── LIEN_GIT.txt
```

### Fichiers produits

| Fichier | Produit par | Contenu |
|---|---|---|
| `data/raw/variant_summary.txt.gz` | download.sh | fichier ClinVar brut |
| `data/raw/variant_summary.txt.gz.md5` | download.sh | somme de contrôle du NCBI, pour revérifier le fichier sans réseau |
| `data/clean/colonnes.txt` | check_data.py | numéro des colonnes utiles |
| `data/clean/variants_filtres.tsv` | filter.sh | variation_id, gene, type, classification, organe |
| `results/comptages_shell.txt` | quick_counts.sh | comptages par type, gène, classification, organe |
| `results/stats.csv` | analyse.py | section, categorie, organe, n_variants, n_pathogenes, proportion_pathogene, part_parmi_pathogenes |
| `results/stats.json` | analyse.py | tests du khi-deux sein vs ovaire |
| `results/figures/figure1_genes.png` | figures.py | proportion de pathogènes par gène |
| `results/figures/figure2_types.png` | figures.py | sein vs ovaire par type de mutation |
| `results/resume.json` | resume.py | chiffres clés, date, options |
| `results/archives/resultats_*.tar.gz` | archive.sh | copie datée des résultats |
| `logs/pipeline_AAAA-MM-JJ.log` | lib.sh (log) | journal horodaté |

### Figures et archives

- `figure1_genes.png` : proportion de pathogènes pour les 15 gènes qui comptent le plus
  de variants pathogènes.
- `figure2_types.png` : proportion de pathogènes par type de mutation, sein contre ovaire.
  Seuls les types d'au moins 10 variants dans chaque organe sont affichés. Avec
  `-o sein` ou `-o ovaire`, cette figure n'est pas produite et l'ancienne est supprimée,
  pour ne pas l'archiver avec des résultats qui ne lui correspondent pas.
- Chaque exécution crée `results/archives/resultats_AAAA-MM-JJ_HHMMSS.tar.gz` et son
  empreinte `.sha256`. Pour vérifier une archive :
  `cd results/archives && sha256sum -c resultats_AAAA-MM-JJ_HHMMSS.tar.gz.sha256`

### Motif recherché (méthode)

- `Assembly` = GRCh38 seulement (sinon chaque variant est compté deux fois) ;
- `PhenotypeList` découpée sur `|` et `;` : on garde le variant si **une** maladie
  contient un mot d'organe (breast, ovarian, ovary) **et** un mot de cancer
  (cancer, carcinoma, neoplasm, tumor, malignant) ;
- les classifications vides, `-`, `not provided`, `no classification…` sont retirées ;
  les autres (Benign, VUS…) sont gardées pour calculer des proportions ;
- un variant est pathogène si sa classification contient le mot `Pathogenic`
  (donc aussi `Likely pathogenic`, mais pas `Conflicting … pathogenicity`).

## 4. Tâche cron et droits

```bash
# remplacer /chemin/vers/Projet_sps dans cron/pipeline.cron, puis :
crontab cron/pipeline.cron     # installer
crontab -l                     # vérifier
```

Le pipeline se lance chaque lundi à 3h00, la sortie va dans `logs/cron.log`.
Les 5 champs `0 3 * * 1` : minute, heure, jour du mois, mois, jour de la semaine.

Droits : scripts exécutables (`755`), données brutes en lecture seule (`444`,
posé par download.sh) pour qu'elles ne soient jamais modifiées par erreur.

## 5. Tests

```bash
source .venv/bin/activate
pytest tests/        # ou : ./run.sh -t
```

- `tests/test_filter.py` (P3) lance `filter.sh` sur `tests/mini_clinvar.tsv`
  (20 lignes inventées dont on connaît le résultat) : nombre de lignes, aucun
  GRCh37, maladies non cancéreuses rejetées, organes, filtre par gène ;
- `tests/test_analyse.py` (P4) vérifie les proportions, le khi-deux, les erreurs
  et que `quick_counts.sh` et `analyse.py` donnent les mêmes nombres.

## 6. Répartition du travail

| Personne | Membre | Rôle | Fichiers |
|---|---|---|---|
| P1 | HABLAL Nazim | Orchestration et robustesse | run.sh, scripts/lib.sh, scripts/resume.py, README |
| P2 | GAHFIF Djihane | Acquisition et intégrité | download.sh, check_data.py, SOURCES.md |
| P3 | ACHERIR Ghilene | Nettoyage et filtrage | filter.sh, tests/ |
| P4 | INEZARENE Hamed Abdelaziz | Analyse statistique | quick_counts.sh, analyse.py, test_analyse.py |
| P5 | SEKAIMA Toufik | Figures, archivage, système | figures.py, archive.sh, cron, requirements.txt |

Le rapport (`rapport.pdf`) a été rédigé en commun.
