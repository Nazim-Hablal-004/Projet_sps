# Partie 5 : figures, archivage, planification

## Installation de l'environnement Python

    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

Les versions exactes sont figées dans `requirements.txt` pour que le projet
se comporte de la même façon sur toutes les machines. Développé avec
Python 3.14 (macOS). Sur un système plus ancien, si une version ne s'installe
pas, assouplir la ligne concernée.

## Figures : `scripts/figures.py`

    python3 scripts/figures.py --input results/stats.csv --outdir results/figures

Lit `stats.csv` (produit par l'étape d'analyse) et crée deux PNG :
- `figure1_genes.png` : proportion de variants pathogènes par gène.
- `figure2_types.png` : proportion de pathogènes par type de mutation, sein
  contre ovaire. Seuls les types ayant au moins 10 variants dans chaque organe
  sont affichés, car une proportion calculée sur très peu de variants est peu
  fiable.

Le script utilise le moteur `Agg` de matplotlib (pas de fenêtre), il fonctionne
donc sans écran, par exemple dans cron. Si le fichier d'entrée est absent, il
affiche une erreur et se termine avec le code 1.

## Archivage : `scripts/archive.sh`

    ./scripts/archive.sh [dossier_results]

Crée `results/archives/resultats_AAAA-MM-JJ_HHMMSS.tar.gz`, avec un fichier
`.sha256` à côté pour vérifier l'intégrité, puis contrôle que l'archive est
lisible. Le dossier `archives/` est exclu de l'archive. Le script peut être
relancé sans nettoyage.

Vérification de l'empreinte :

    cd results/archives
    sha256sum -c resultats_...tar.gz.sha256     (shasum -a 256 -c sous macOS)

## Planification : `cron/pipeline.cron`

Une ligne qui lance le pipeline chaque lundi à 03h00 et écrit la sortie et les
erreurs dans `logs/cron.log`. Le chemin après `cd` est celui de la machine de
développement : l'adapter avant l'installation.

    crontab cron/pipeline.cron     # installer
    crontab -l                     # vérifier

## Droits d'accès

    chmod +x run.sh scripts/*.sh scripts/*.py

Les scripts doivent être exécutables (`755`). Les données brutes restent en
lecture seule (`444`) pour éviter toute modification accidentelle.
