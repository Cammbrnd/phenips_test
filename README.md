# Mise en prod mutualisée — PHENIPS, CHAPY, STENO, Monochamus

## Arborescence finale du repo (remplace complètement l'ancienne)

```
.github/workflows/
  modeles-update.yml          <- remplace phenips-update.yml (a supprimer)

R/
  safran_common.R             <- fonctions partagees (telechargement/cache/assemblage)
  fetch_safran.R               <- etape unique de recuperation SAFRAN

  phenips_functions.R          <- modele + export GeoJSON PHENIPS (deja allege)
  pipeline_phenips.R           <- driver PHENIPS (charge le jeu partage)

  chapy_functions.R            <- modele + export GeoJSON CHAPY (deja allege)
  pipeline_chapy.R             <- driver CHAPY

  steno_functions.R            <- modele + export GeoJSON STENO (deja allege)
  pipeline_steno.R             <- driver STENO (grille desactivee -> 2 fichiers)

  monochamus_functions.R       <- modele + export GeoJSON Monochamus (allege +
                                   patche : bivoltisme A ajoute au fichier stade)
  pipeline_monochamus.R        <- driver Monochamus

app/
  app.R                        <- app Shiny multi-modeles (selecteur)

data/                          <- genere automatiquement par le workflow
  cache_safran/                <- cache CSV partage (a mettre dans .gitignore)
  cache_meteo/                 <- rasters nationaux assembles (a mettre dans .gitignore)
  phenips_generations_latest.geojson
  phenips_onset_latest.geojson
  chapy_generations_latest.geojson
  chapy_envol_latest.geojson
  steno_generations_latest.geojson
  steno_essaimage_latest.geojson
  mono_stade_latest.geojson          (porte aussi le bivoltisme A)
  mono_bivoltB_gen2_latest.geojson   (bivoltisme B)
```

## .gitignore a ajouter/completer a la racine

```
data/cache_safran/
data/cache_meteo/
```

## Etapes de mise en route

1. **Pousser tous les fichiers ci-dessus** dans le repo `phenips_test`
   (remplacer l'ancien `.github/workflows/phenips-update.yml` par
   `modeles-update.yml`, supprimer l'ancien si tu veux faire le menage).

2. **Verifier chaque `xxx_functions.R` en le sourçant en local** avant de
   lancer sur GitHub Actions :
   ```r
   source("R/safran_common.R")
   source("R/phenips_functions.R")     # doit se sourcer sans erreur, sans rien lancer
   source("R/chapy_functions.R")
   source("R/steno_functions.R")
   source("R/monochamus_functions.R")
   ```
   Si l'un des quatre plante au sourcing, c'est probablement une fonction
   encore appelee mais supprimee par erreur — regarde le message d'erreur,
   il indiquera le nom de la fonction introuvable.

3. **Test complet en local** (optionnel mais rassurant avant de lancer sur
   Actions) :
   ```r
   source("R/safran_common.R")
   source("R/fetch_safran.R")          # telecharge + sauvegarde data/cache_meteo/meteo_national.rds
   source("R/phenips_functions.R"); source("R/pipeline_phenips.R")
   source("R/chapy_functions.R");   source("R/pipeline_chapy.R")
   source("R/steno_functions.R");   source("R/pipeline_steno.R")
   source("R/monochamus_functions.R"); source("R/pipeline_monochamus.R")
   ```
   Tu dois voir apparaitre les 8 GeoJSON dans `data/`.

4. **Sur GitHub** : onglet **Actions** -> workflow "Mise a jour des 4 modeles
   (GeoJSON)" -> **Run workflow** (test manuel avant de compter sur le cron).
   Prevoir un premier run plus long (~30-60 min : telechargement complet +
   4 modelisations), les suivants seront plus rapides grace aux caches.

5. **App Shiny** : `UTILISATEUR_GH`/`REPO_GH`/`BRANCHE_GH` sont deja remplis
   avec tes valeurs (`Cammbrnd`/`phenips_test`/`main`). Teste en local :
   ```r
   shiny::runApp("app")
   ```
   Puis deploie :
   ```r
   rsconnect::deployApp("app")
   ```

## Points de vigilance

- **STENO** : ses GeoJSON portent un champ `couleur` pre-calcule (comme
  PHENIPS/CHAPY) ; l'app l'utilise tel quel et en tire la legende.
- **Le patch bivoltisme A** dans `monochamus_functions.R` extrait les valeurs
  par centroide de polygone (pas cellule par cellule) : sur de tres grandes
  zones dissoutes, une petite perte de precision locale est possible aux
  bordures. Acceptable pour un indicateur de risque annuel, mais a garder en
  tete si tu compares avec le raster brut.
- **Premier run plus lent** que les suivants : le cache SAFRAN et le cache de
  packages R doivent se remplir une premiere fois.

## PHENIPS v2 : scenarios, controle qualite, sources

### Fichiers produits (data/)

| Fichier | Contenu |
|---|---|
| `phenips_generations_{sunny,shaded}_latest.geojson` | classes de generations (`-`, `1`, `1s`, `2`, `2s`...) |
| `phenips_indice_{sunny,shaded}_latest.geojson` | indice continu de developpement cumule (classes de 0.25 generation) |
| `phenips_onset_{sunny,shaded}_latest.geojson` | date du premier essaimage (`date_iso`, `date_onset`) |
| `phenips_meta.json` | date des donnees, date de generation, version de barrks, mailles comblees, source |

- `sunny` = ecorce exposee au sud, `shaded` = ecorce a l'ombre (scenario `max`
  de barrks dans les deux cas). Aucun des deux n'est systematiquement le pire cas :
  au-dela de 38.9 C le developpement PHENIPS s'arrete, et une ecorce exposee au sud
  depasse souvent ce seuil en ete. Sur les donnees d'octobre 2026, le scenario a l'ombre
  donne en general une generation de plus que le scenario au soleil.
- `COMPAT_ANCIENS_NOMS <- TRUE` (phenips_functions.R) recopie aussi le scenario
  `sunny` vers les anciens noms `phenips_generations_latest.geojson` et
  `phenips_onset_latest.geojson`, le temps de redeployer l'app. A passer a `FALSE`
  ensuite, puis supprimer ces deux fichiers.
- Les anciens fichiers dates `phenips_*_AAAAMMJJ.geojson` sont purges : on garde
  les 4 dates les plus recentes par type de fichier. Un `*_latest` est toujours
  la copie du fichier ecrit pendant le run (jamais d'un fichier plus ancien) ;
  si la couche n'a pas pu etre produite, le `latest` est supprime.
- **Que compte une « generation » ?** `get_generations_rst(threshold = 0)` compte les generations
  *presentes* (documentation barrks : « how many generations are present »), c'est-a-dire dont le
  developpement a commence ; elles ne sont pas forcement achevees. Les libelles sont donc `Gen. 3`,
  `Gen. 3 + couvee soeur`, sans « achevee ». L'avancement reel est donne par l'indice continu
  (somme des developpements des generations filiales, plafonnes a 1).
- Les cartes montrent un **potentiel thermique de developpement**, sans masque
  d'epicea ni masque du territoire.

### Controle qualite SAFRAN (R/safran_common.R, R/fetch_safran.R)

1. Telechargement par dalles de 100 km avec 16 km de recouvrement ; une dalle en
   erreur est retentee, puis fait echouer le run.
2. Les mailles absentes du cache sont detectees par comparaison avec leurs
   voisines, puis retelechargees par bandes etroites.
3. Les trous residuels de 2 mailles de large au plus sont interpoles.
4. Tout trou restant fait echouer le run (rien n'est publie, ticket GitHub ouvert).
   Pour contourner temporairement : `SAFRAN_QC_STRICT: "false"` dans le workflow
   (les trous sont alors signales par `::warning::`).

Le cache est invalide par le prefixe `safran_v2_` (et la cle `safran-v2-` du workflow) :
le premier run retelecharge tout (~121 requetes).

### Versions des packages (renv)

Lancer **une fois** en local `source("R/init_renv.R")`, puis commiter `renv.lock`,
`renv/activate.R`, `.Rprofile` et `renv/settings.json`. Le workflow utilise
alors exactement ces versions ; sans `renv.lock` il installe les dernieres versions
et affiche un avertissement.

### Sources et mentions

- Donnees meteorologiques : **Meteo-France** (reanalyse SAFRAN-ISBA), distribuees par
  **GeoSAS (INRAE)**. La mention figure dans l'app, dans `phenips_meta.json` et dans le
  champ `source` des GeoJSON PHENIPS.
- Modele : PHENIPS (Baier et al. 2007), implementation du package R **barrks**.
- Licence du depot : **Licence Ouverte 2.0 (Etalab)**, concedant : Departement de la sante
  des forets (DSF) ; voir le fichier `LICENSE`. La mention « Source : Meteo-France » doit
  etre conservee en cas de reutilisation.

### STENO-DJ : trois scenarios (regles de Pineau et al. 2016)

`R/steno_functions.R` et `R/pipeline_steno.R` appliquent les regles de Pineau et al.
(2016) : temperature de l'air, 506 DJ au-dessus de 10.9 C, debut d'activite quand
Tmax > 18 C deux jours de suite, ni maturation ni arret photoperiodique. Les
generations sont des **maxima thermiques**, pas une prediction validee.
Trois scenarios reprennent les marges d'erreur de la these :

| Scenario | Seuil | Constante | Fichiers `latest` |
|---|---|---|---|
| `moyen` (reference) | 10.9 C | 506 DJ | `steno_generations_latest.geojson`, `steno_essaimage_latest.geojson` |
| `pire` (rapide) | 10.7 C | 499 DJ | `steno_*_pire_latest.geojson` |
| `moins_pire` (lent) | 11.1 C | 513 DJ | `steno_*_moins_pire_latest.geojson` |

L'app propose ces trois scenarios dans un selecteur. Les fichiers `pire` et `moins_pire`
n'existent qu'apres le premier run du workflow.

### CHAPY : classes de generations

`-` (aucune generation achevee) est lu comme le code 0 et non plus comme NA : le polygone
gris « Generation NA » disparait. Le code 4.5 (couvee soeur de la generation 4) a son
libelle et sa couleur. Un code inconnu ressort en magenta avec un `::warning::`.

## v3 : masques d'essence (BD Forêt IGN)

Depuis la v3, les modèles ne sont calculés que sur les mailles SAFRAN (8 km) qui **croisent
l'essence-hôte** :

| Modèle | Essence | Masque |
|---|---|---|
| PHENIPS-Clim, CHAPY | sapin / épicéa | `data/masques/masque_epicea.csv` |
| STENO-DJ | tous les pins | `data/masques/masque_pin.csv` |
| PHENTHAU | chênes décidus (« Chênes décidus », hors sempervirents) | `data/masques/masque_chene.csv` |

- Source : **IGN, BD Forêt® V2** (WFS Géoplateforme, couche
  `LANDCOVER.FORESTINVENTORY.V2:formation_vegetale`, champ `essence`). Épicéa = classe
  « Sapin, épicéa » ; pin = toutes les classes commençant par « Pin » (sylvestre, maritime,
  laricio/noir, d'Alep, à crochets/cembro, pins mélangés, autre pin).
- Une maille est gardée dès qu'elle croise au moins un peuplement de l'essence
  (`SEUIL_PART_MASQUE <- 0` dans `R/masques.R` ; le CSV contient aussi `part`, la part de la
  maille couverte, pour monter ce seuil si besoin).
- **Limite** : la BD Forêt V2 attribue une seule essence par peuplement. « Sapin, épicéa »
  désigne les peuplements purs ; l'épicéa présent dans un peuplement classé « Mélange de
  conifères » ou « Conifères » n'est pas visible. Les masques sont donc plutôt prudents.
- Les masques sont **construits à la demande** (la BD Forêt change peu) : onglet Actions →
  « Mise à jour des 4 modèles » → *Run workflow* avec `construire_masques` coché. Le workflow
  commite `data/masques/`. Les pipelines s'arrêtent avec un message clair si un masque manque.
- Monochamus n'est pas masqué.

## Historique 1960-2025 : générations et premier envol par maille

Table annuelle, pour chaque maille SAFRAN (8 km) masquée par l'essence-hôte (cf. v3), du
**nombre de générations** et de la **date du premier envol** de 1960 à 2025.

| Modèle | Fichier | Colonnes (hors `annee`, `x`, `y`) |
|---|---|---|
| PHENIPS-Clim | `data/historique/phenips_annuel.csv.gz` | `gen_sunny`, `doy_sunny`, `gen_shaded`, `doy_shaded` |
| CHAPY | `data/historique/chapy_annuel.csv.gz` | `gen`, `doy_envol` |
| STENO-DJ | `data/historique/steno_annuel.csv.gz` | `gen_pire`, `doy_pire`, `gen_moyen`, `doy_moyen`, `gen_moins_pire`, `doy_moins_pire` |

`data/historique/meta.json` donne la date de construction, la source et, par modèle, les années
couvertes et le nombre de mailles.

- **x, y** : centre de la maille en Lambert-93 (EPSG:2154). Maille de 8 km (x ± 4000 m, y ± 4000 m).
- **gen\*** : nombre de générations au 31 décembre. PHENIPS-Clim et CHAPY : une génération est
  comptée dès qu'elle a commencé (+0,5 pour une couvée sœur) ; STENO-DJ : générations achevées.
- **doy\*** : jour de l'année du premier envol (PHENIPS-Clim : premier essaimage ; CHAPY et
  STENO-DJ : premier envol). Vide si le modèle n'atteint pas l'envol cette année-là.
- **gen = 0 avec une date d'envol renseignée** : le vol de printemps a lieu mais aucune génération
  n'est comptée cette année-là (surtout les années froides : ~3 % des lignes CHAPY dans les années
  1960-70, < 1 % depuis 2020).
- **STENO-DJ, dates précoces** : ~5 % des lignes ont un premier envol avant le jour 60 (journées
  douces de fin d'hiver dépassant les seuils du modèle dans les mailles les plus chaudes) ; à lire
  avec prudence, de préférence via la médiane ou la moyenne sur plusieurs années.
- Chaque année est calculée de janvier à décembre, **sans lien avec l'année précédente** (même
  hypothèse que les cartes de l'année en cours).
- Mêmes scénarios, mêmes masques et même contrôle qualité SAFRAN que les cartes de l'année en
  cours. L'année en cours n'est pas incluse (elle est incomplète).
- Source météo : SAFRAN (Météo-France) via GeoSAS INRAE ; la série GeoSAS démarre en août 1958,
  d'où le début à 1960 (années civiles complètes).

### Reconstruire la table

Onglet Actions → « Mise à jour des 4 modèles » → *Run workflow* avec `historique` coché
(`annee_debut`, `annee_fin`, `parallele` réglables). Un job par année (5 en parallèle par défaut,
pour ne pas surcharger GeoSAS), puis un job d'assemblage qui commite `data/historique/`.
Les années recalculées remplacent les anciennes lignes ; les autres sont conservées.
Le calcul complet 1960-2025 prend environ 6 h (≈ 16 min par année, 5 en parallèle ; la première
exécution est plus longue car elle installe les packages R, ensuite mis en cache). Si un job reste
bloqué (déjà vu une fois à l'installation de R), « Re-run failed jobs » relance uniquement l'année
en cause : les autres résultats sont conservés. Le job
`historique_diag` (case `historique_diag`) sonde GeoSAS (étendue temporelle, durée des requêtes)
sans rien calculer.

Pour l'année en cours, `R/valider_historique.R` compare l'extraction par maille aux GeoJSON de
production (annotations `VALIDATION` du job).

### Onglet « Historique 1960-2025 » de l'app

Carte d'une année (curseur animable), moyenne sur une période, ou tendance par décennie
(régression linéaire maille par maille, grisage des tendances non significatives à p > 0,05) ;
clic sur une maille pour afficher sa série temporelle face à la moyenne de toutes les mailles.
L'app ne fait que lire les tables du dépôt (aucun calcul lourd).

## PHENTHAU : processionnaire du chêne (*Thaumetopoea processionea*)

5e modèle, package CRAN **PHENTHAUproc** (Halbig et al. 2024, *Forest Ecol. Manag.* 552:121525),
mode régional `dailymeanminmax` : températures **journalières** SAFRAN (`T_Q`, `TINF_H_Q`,
`TSUP_H_Q`), sans reconstruction horaire.

- **Calcul** : `R/pipeline_phenthau.R` (fonctions dans `R/phenthau_functions.R`), à la suite des
  4 autres modèles dans `modeles-update.yml` (étapes 7 à 9, `continue-on-error` : un échec
  PHENTHAU ne bloque pas les 4 modèles, déjà publiés, et ouvre le ticket d'échec).
- **Données amont** : le package compte les jours de froid depuis le **1er septembre N-1**.
  `R/fetch_safran_amont.R` récupère septembre-décembre N-1 (même contrôle qualité que le fetch
  partagé, cache `safran_v2_*_<N-1>.csv`, sortie `data/cache_meteo/meteo_amont.rds`) ; le fetch
  partagé (1er janvier N) n'est pas modifié.
- **Sorties** (`data/`) : `phenthau_{stades,l3,eclosion,debourrement,mortalite}_<date>.geojson`
  (+ `_latest`), `phenthau_meta.json`. GeoTIFF L93 (17 bandes de jour de l'année par maille) et PNG
  nationaux : artefact `phenthau-rasters` du run (30 jours), non commités. Champs GeoJSON :
  `espece`, `modele = "PHENTHAU"`, `couleur`, `label`, `date_iso`, `avertissement`.
- **Éclosion** : `DEF_ECLOSION = "first"` (la plus précoce des 3 modèles : Custers, Wagenhoff,
  Meurisse). L'option `mean` du package n'a pas d'effet utile (testée sur 2 zones × 3 ans :
  mêmes dates de stades).
- **Mortalité des L1** : publiée seulement une fois le débourrement atteint (avant, le package
  renvoie 0).
- **Limite d'applicabilité** : calibré sur *Quercus robur* en Allemagne, jusqu'à 800 m. Hors de ce
  domaine (sud de la France, montagne), c'est une extrapolation : au sud-ouest, l'éclosion
  (Meurisse) précède de 3 à 4 semaines le débourrement et la mortalité L1 dépasse 50 % sur presque
  toutes les mailles. À lire comme un indicateur, pas comme une prévision validée.
- **Tests** (entrées *Run workflow*) : `diagnostic_phenthau`, `test_phenthau` (2 petites emprises,
  2023-2025, `def_eclosion` first/mean), `test_amont` (fetch amont national, bilan dans
  `data/diag/`), `masques=chene` avec `construire_masques`.
