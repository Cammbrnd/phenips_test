# Phénologie des ravageurs forestiers : cartes SAFRAN

Cartes nationales (mailles SAFRAN de 8 km, Lambert-93) du développement thermique de cinq
ravageurs. Les données météo viennent de SAFRAN (Météo-France) via GeoSAS (INRAE), avec un
décalage d'environ 5 jours. Un workflow GitHub Actions calcule les cartes chaque lundi
(05h00 UTC) et publie des GeoJSON dans `data/` ; l'app Shiny/leaflet (`app/app.R`) les lit
directement sur GitHub.

| Modèle | Espèce | Masque d'essence | Implémentation |
|---|---|---|---|
| PHENIPS-Clim | *Ips typographus* | épicéa | package **barrks** (`phenology()`) |
| CHAPY | *Pityogenes chalcographus* | épicéa | package **barrks** (`phenology()`) |
| STENO-DJ | *Ips sexdentatus* | pins | interne (règles de Pineau et al.) |
| Monochamus | *Monochamus galloprovincialis* | aucun | interne |
| PHENTHAU | *Thaumetopoea processionea* | chênes décidus | package **PHENTHAUproc** |

## Arborescence

```
.github/workflows/modeles-update.yml   workflow principal (cron hebdomadaire + lancements manuels)
.github/workflows/cloture-annuelle.yml  lance l'historique de l'année écoulée le 8 janvier

R/
  safran_common.R          période, contrôle qualité et cache SAFRAN, constantes partagées
  fetch_safran.R           récupération SAFRAN de l'année (-> data/cache_meteo/meteo_national.rds)
  fetch_safran_amont.R     SAFRAN septembre-décembre N-1 (PHENTHAU)
  fichiers_dates.R         copie vers les noms stables *_latest et purge des fichiers datés
  masques.R, construire_masques.R   masques d'essence (BD Forêt IGN)
  <modele>_functions.R     modèle + export GeoJSON (phenips, chapy, steno, monochamus, phenthau)
  pipeline_<modele>.R      pilote de chaque modèle
  historique_*.R, assembler_historique.R, valider_historique.R   historique 1960+
  diag_*.R, init_renv.R    diagnostics ; figement des versions de packages (renv)

app/app.R                  app Shiny multi-modèles

data/                      généré par le workflow
  <modele>_..._latest.geojson      noms stables lus par l'app
  <modele>_..._AAAAMMJJ.geojson    fichiers datés (les 4 dernières dates)
  masques/  historique/  diag/
  cache_safran/  cache_meteo/      dans .gitignore
```

## Fonctionnement

- **Calendrier** : cron le lundi à 05h00 UTC, ou *Run workflow* à la main (onglet Actions,
  workflow « Mise à jour des modèles (GeoJSON) »). Le premier run est long (téléchargement
  complet de SAFRAN, timeout du job : 180 min) ; les suivants s'appuient sur le cache.
- **Année traitée** : celle du dernier jour SAFRAN disponible (J-5), du 1er janvier à J-5. La
  variable `ANNEE_SAFRAN` force une année passée (historique). Voir « Clôture annuelle ».
- **Un processus par modèle**, lancé depuis la racine du dépôt, comme dans le workflow. Chaque
  `pipeline_<modele>.R` charge `R/safran_common.R` puis son fichier `<modele>_functions.R`.
- **Fichiers de sortie** : chaque pipeline écrit un fichier daté `nom_AAAAMMJJ.geojson`, le
  copie vers `nom_latest.geojson` (le fichier de date la plus récente *dans son nom*, pas la
  date de modification), puis ne garde que les **4 dates les plus récentes** par type de
  fichier. Les suppressions sont publiées avec le commit du workflow. Cas particuliers :
  Monochamus ajoute `_WGS84` au nom daté ; l'essaimage STENO est un fichier par année
  (`steno_essaimage_AAAA.geojson`, non purgé).
- **Publication** : un commit unique pour PHENIPS, CHAPY, STENO et Monochamus (étape 6) ; PHENTHAU
  vient après, dans un commit séparé et en `continue-on-error` (étapes 7 à 9). Un échec ouvre un
  ticket GitHub automatique.

### Reprise en local

```
Rscript R/fetch_safran.R            # SAFRAN -> data/cache_meteo/meteo_national.rds
Rscript R/pipeline_phenips.R        # puis pipeline_chapy.R, pipeline_steno.R, pipeline_monochamus.R
Rscript R/fetch_safran_amont.R      # PHENTHAU : SAFRAN sept-déc N-1
Rscript R/pipeline_phenthau.R
shiny::runApp("app")                # l'app lit les fichiers sur GitHub (UTILISATEUR_GH / REPO_GH / BRANCHE_GH)
```

Les masques d'essence (`data/masques/`) doivent exister (voir « v3 : masques d'essence »).
Déploiement de l'app : `rsconnect::deployApp("app")`.

## Clôture annuelle

Au 1er janvier, les cartes repartent de zéro (chaque année est calculée de janvier à décembre,
sans lien avec la précédente). L'**année écoulée N-1 est ajoutée automatiquement à l'historique le
8 janvier** (06h17 UTC) par `cloture-annuelle.yml`, qui lance le workflow principal avec
`historique` coché et `annee_debut` = `annee_fin` = N-1. Les autres années de l'historique sont
conservées, et l'app règle ses curseurs d'après les années présentes dans la table (elle garde les
tables en mémoire jusqu'au redémarrage de son processus). En cas d'échec (SAFRAN en retard, délai
dépassé), un ticket GitHub est ouvert.

À la main (ou pour reprendre après un échec) : onglet Actions → « Mise à jour des modèles
(GeoJSON) » → *Run workflow*, case `historique` cochée, `annee_debut` = `annee_fin` = N-1, une fois
que SAFRAN dispose du 31 décembre (vers le 5 janvier : décalage de 5 jours). Le workflow
« Cloture annuelle de l'historique » accepte aussi un lancement manuel, avec l'année en option.

Autour du changement d'année, le run hebdomadaire reste sur l'année qui se termine tant que J-5
tombe dans cette année : le run du lundi 4 janvier 2027 traite encore 2026 (données jusqu'au
30 décembre) ; les cartes passent à 2027 dès que J-5 est en 2027 (run du 11 janvier). Un retard
de SAFRAN de plus de 5 jours (rare) peut faire échouer le contrôle qualité : relancer le workflow
quelques jours plus tard.

## Points de vigilance

- **Monochamus, bivoltisme A** : `mono_stade_*.geojson` porte `risque_bivoltA` et `dj10_projete`,
  extraits **au centroïde des polygones de stade**. Ces polygones regroupent les mailles par classe
  de stade : on obtient une valeur par polygone (4 polygones, donc 4 valeurs), pas une valeur par
  maille. Ces champs ne sont pas exploitables maille par maille. L'app nomme la couche « Stade +
  bivoltisme A » mais n'en affiche que le stade (légende et survol sur le champ `label`). La
  correction (extraction cellule par cellule) est à faire. Le bivoltisme B
  (`mono_bivoltB_gen2_*`) n'est pas concerné.
- **Longueur du jour** : `calculer_daylength()` (`R/safran_common.R`) appelle
  `barrks::create_daylength_rst()` (`geosphere::daylength`, avec correction de réfraction), la même
  fonction que barrks pour la diapause (seuil de 13,6 h pour CHAPY ; pour PHENIPS-Clim en scénario
  `max`, diapause thermique fonction de la durée du jour et de Tmax). `geosphere` est donc requis (il
  est installé par le workflow). La formule maison utilisée avant était 7 à 14 min plus courte : le
  seuil de CHAPY était franchi 4 à 5 jours trop tôt (pour PHENIPS-Clim, l'écart équivaut à 0,2-0,4 °C sur Tmax). Les
  tables de l'historique calculées avec l'ancienne formule sont à recalculer pour PHENIPS et CHAPY
  (STENO n'est pas concerné : `DL_INIT_MIN = 0`), sinon l'année en cours n'est pas comparable.
- **STENO** : les GeoJSON portent un champ `couleur` pré-calculé, utilisé tel quel par l'app.
- **Premier run plus lent** : le cache SAFRAN et le cache de packages R se remplissent une
  première fois.
- **Redéploiement de l'app** : à refaire après toute modification de `app/app.R`. Ne passer
  `COMPAT_ANCIENS_NOMS` à `FALSE` qu'ensuite (voir PHENIPS).

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
  PHENIPS-Clim lit des taux de developpement tabules selon la temperature moyenne de
  l'ecorce et l'amplitude journaliere (`params("phenips-clim")$dev_rates`), nuls au-dela
  d'environ 30 a 36 C selon l'amplitude. Une ecorce exposee au sud atteint souvent ces
  valeurs en ete, si bien que `shaded` peut donner plus de generations que `sunny` : c'est
  le cas en 2026 (environ 3,5 contre 2,6 generations, moyenne ponderee par la surface,
  donnees du 18/11/2026). Sur 1960-2025, `sunny` donne autant ou plus dans 94 % des
  mailles-annees. (Le seuil de 38.9 C est celui de PHENIPS, pas de PHENIPS-Clim.)
- `COMPAT_ANCIENS_NOMS <- TRUE` (phenips_functions.R) recopie aussi le scenario
  `sunny` vers les anciens noms `phenips_generations_latest.geojson` et
  `phenips_onset_latest.geojson`, le temps de redeployer l'app. Apres le redeploiement de
  l'app (`app/app.R` ne les lit plus), passer a `FALSE`, puis supprimer ces deux fichiers.
- Les anciens fichiers dates `phenips_*_AAAAMMJJ.geojson` sont purges : on garde
  les 4 dates les plus recentes par type de fichier (meme regle pour tous les modeles,
  voir « Fonctionnement »). Un `*_latest` est toujours la copie du fichier ecrit pendant
  le run (jamais d'un fichier plus ancien) ; si la couche n'a pas pu etre produite, le
  `latest` est supprime.
- **Que compte une « generation » ?** `get_generations_rst(threshold = 0)` compte les generations
  *presentes* (documentation barrks : « how many generations are present »), c'est-a-dire dont le
  developpement a commence ; elles ne sont pas forcement achevees. Les libelles sont donc `Gen. 3`,
  `Gen. 3 + couvee soeur`, sans « achevee ». L'avancement reel est donne par l'indice continu
  (somme des developpements des generations filiales, plafonnes a 1).
- Les cartes montrent un **potentiel thermique de developpement**, calcule sur les mailles qui
  croisent de l'epicea (masque v3, voir plus bas).

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
- Modele : PHENIPS-Clim, implementation du package R **barrks** (Jentschke 2025). PHENIPS-Clim
  n'est pas encore publie et a ete parametre pour le sud de l'Allemagne (documentation de barrks) ;
  les equations de temperature d'ecorce viennent de Baier et al. 2007 (PHENIPS).
- Licence du depot : **Licence Ouverte 2.0 (Etalab)**, concedant : Departement de la sante
  des forets (DSF) ; voir le fichier `LICENSE`. La mention « Source : Meteo-France » doit
  etre conservee en cas de reutilisation.

### STENO-DJ : trois scenarios (regles de Pineau et al. 2016)

`R/steno_functions.R` et `R/pipeline_steno.R` appliquent les regles de Pineau et al.
(2016) : temperature de l'air, 506 DJ au-dessus de 10.9 C, debut d'activite quand
Tmax > 18 C deux jours de suite, ni maturation ni arret photoperiodique. Les
generations sont des **maxima thermiques**, pas une prediction validee.
Trois scenarios reprennent les marges d'erreur de la these de Pineau :

| Scenario | Seuil | Constante | Fichiers `latest` |
|---|---|---|---|
| `moyen` (reference) | 10.9 C | 506 DJ | `steno_generations_latest.geojson`, `steno_essaimage_latest.geojson` |
| `pire` (rapide) | 10.7 C | 499 DJ | `steno_generations_pire_latest.geojson` |
| `moins_pire` (lent) | 11.1 C | 513 DJ | `steno_generations_moins_pire_latest.geojson` |

L'app propose ces trois scenarios dans un selecteur (« marges d'erreur de la these »).

- **Constante de 506 DJ** : valeur de la these. Le resume de l'article (Agric. Forest Entomol.
  19:23-33) donne 517 DJ avec le meme seuil de 10.9 C ; les deux valeurs sont publiees, le depot
  retient la these (506).
- **Essaimage** : le critere de vol (Tmax > 18 C deux jours de suite) ne depend pas des seuils de
  developpement ; les trois scenarios donnent donc les memes dates (`doy_*` identiques dans
  l'historique). La couche « Essaimage » est unique (`steno_essaimage_latest.geojson`,
  et `steno_essaimage_AAAA.geojson` par annee), commune aux trois scenarios.

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
  « Mise à jour des modèles (GeoJSON) » → *Run workflow* avec `construire_masques` coché. Le workflow
  commite `data/masques/`. Les pipelines s'arrêtent avec un message clair si un masque manque.
- Monochamus n'est pas masqué.

## Historique (1960 → année écoulée) : générations et premier envol par maille

Table annuelle, pour chaque maille SAFRAN (8 km) masquée par l'essence-hôte (cf. v3), du
**nombre de générations** et de la **date du premier envol**, de 1960 à la dernière année
clôturée (voir « Clôture annuelle »).

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
  cours. L'année en cours n'est pas incluse (elle est incomplète) : elle est ajoutée à la
  clôture annuelle.
- Source météo : SAFRAN (Météo-France) via GeoSAS INRAE ; la série GeoSAS démarre en août 1958,
  d'où le début à 1960 (années civiles complètes).

### Reconstruire la table

Onglet Actions → « Mise à jour des modèles (GeoJSON) » → *Run workflow* avec `historique` coché
(`annee_debut`, `annee_fin`, `parallele` réglables). Un job par année (5 en parallèle par défaut,
pour ne pas surcharger GeoSAS), puis un job d'assemblage qui commite `data/historique/`.
Les années recalculées remplacent les anciennes lignes ; les autres sont conservées.
Le calcul complet 1960-2025 prend environ 3 h 30 en théorie (≈ 16 min par année × 66 années ÷ 5 en
parallèle) ; la première exécution est plus longue car elle installe les packages R, ensuite mis en
cache. Ajouter une seule année prend le temps d'un job (≈ 16 min). Si un job reste
bloqué (déjà vu une fois à l'installation de R), « Re-run failed jobs » relance uniquement l'année
en cause : les autres résultats sont conservés. Le job
`historique_diag` (case `historique_diag`) sonde GeoSAS (étendue temporelle, durée des requêtes)
sans rien calculer.

Pour l'année en cours, `R/valider_historique.R` compare l'extraction par maille aux GeoJSON de
production (annotations `VALIDATION` du job).

### Onglet « Historique » de l'app

Carte d'une année (curseur animable), moyenne sur une période, ou tendance par décennie
(régression linéaire maille par maille, grisage des tendances non significatives à p > 0,05) ;
clic sur une maille pour afficher sa série temporelle face à la moyenne de toutes les mailles.
L'app ne fait que lire les tables du dépôt (aucun calcul lourd).

## PHENTHAU : processionnaire du chêne (*Thaumetopoea processionea*)

Cinquième modèle, package CRAN **PHENTHAUproc** (Halbig et al. 2024, *Forest Ecol. Manag.* 552:121525),
mode régional `dailymeanminmax` : températures **journalières** SAFRAN (`T_Q`, `TINF_H_Q`,
`TSUP_H_Q`), sans reconstruction horaire.

- **Calcul** : `R/pipeline_phenthau.R` (fonctions dans `R/phenthau_functions.R`), à la suite des
  autres modèles dans `modeles-update.yml` (étapes 7 à 9, `continue-on-error` : un échec
  PHENTHAU ne bloque pas les autres modèles, déjà publiés, et ouvre le ticket d'échec).
- **Données amont** : le package compte les jours de froid depuis le **1er septembre N-1**.
  `R/fetch_safran_amont.R` récupère septembre-décembre N-1 (même contrôle qualité que le fetch
  partagé, cache `safran_v2_*_<N-1>.csv`, sortie `data/cache_meteo/meteo_amont.rds`) ; le fetch
  partagé (1er janvier N) n'est pas modifié.
- **Sorties** (`data/`) : `phenthau_{stades,l3,eclosion,debourrement,mortalite}_<date>.geojson`
  (+ `_latest`), `phenthau_meta.json`. GeoTIFF L93 (22 bandes par maille : jours de l'année des stades et mortalités) et PNG
  nationaux : artefact `phenthau-rasters` du run (30 jours), non commités. Champs GeoJSON :
  `espece`, `modele = "PHENTHAU"`, `couleur`, `label`, `date_iso`, `avertissement`.
- **Éclosion** : `DEF_ECLOSION = "first"` (la plus précoce des 3 modèles : Custers, Wagenhoff,
  Meurisse). L'option `mean` du package n'a pas d'effet utile (testée sur 2 zones × 3 ans :
  mêmes dates de stades).
- **Mortalité des L1** : publiée seulement une fois le débourrement atteint (avant, le package
  renvoie 0). Le package ne donne qu'**une** valeur, celle de l'éclosion la plus précoce
  (`first`). On la recalcule avec **chacun des 3 modèles d'éclosion** : la couche `mortalite`
  porte donc `mort_custers_pct`, `mort_wagenhoff_pct`, `mort_meurisse_pct`, la **fourchette**
  `mort_min_pct` – `mort_max_pct` (texte `fourchette`) et la source (`source_mortalite`). La
  mortalité `first` (classes de couleur) est la **borne haute** de la fourchette (éclosion plus
  précoce = écart plus long au débourrement = mortalité plus forte) ; le calcul vérifie cette
  égalité et n'ajoute pas la fourchette si elle échoue. Les mortalités du package sont des
  paliers (0, 5, 10, 25, 50, 75, 100 %) déduits de la somme de températures entre l'éclosion et
  le débourrement.
  Sources de la fourchette :
  - mortalité par famine et chaîne de calcul : Halbig P., Stelzer A.S., Baier P., Pennerstorfer J.,
    Delb H., Schopf A. (2024) PHENTHAUproc, *Forest Ecol. Manag.* 552:121525,
    doi:10.1016/j.foreco.2023.121525 (voir aussi Wagenhoff et al. 2013, *J. Pest Sci.* 86:193-202,
    doi:10.1007/s10340-012-0457-7, cité par le package pour la mortalité) ;
  - éclosion, modèle Custers : Custers C. (2003) *Climate change and trophic synchronisation*,
    mémoire de master, Wageningen UR ;
  - éclosion, modèle Wagenhoff : Wagenhoff E., Wagenhoff A., Blum R., Veit H., Zapf D., Delb H.
    (2014) *Eur. J. Entomol.* 111:207-215, doi:10.14411/eje.2014.030 ;
  - éclosion, modèle Meurisse : Meurisse N., Hoch G., Schopf A., Battisti A., Grégoire J.-C.
    (2012) *Agric. Forest Entomol.* 14:239-250, doi:10.1111/j.1461-9563.2011.00562.x ;
  - débourrement : *Quercus robur*, paramétrisation `quercus_robur_clone256_type1` (défaut du
    package).
  Les trois modèles sont calibrés en Allemagne (≤ 800 m) : la fourchette mesure le désaccord
  entre modèles, **pas** une incertitude statistique. Détail technique : `PHENTHAUproc::mortality()`
  échoue avec `custers` ou `wagenhoff` seuls (vecteurs temps différents) ; le pipeline appelle donc
  la fonction interne `calc_mortality()` avec les sorties de `phenthau()`.
- **Limite d'applicabilité** : calibré sur *Quercus robur* en Allemagne, jusqu'à 800 m. Hors de ce
  domaine (sud de la France, montagne), c'est une extrapolation. Sur la zone test de Toulouse,
  l'éclosion (Meurisse) précède de 3 à 4 semaines le débourrement et la mortalité L1 atteint ou
  dépasse 50 % sur 56 à 59 mailles sur 64 selon le modèle d'éclosion. À l'échelle nationale (données
  du 30/09/2026), environ 37 % de la surface est dans les classes de mortalité 60-90 % et ≥ 90 %.
  À lire comme un indicateur, pas comme une prévision validée.
- **Tests** (entrées *Run workflow*) : `diagnostic_phenthau`, `test_phenthau` (2 petites emprises,
  2023-2025, `def_eclosion` first/mean), `test_amont` (fetch amont national, bilan dans
  `data/diag/`), `masques=chene` avec `construire_masques`.
