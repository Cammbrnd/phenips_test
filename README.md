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
  de barrks dans les deux cas). Le scenario defavorable n'est pas le meme partout :
  au-dela de 38.9 C le developpement PHENIPS s'arrete, si bien qu'en region
  mediterraneenne l'ombre peut donner plus de generations que le plein soleil.
- `COMPAT_ANCIENS_NOMS <- TRUE` (phenips_functions.R) recopie aussi le scenario
  `sunny` vers les anciens noms `phenips_generations_latest.geojson` et
  `phenips_onset_latest.geojson`, le temps de redeployer l'app. A passer a `FALSE`
  ensuite, puis supprimer ces deux fichiers.
- Les anciens fichiers dates `phenips_*_AAAAMMJJ.geojson` sont purges : on garde
  les 4 dates les plus recentes par type de fichier. Un `*_latest` est toujours
  la copie du fichier ecrit pendant le run (jamais d'un fichier plus ancien) ;
  si la couche n'a pas pu etre produite, le `latest` est supprime.
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
