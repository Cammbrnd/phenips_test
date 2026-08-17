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

- **STENO n'a pas de champ `couleur` pre-calcule** dans ses GeoJSON (a la
  difference de PHENIPS/CHAPY) : l'app genere une couleur automatiquement a
  partir du premier champ de classification. Si le rendu ne te convient pas,
  dis-le moi et j'ajoute un champ couleur cote export, comme pour PHENIPS.
- **Le patch bivoltisme A** dans `monochamus_functions.R` extrait les valeurs
  par centroide de polygone (pas cellule par cellule) : sur de tres grandes
  zones dissoutes, une petite perte de precision locale est possible aux
  bordures. Acceptable pour un indicateur de risque annuel, mais a garder en
  tete si tu compares avec le raster brut.
- **Premier run plus lent** que les suivants : le cache SAFRAN et le cache de
  packages R doivent se remplir une premiere fois.
