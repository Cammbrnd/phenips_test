# =============================================================================
# CHAPY — France metropolitaine complete
# v4.1 — Documents de sortie alignes sur phenips_france_safran_v3.R
# -----------------------------------------------------------------------------
# Espece cible  : Pityogenes chalcographus (Chalcographe / graveur a 6 dents)
# Modele barrks : CHAPY (Ogris et al. 2020, Ecol. Model. 430:109137)
#
# Source meteo  : API OGC EDR GeoSAS INRAE
#   URL         : https://api.geosas.fr/edr/collections/safran-isba/cube
#   Variables   : T_Q (Tmoy), TINF_H_Q (Tmin), TSUP_H_Q (Tmax)
#   Resolution  : grille SAFRAN 8x8 km, France metropolitaine, 1958-2026
#   Avantage    : requete bbox L93 directe, pas de point-par-point
#
# Mention obligatoire : "Source : Meteo France"
# Ref. SAFRAN : Quintana-Segui et al. (2008) J. Appl. Meteor. Climatol.
#
# =============================================================================
# NOUVEAUTES v4 — PARITE DES SORTIES AVEC PHENIPS
# =============================================================================
# Objectif : produire pour CHAPY exactement la meme famille de documents que
# phenips_france_safran_v3.R, avec le meme nommage, afin que les deux modeles
# soient superposables dans QGIS et dans un visualiseur web.
#
#   PHENIPS                              CHAPY v4
#   -----------------------------------  ---------------------------------------
#   phenips_onset_{annee}_L93.tif        chapy_onset_{annee}_L93.tif
#   phenips_btmean_{date}_L93.tif        chapy_tmax_{date}_L93.tif          (*)
#   phenips_gen1_dev_{date}_L93.tif      chapy_gen1_dev_{date}_L93.tif
#   phenips_generations_{date}_L93.tif   chapy_generations_{date}_L93.tif
#   phenips_generations_{date}.geojson   chapy_generations_{date}.geojson
#   phenips_onset_{date}.geojson         chapy_onset_{date}.geojson
#   phenips_onset_{date}.png             chapy_onset_{date}.png
#   phenips_btmean_{date}.png            chapy_tmax_{date}.png             (*)
#   phenips_gen1_{date}.png              chapy_gen1_{date}.png
#   phenips_generations_{date}.png       chapy_generations_{date}.png
#   (absent de PHENIPS)                  chapy_premier_envol_{date}_doy_L93.tif
#
#   (*) CHAPY ne modelise PAS de temperature sous ecorce : le modele ne prend
#       pas le rayonnement en entree et travaille sur la temperature de l'air.
#       L'equivalent fonctionnel de btmean est donc Tmax, variable qui pilote
#       le seuil de vol printanier. Fichier nomme tmax et non btmean pour
#       eviter toute confusion d'interpretation biologique.
#
# Sorties propres a CHAPY, conservees (sans equivalent PHENIPS) :
#   chapy_france_{date}.png        cartes mensuelles de generations
#   chapy_france_onset.png         courbe de progression du vol printanier
#   chapy_france_tmax.png          courbe Tmax / Tmean nationale
#   chapy_france_photoperiode.png  courbe photoperiode et seuil de diapause
#   chapy_france_gen1.png          courbe de developpement de la generation 1
#   diagramme_*.png                diagrammes de developpement par station
#
# -----------------------------------------------------------------------------
# ECARTS ASSUMES vs phenips_v3 (corrections deliberees, pas des oublis)
#
#   [A] Carte onset : PHENIPS cartographie l'INDICE de couche temporelle en
#       l'etiquetant "jour de l'annee". Exact seulement si la serie demarre au
#       01/01 ET sans jour manquant. Ici : vrai DOY via format(date, "%j").
#
#   [B] max.col() : PHENIPS l'applique a une matrice pouvant contenir des NA
#       (propagation NA -> indice faux, maille perdue). Ici : NA forces a FALSE.
#
#   [C] Vectorisation : dissolve = TRUE (comme PHENIPS) au lieu du
#       dissolve = FALSE de la v3. Sur grille 8 km France entiere :
#       ~14 000 polygones -> quelques dizaines. Fallback automatique sur
#       dissolve = FALSE si terra refuse.
#
#   [D] Attributs GeoJSON : ajout de "espece" et "modele" (absents de PHENIPS),
#       necessaires des lors que les couches PHENIPS et CHAPY cohabitent dans
#       le meme visualiseur. Recommande d'ajouter les memes champs cote PHENIPS.
#
#   [E] Generations : normalisation des deux conventions barrks (numerique
#       1 / 1.5 et categorielle "1" / "1s") via normaliser_generation().
#
#   [F] catg() : glue() applique .trim = TRUE et supprime les lignes blanches
#       de TETE et de QUEUE du template. Les cat(glue("...\n")) de la v3 et de
#       PHENIPS n'emettent donc aucun retour a la ligne, d'ou les lignes de
#       console collees les unes aux autres. catg() le garantit.
#       ATTENTION : catg() est defini en section 0, avant tout code de niveau
#       global qui l'utilise (section 1). Le deplacer plus bas casse le
#       chargement du script.
#
#   [G] premier envol : terra::app() (v3) applique une fonction R cellule par
#       cellule (~14 000 appels). Remplace par max.col() sur la matrice
#       entiere, en une passe.
#
#   [H] catg() doit recevoir .envir = parent.frame() : glue() interpole sinon
#       dans le frame de catg(), rendant invisible toute variable LOCALE a la
#       fonction appelante. Bug reel remonte au premier run reel (v4.0) :
#       "objet 'nom' introuvable" dans exporter_rasters_chapy().
#
#   [I] create_daylength_rst() refuse un template sans valeurs. RST_FRANCE
#       etant une grille vide, la correction v3 [1] etait morte : le fallback
#       manuel se declenchait a chaque run ("[*] raster has no values").
#       Le template est desormais initialise avant l'appel.
#
#   [J] options(datatable.showProgress = FALSE) : coupe les messages
#       "Traitement de N groupes..." qui se melaient aux traces du pipeline.
#
# Fonctions v3 devenues obsoletes et supprimees :
#   calculer_premier_envol(), exporter_geojson_depuis_raster(),
#   exporter_geojson_chapy(), raster_est_vide() (alias conserve vers rst_vide())
#
# -----------------------------------------------------------------------------
# HISTORIQUE v3 (corrections conservees telles quelles)
#   [1] create_daylength_rst() barrks natif au lieu d'une formule maison
#   [2] gen_vide : suppression du || is.numeric(vals) qui bloquait les cartes
#   [3] synthese_console() via get_onset_rst() + comptage v > 0
#   [4] assembler_rasters_nationaux() : pivot large data.table (matriciel)
#   [5] dates extraites explicitement depuis dt avant rm(dt)
#   v3.1 isTRUE() sur vecteur supprime (retournait toujours FALSE)
#   v3.2 mat_large[[i]] extrayait une colonne, pas une ligne -> as.matrix()
#   v3.4 signature create_daylength_rst(template=, dates=) confirmee
#
# Frequence de mise a jour : 1 fois/mois apres le 5 du mois
#   (GeoSAS met a jour SAFRAN debut de chaque mois avec le mois precedent)
#
# Auteur : CNPF / CRPF AuRA
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

# install.packages(c("barrks", "terra", "httr2", "ggplot2",
#                    "tidyterra", "lubridate", "glue", "sf",
#                    "scales", "data.table"))

library(barrks)
library(terra)
library(httr2)
library(ggplot2)
library(tidyterra)
library(lubridate)
library(glue)
library(sf)
library(scales)
library(data.table)


# --- Utilitaire d'affichage console -----------------------------------------
# Defini ICI, en section 0, et non aupres des autres helpers : catg() est
# utilise des la section 1 par du code de NIVEAU GLOBAL, execute au moment du
# source(). Une definition plus bas dans le fichier produit
# "impossible de trouver la fonction catg". Les fonctions, elles, resolvent
# leurs appels a l'execution : leur ordre de definition est indifferent.

#' cat + glue avec mise en forme litterale et retour a la ligne garanti.
#'
#' [F] glue() applique .trim = TRUE par defaut : suivant les regles des
#' docstrings Python, les lignes blanches de TETE et de QUEUE du template sont
#' supprimees. cat(glue("texte\n")) n'emet donc aucun retour a la ligne, et
#' cat(glue("\n>> titre")) perd sa ligne vide de separation — d'ou l'affichage
#' compact des consoles v3 et PHENIPS. .trim = FALSE preserve le template.
#'
#' [H] .envir est OBLIGATOIRE : glue() interpole par defaut dans
#' parent.frame(), qui vaut ici le frame de catg() et non celui de l'appelant.
#' Sans ce parametre, toute variable LOCALE a la fonction appelante
#' ("objet 'nom' introuvable") est invisible ; seules les variables globales
#' passent, par heritage lexical — d'ou un bug qui ne se voit pas en test
#' avec des variables de haut niveau. Le defaut parent.frame() est evalue
#' paresseusement dans le frame de catg(), il designe donc bien l'appelant.
catg <- function(..., .envir = parent.frame()) {
  cat(glue(..., .trim = FALSE, .envir = .envir), "\n", sep = "")
}


# =============================================================================
# 1. PARAMETRES GLOBAUX
# =============================================================================

CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))

# SAFRAN : GeoSAS met a jour debut de chaque mois avec le mois precedent.
# On prend le dernier jour du mois precedent comme date fin sure.
# Ajustable manuellement : DATE_FIN <- as.Date("2026-03-31")
DATE_FIN <- floor_date(Sys.Date(), "month") - 1

# --- Grille SAFRAN France metropolitaine (L93) ---
# IMPORTANT : coordonnees en entiers via formatC() -> evite "6e+06" dans URL
FRANCE_XMIN <-  100000L
FRANCE_XMAX <- 1200000L
FRANCE_YMIN <- 6050000L
FRANCE_YMAX <- 7150000L
SAFRAN_RES  <-    8000L

# Decoupage en dalles pour eviter les timeouts API (200 km = ~25 tuiles/cote)
DALLE_KM <- 200L
DALLE_M  <- DALLE_KM * 1000L

# API GeoSAS
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
# CHAPY : T_Q (tmoy), TINF_H_Q (tmin), TSUP_H_Q (tmax) — pas de rayonnement
GEOSAS_VARS <- "T_Q,TINF_H_Q,TSUP_H_Q"

# Seuils biologiques du Chalcographe
SEUIL_VOL_TMAX <- 15.6   # Tmax >= 15.6 C -> vol printanier
SEUIL_DIAPAUSE <- 13.6   # Photopériode < 13.6 h -> diapause

# v4 [J] : coupe les messages "Traitement de N groupes..." emis par data.table
# lors des agregations GForce de assembler_rasters_nationaux(), qui se melangent
# aux traces du pipeline et cassent la lisibilite de la console.
options(datatable.showProgress = FALSE)

DIR_OUT   <- "data"
DIR_CACHE <- file.path(DIR_OUT, "cache_safran")
dir.create(DIR_OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

cat("=== CHAPY v4.1 — Pityogenes chalcographus — France entiere ===\n")
cat("    Source meteo : API GeoSAS INRAE / SAFRAN\n")
catg("    Periode      : {DATE_DEB} -> {DATE_FIN}")
catg("    Resolution   : {SAFRAN_RES/1000} km | dalles {DALLE_KM}x{DALLE_KM} km")
catg("    Variables    : {GEOSAS_VARS}")
catg("    Sorties      : {DIR_OUT}/")
catg("    Seuil vol    : Tmax >= {SEUIL_VOL_TMAX} C")
catg("    Diapause     : photopériode < {SEUIL_DIAPAUSE} h\n")
cat("  [INFO] Mise a jour GeoSAS : debut de chaque mois (mois precedent).\n")
cat("  [INFO] Relancer 1 fois/mois apres le 5 du mois.\n\n")


# =============================================================================
# 2. GRILLE DE DALLES ET RASTER DE REFERENCE
# =============================================================================

# =============================================================================
# 3. TELECHARGEMENT SAFRAN PAR DALLE (avec cache incremental)
# =============================================================================

# =============================================================================
# 4. LECTURE ET ASSEMBLAGE DES CSV (inspire de Monochamus)
# =============================================================================

#' Detecte une colonne par patterns regex (insensible a la casse).
#' Retourne NULL si aucun pattern ne correspond.
# =============================================================================
# 5. CONSTRUCTION DES RASTERS BARRKS (approche matricielle via pivot large)
# =============================================================================

#' Assemble les rasters nationaux depuis le data.table SAFRAN.
#' CORRECTION v3 [4] : vraie approche matricielle via dcast() + affectation
#' directe des valeurs -> evite la boucle rasterize() jour par jour (~10x plus rapide).
#'
#' Principe :
#'   1. Creer un index de cellule L93 pour chaque point SAFRAN (cellFromXY)
#'   2. dcast(date ~ cell_idx, value.var = variable) -> matrice [njours x ncell]
#'   3. Affecter la matrice transposee en une seule fois avec values(rst) <- mat
#'
#' CORRECTION v3.2 (bug de production) : la version precedente utilisait
#' mat_large[[idx_dates[j]]] dans une boucle pour extraire "la ligne j" -- mais
#' [[ ]] sur un data.table extrait une COLONNE (comme sur une liste), pas une
#' ligne. Ca renvoyait un vecteur de longueur "nombre de jours" au lieu de
#' "nombre de cellules", d'ou l'erreur "nombre d'objets a remplacer n'est pas
#' multiple de la taille du remplacement". Fix : conversion en matrice reelle
#' via as.matrix() puis extraction de LIGNES par indexation numerique standard,
#' vectorisee en un seul bloc (plus rapide qu'une boucle en plus d'etre juste).
# =============================================================================
# 6. DUREE DU JOUR (via barrks natif — correction v3 [1])
# =============================================================================

#' CORRECTION v3 [1] : remplace calculer_daylength() (formule maison) par
#' create_daylength_rst() de barrks, qui garantit la coherence exacte avec
#' les calculs internes de CHAPY (meme algorithme astronomique).
#'
#' CORRECTION v3.4 (signature CONFIRMEE par l'utilisateur en local) :
#' args(create_daylength_rst) donne :
#'   function(template, dates = terra::time(template), crs = "EPSG:4258",
#'            .quiet = FALSE)
#' Ca confirme le diagnostic v3.3 : le raster de reference est bien le
#' PREMIER argument (nomme `template`, pas `rst`), `dates` est le second.
#' Appel desormais fait par arguments nommes pour plus de robustesse et de
#' lisibilite (protege contre un futur changement d'ordre dans barrks).
#'
#' Note sur `crs` (defaut "EPSG:4258", ETRS89) : ce parametre precise le
#' systeme geographique utilise en interne pour calculer les latitudes
#' (reprojection automatique depuis le CRS de `template`, ici L93). ETRS89
#' et WGS84 donnent des latitudes quasi identiques (ecart de l'ordre du
#' centimetre) -> aucun impact mesurable sur la photoperiode calculee, pas
#' besoin de le surcharger pour la France.
#'
#' Le fallback ci-dessous est conserve par securite (defense en profondeur)
#' mais ne devrait plus jamais se declencher avec la signature confirmee.
calculer_daylength_manuel <- function(rst_ref, dates) {
  lats <- crds(rst_ref) |>
    as.data.frame() |>
    vect(geom = c("x", "y"), crs = CRS_L93) |>
    project(CRS_WGS) |>
    crds() |>
    (\(m) m[, 2])()
  
  duree <- function(lat_deg, doy) {
    decl   <- -asin(0.39779 * cos(pi/180 * (0.98563*(doy+10) +
                                              1.914*sin(pi/180 * 0.98563*(doy-2)))))
    cos_ha <- pmin(pmax(-tan(lat_deg*pi/180) * tan(decl), -1), 1)
    2 * acos(cos_ha) * 12 / pi
  }
  
  couches <- lapply(seq_along(dates), function(i) {
    r         <- rst_ref
    values(r) <- duree(lats, yday(dates[i]))
    names(r)  <- as.character(dates[i])
    r
  })
  
  rst       <- rast(couches)
  time(rst) <- dates
  rst
}


# =============================================================================
# 7. MODELISATION CHAPY
# =============================================================================

modeliser_chapy <- function(meteo, mode = "max") {
  catg("\n>> Modelisation CHAPY [mode={mode}]...")
  cat("   Espece : Pityogenes chalcographus (Chalcographe)\n")
  cat("   Ref.   : Ogris et al. (2020) Ecol. Model. 430:109137\n")
  
  pheno <- phenology(
    "chapy",
    tmin      = meteo$tmin,
    tmean     = meteo$tmean,
    tmax      = meteo$tmax,
    daylength = meteo$daylength,
    mode      = mode
  )
  cat("   Modelisation terminee.\n")
  pheno
}


# =============================================================================
# 8. HELPERS D'ACCES BARRKS ET NOMENCLATURE DES GENERATIONS
# =============================================================================

# Note : catg() est defini en section 0 (contrainte d'ordre, voir commentaire).

#' TRUE si le raster est NULL, vide ou integralement NA.
rst_vide <- function(rst) {
  if (is.null(rst)) return(TRUE)
  tryCatch({
    v <- values(rst)
    length(v) == 0 || all(is.na(v))
  }, error = function(e) TRUE)
}

#' Alias de compatibilite : conserve pour les fonctions v3 non modifiees
#' (cartographier_generations() notamment).
raster_est_vide <- rst_vide

#' Vecteur des dates de la serie modelisee.
#' Equivalent de pheno$dates cote PHENIPS, que barrks ne construit pas pour
#' CHAPY. Source de verite : time(pheno$onset).
dates_chapy <- function(pheno) {
  d <- tryCatch(as.Date(time(pheno$onset)), error = function(e) NULL)
  if (is.null(d) || all(is.na(d))) {
    stop("time(pheno$onset) indisponible — verifier assembler_rasters_nationaux().")
  }
  d
}

#' Derniere date effectivement couverte par le modele.
#' Peut differer de DATE_FIN si SAFRAN n'a pas publie le mois complet.
date_fin_chapy <- function(pheno) max(dates_chapy(pheno))

#' Raster de developpement de la generation 1 (NULL si non demarree).
get_gen1_chapy <- function(pheno) {
  if (!is.null(pheno[["development"]][["gen_1"]]))        return(pheno$development$gen_1)
  if (!is.null(pheno[["development"]][["generation_1"]])) return(pheno$development$generation_1)
  NULL
}

#' Raster categoriel des generations a une date donnee (NULL si aucune complete).
get_gen_rst_chapy <- function(pheno, date_ref = NULL) {
  if (is.null(date_ref)) date_ref <- date_fin_chapy(pheno)
  tryCatch(
    get_generations_rst(pheno, as.character(date_ref)),
    error = function(e) {
      message(glue("  [INFO] get_generations_rst({date_ref}) : {e$message}"))
      message("  -> Aucune generation complete (normal avant l'ete).")
      NULL
    }
  )
}

# --- Nomenclature des generations -------------------------------------------
# [E] barrks stocke les generations soit en numerique (1, 1.5, 2...), soit en
# categoriel ("1", "1s", "2"...) selon la version et le modele. Les deux
# conventions sont normalisees vers une valeur numerique unique.

#' "1" / 1 -> 1.0 ; "1s" / 1.5 -> 1.5 (couvee soeur / sister brood)
normaliser_generation <- function(x) {
  x_chr <- as.character(x)
  soeur <- grepl("s$", x_chr, ignore.case = TRUE)
  base  <- suppressWarnings(as.numeric(sub("s$", "", x_chr, ignore.case = TRUE)))
  ifelse(soeur & !is.na(base), base + 0.5, base)
}

LABELS_GENERATIONS <- c(
  "1"   = "Generation 1 complete",
  "1.5" = "Couvee soeur (gen 1.5)",
  "2"   = "Generation 2 complete",
  "2.5" = "Couvee soeur (gen 2.5)",
  "3"   = "Generation 3 complete",
  "3.5" = "Couvee soeur (gen 3.5)",
  "4"   = "Generation 4 complete"
)

# Palette alignee sur barrks_colors() / viridis, identique a phenips_v3
PALETTE_GENERATIONS <- c(
  "1"   = "#35B779FF",
  "1.5" = "#FDE725FF",
  "2"   = "#31688EFF",
  "2.5" = "#FFA500FF",
  "3"   = "#440154FF",
  "3.5" = "#C0392BFF",
  "4"   = "#7B241CFF"
)

#' Libelle lisible pour popup web.
libelle_generation <- function(gen_num) {
  cle <- as.character(gen_num)
  lab <- unname(LABELS_GENERATIONS[match(cle, names(LABELS_GENERATIONS))])
  lab[is.na(lab)] <- paste0("Generation ", cle[is.na(lab)])
  lab
}

#' Couleur hex pour rendu MapLibre / Leaflet.
couleur_generation <- function(gen_num) {
  cle <- as.character(gen_num)
  col <- unname(PALETTE_GENERATIONS[match(cle, names(PALETTE_GENERATIONS))])
  col[is.na(col)] <- "#CCCCCC"
  col
}

#' Palette nommee pour les niveaux d'un raster categoriel.
#' Priorite : barrks_colors() natif -> palette interne -> gris neutre.
palette_pour_niveaux <- function(niveaux) {
  niveaux <- as.character(niveaux)
  cols    <- setNames(rep(NA_character_, length(niveaux)), niveaux)
  
  bc <- tryCatch(barrks_colors(), error = function(e) NULL)
  if (!is.null(bc) && !is.null(names(bc))) {
    m  <- match(niveaux, names(bc))
    ok <- !is.na(m)
    cols[ok] <- unname(bc[m[ok]])
  }
  
  idx_reste <- which(is.na(cols))
  if (length(idx_reste) > 0) {
    cle <- as.character(normaliser_generation(niveaux[idx_reste]))
    m2  <- match(cle, names(PALETTE_GENERATIONS))
    ok2 <- !is.na(m2)
    cols[idx_reste[ok2]] <- unname(PALETTE_GENERATIONS[m2[ok2]])
  }
  
  cols[is.na(cols)] <- "#CCCCCC"
  cols
}


# =============================================================================
# 9. DATE DE PREMIER ENVOL PAR MAILLE
# =============================================================================
# [G] Remplace calculer_premier_envol() de la v3, qui passait par terra::app()
#     et appelait une fonction R pour chacune des ~14 000 cellules. max.col()
#     traite la matrice entiere en une passe.
# [B] Les NA sont forces a FALSE avant max.col() : sinon une maille presentant
#     un trou SAFRAN avant son envol remonte un indice arbitraire ou NA et
#     disparait de la carte.

#' Indice (1..nlyr) du premier jour ou onset > 0 ; NA si aucun envol.
indice_premier_envol <- function(onset_stack) {
  vals     <- values(onset_stack)                 # matrice ncell x nlyr
  a_envole <- rowSums(vals > 0, na.rm = TRUE) > 0
  
  idx <- rep(NA_integer_, nrow(vals))
  if (any(a_envole, na.rm = TRUE)) {
    sous_bin <- vals[a_envole, , drop = FALSE] > 0
    sous_bin[is.na(sous_bin)] <- FALSE            # [B]
    idx[a_envole] <- max.col(sous_bin, ties.method = "first")
  }
  idx
}

#' Raster du jour de l'annee (DOY) du premier envol.
#' Template pris sur pheno$onset et non sur RST_FRANCE : garantit une geometrie
#' identique meme si l'emprise du modele a ete recadree en amont.
raster_premier_envol_doy <- function(pheno) {
  onset <- pheno$onset
  dates <- dates_chapy(pheno)
  
  idx <- indice_premier_envol(onset)
  ok  <- !is.na(idx)
  
  doy     <- rep(NA_real_, length(idx))
  doy[ok] <- as.integer(format(dates[idx[ok]], "%j"))   # [A] vrai DOY
  
  rst <- rast(onset[[1]])
  values(rst) <- doy
  names(rst)  <- "doy_premier_envol"
  
  list(rst = rst, idx = idx, dates = dates, n_envol = sum(ok))
}


# =============================================================================
# 10. EXPORT GEOTIFF (parite PHENIPS : LZW + TILED + BIGTIFF)
# =============================================================================

exporter_rasters_chapy <- function(pheno, meteo = NULL, dir_out = DIR_OUT) {
  catg(">> Export GeoTIFF (L93, LZW)...")
  
  date_ref <- date_fin_chapy(pheno)
  date_j   <- format(date_ref, "%Y%m%d")
  
  ecrire <- function(rst, nom) {
    ch <- file.path(dir_out, nom)
    writeRaster(rst, ch, overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    catg("   {nom}")
  }
  
  # --- 1. Pile onset complete de l'annee ---
  # La v3 n'exportait que la derniere couche : rejouer une date anterieure
  # imposait de relancer tout le pipeline.
  ecrire(pheno$onset, glue("chapy_onset_{ANNEE}_L93.tif"))
  
  # --- 2. DOY du premier envol (exploitable directement dans QGIS) ---
  envol <- raster_premier_envol_doy(pheno)
  if (!rst_vide(envol$rst)) {
    ecrire(envol$rst, glue("chapy_premier_envol_{date_j}_doy_L93.tif"))
    catg("   -> {envol$n_envol} mailles avec envol detecte")
  } else {
    catg("   [INFO] Aucun envol detecte sur la periode.")
  }
  
  # --- 3. Variable thermique pilote : Tmax (equivalent fonctionnel de btmean) ---
  if (!is.null(meteo) && !is.null(meteo$tmax)) {
    rst_t <- meteo$tmax[[nlyr(meteo$tmax)]]
    names(rst_t) <- glue("tmax_{date_j}")
    ecrire(rst_t, glue("chapy_tmax_{date_j}_L93.tif"))
  } else {
    catg("   [INFO] meteo absent -> export tmax ignore.")
  }
  
  # --- 4. Developpement de la generation 1, en pourcentage ---
  gen1 <- get_gen1_chapy(pheno)
  if (!is.null(gen1)) {
    rst_g1 <- gen1[[nlyr(gen1)]] * 100
    rst_g1[rst_g1 <= 0] <- NA          # masque l'absence de developpement
    names(rst_g1) <- glue("gen1_dev_pct_{date_j}")
    ecrire(rst_g1, glue("chapy_gen1_dev_{date_j}_L93.tif"))
  } else {
    catg("   [INFO] gen1 non disponible (normal avant mi-saison).")
  }
  
  # --- 5. Generations completes ---
  rst_gen <- get_gen_rst_chapy(pheno, date_ref)
  if (!rst_vide(rst_gen)) {
    names(rst_gen) <- glue("generations_{date_j}")
    ecrire(rst_gen, glue("chapy_generations_{date_j}_L93.tif"))
  }
  
  catg("\n  Rasters exportes dans {dir_out}/\n")
  invisible(TRUE)
}


# =============================================================================
# 11. EXPORT GEOJSON (WGS84, dissous, pret pour MapLibre / Leaflet)
# =============================================================================
# [C] dissolve = TRUE avec repli automatique sur dissolve = FALSE.
#     RFC 7946 impose WGS84 pour le GeoJSON -> reprojection systematique,
#     method = "near" pour preserver les valeurs discretes.

#' Vectorise un raster en polygones dissous, avec repli si terra refuse.
vectoriser_dissous <- function(rst_wgs, contexte = "") {
  tryCatch(
    st_as_sf(as.polygons(rst_wgs, dissolve = TRUE)),
    error = function(e) {
      message("  dissolve=TRUE refuse ", contexte, " (", e$message,
              ") -> repli sur dissolve=FALSE")
      tryCatch(
        st_as_sf(as.polygons(rst_wgs, dissolve = FALSE, na.rm = TRUE)),
        error = function(e2) NULL
      )
    }
  )
}

produire_geojson_generations_chapy <- function(pheno, dir_out = DIR_OUT) {
  catg(">> Vectorisation generations -> GeoJSON...")
  
  date_ref <- date_fin_chapy(pheno)
  date_j   <- format(date_ref, "%Y%m%d")
  
  rst_gen <- get_gen_rst_chapy(pheno, date_ref)
  if (rst_vide(rst_gen)) {
    catg("  Pas de generation complete disponible (normal avant juillet).\n")
    return(invisible(NULL))
  }
  
  rst_wgs  <- project(rst_gen, CRS_WGS, method = "near")
  vect_gen <- vectoriser_dissous(rst_wgs, "(generations)")
  
  if (is.null(vect_gen) || nrow(vect_gen) == 0) {
    catg("  Aucun polygone produit.\n")
    return(invisible(NULL))
  }
  
  # [E] normalisation des deux conventions barrks
  names(vect_gen)[1]  <- "generation_brute"
  vect_gen$generation <- normaliser_generation(vect_gen$generation_brute)
  vect_gen$generation_brute <- NULL
  
  vect_gen$label       <- libelle_generation(vect_gen$generation)
  vect_gen$couleur     <- couleur_generation(vect_gen$generation)
  vect_gen$date_calcul <- as.character(date_ref)
  vect_gen$annee       <- ANNEE
  vect_gen$espece      <- "Pityogenes chalcographus"   # [D]
  vect_gen$modele      <- "CHAPY"                      # [D]
  
  # Simplification : tolerance en degres, ~500 m aux latitudes francaises
  vect_gen <- st_simplify(vect_gen, dTolerance = 0.005, preserveTopology = TRUE)
  vect_gen <- vect_gen[!st_is_empty(vect_gen), ]
  
  chemin <- file.path(dir_out, glue("chapy_generations_{date_j}.geojson"))
  st_write(vect_gen, chemin, delete_dsn = TRUE, quiet = TRUE,
           layer_options = "COORDINATE_PRECISION=5")
  
  catg("  {nrow(vect_gen)} polygones | {round(file.size(chemin)/1024, 0)} Ko",
       " -> {basename(chemin)}\n")
  invisible(vect_gen)
}


produire_geojson_envol_chapy <- function(pheno, dir_out = DIR_OUT) {
  catg(">> Vectorisation premier envol -> GeoJSON...")
  
  date_ref <- date_fin_chapy(pheno)
  date_j   <- format(date_ref, "%Y%m%d")
  
  envol <- raster_premier_envol_doy(pheno)
  if (rst_vide(envol$rst)) {
    catg("  Aucun envol detecte sur la periode.\n")
    return(invisible(NULL))
  }
  
  rst_wgs    <- project(envol$rst, CRS_WGS, method = "near")
  vect_envol <- vectoriser_dissous(rst_wgs, "(premier envol)")
  
  if (is.null(vect_envol) || nrow(vect_envol) == 0) {
    catg("  Aucun polygone produit.\n")
    return(invisible(NULL))
  }
  
  names(vect_envol)[1] <- "doy"
  vect_envol$doy <- as.integer(vect_envol$doy)
  
  origine <- as.Date(glue("{ANNEE}-01-01"))
  d_env   <- origine + vect_envol$doy - 1L
  
  vect_envol$date_envol  <- format(d_env, "%d/%m/%Y")   # popup web
  vect_envol$date_iso    <- format(d_env, "%Y-%m-%d")   # tri / filtrage
  vect_envol$semaine     <- as.integer(format(d_env, "%V"))
  vect_envol$date_calcul <- as.character(date_ref)
  vect_envol$annee       <- ANNEE
  vect_envol$espece      <- "Pityogenes chalcographus"  # [D]
  vect_envol$modele      <- "CHAPY"                     # [D]
  vect_envol$seuil_vol   <- SEUIL_VOL_TMAX
  vect_envol$doy         <- NULL
  
  vect_envol <- st_simplify(vect_envol, dTolerance = 0.005, preserveTopology = TRUE)
  vect_envol <- vect_envol[!st_is_empty(vect_envol), ]
  
  chemin <- file.path(dir_out, glue("chapy_onset_{date_j}.geojson"))
  st_write(vect_envol, chemin, delete_dsn = TRUE, quiet = TRUE,
           layer_options = "COORDINATE_PRECISION=5")
  
  catg("  {nrow(vect_envol)} polygones | ",
       "{length(unique(vect_envol$date_iso))} dates uniques | ",
       "{round(file.size(chemin)/1024, 0)} Ko -> {basename(chemin)}\n")
  invisible(vect_envol)
}


# =============================================================================
# 12. CARTES NATIONALES (equivalent carte_france() de phenips_v3)
# =============================================================================

carte_france_chapy <- function(rst_l93, titre, legende, palette,
                               breaks     = NULL,
                               type       = c("gradient", "categoriel"),
                               date_label = NULL,
                               fichier    = NULL) {
  type <- match.arg(type)
  if (is.null(date_label)) date_label <- format(DATE_FIN, "%d/%m/%Y")
  
  rst_wgs <- project(rst_l93, CRS_WGS,
                     method = if (type == "categoriel") "near" else "bilinear")
  if (type == "categoriel") rst_wgs <- as.factor(rst_wgs)
  
  p <- ggplot() +
    geom_spatraster(data = rst_wgs) +
    {
      if (type == "gradient") {
        scale_fill_gradientn(
          colours  = palette,
          values   = if (!is.null(breaks)) rescale(breaks) else NULL,
          na.value = "grey92",
          name     = legende,
          guide    = guide_colorbar(barwidth = 14, barheight = 0.8,
                                    title.position = "top")
        )
      } else {
        scale_fill_manual(
          values   = palette,
          na.value = "grey92",
          name     = legende,
          drop     = FALSE
        )
      }
    } +
    fond_france() +
    coord_sf(xlim = c(-5.2, 9.7), ylim = c(41.2, 51.2)) +
    labs(
      title    = titre,
      subtitle = glue("CHAPY (Ogris et al. 2020) | SAFRAN 8 km — GeoSAS INRAE",
                      " | {date_label} | {ANNEE}"),
      caption  = "Source : Meteo France — CNPF/CRPF AuRA",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.caption    = element_text(size = 8, color = "grey50"),
      legend.position = "bottom"
    )
  
  if (!is.null(fichier)) {
    ggsave(fichier, plot = p, width = 10, height = 9, dpi = 150)
    catg("  Carte : {basename(fichier)}")
  }
  invisible(p)
}


#' Genere les 4 cartes nationales de synthese, format PHENIPS.
generer_cartes_nationales_chapy <- function(pheno, meteo = NULL, dir_out = DIR_OUT) {
  catg(">> Cartes nationales (format PHENIPS)...")
  
  date_ref   <- date_fin_chapy(pheno)
  date_j     <- format(date_ref, "%Y%m%d")
  date_label <- format(date_ref, "%d/%m/%Y")
  
  # --- 1. Date de premier envol (DOY) ---
  envol <- raster_premier_envol_doy(pheno)
  if (!rst_vide(envol$rst)) {
    carte_france_chapy(
      envol$rst,
      titre      = "Date de premier envol — Pityogenes chalcographus — France",
      legende    = glue("Jour de l'annee (1 = 01/01/{ANNEE})"),
      palette    = c("#FFF3B0", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"),
      type       = "gradient",
      date_label = date_label,
      fichier    = file.path(dir_out, glue("chapy_onset_{date_j}.png"))
    )
  } else {
    catg("  Carte envol : aucun vol printanier declenche.")
  }
  
  # --- 2. Tmax a la derniere date, avec le seuil de vol dans la legende ---
  if (!is.null(meteo) && !is.null(meteo$tmax)) {
    carte_france_chapy(
      meteo$tmax[[nlyr(meteo$tmax)]],
      titre      = glue("Temperature maximale — {date_label}"),
      legende    = glue("Tmax (C) — seuil vol {SEUIL_VOL_TMAX} C"),
      palette    = c("#E8EDF2", "#A8C6DF", "#FDE725FF", "#F39C12", "#C0392B"),
      breaks     = c(-5, 5, SEUIL_VOL_TMAX, 25, 38),
      type       = "gradient",
      date_label = date_label,
      fichier    = file.path(dir_out, glue("chapy_tmax_{date_j}.png"))
    )
  }
  
  # --- 3. Developpement de la generation 1, en pourcentage ---
  gen1 <- get_gen1_chapy(pheno)
  if (!is.null(gen1)) {
    rst_g1 <- gen1[[nlyr(gen1)]] * 100
    rst_g1[rst_g1 <= 0] <- NA
    if (!rst_vide(rst_g1)) {
      carte_france_chapy(
        rst_g1,
        titre      = glue("Developpement generation 1 — P. chalcographus — {date_label}"),
        legende    = "Developpement gen 1 (%)",
        palette    = c("#FFF7BC", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"),
        breaks     = c(0, 25, 50, 75, 100),
        type       = "gradient",
        date_label = date_label,
        fichier    = file.path(dir_out, glue("chapy_gen1_{date_j}.png"))
      )
    }
  } else {
    catg("  Carte gen1 : generation 1 non demarree.")
  }
  
  # --- 4. Generations completes (categoriel) ---
  rst_gen <- get_gen_rst_chapy(pheno, date_ref)
  if (!rst_vide(rst_gen)) {
    rst_cat <- as.factor(rst_gen)
    niveaux <- tryCatch(levels(rst_cat)[[1]][, 2], error = function(e) NULL)
    if (is.null(niveaux)) niveaux <- sort(unique(na.omit(values(rst_gen))))
    
    carte_france_chapy(
      rst_gen,
      titre      = glue("Generations P. chalcographus — France — {date_label}"),
      legende    = "Generation",
      palette    = palette_pour_niveaux(niveaux),
      type       = "categoriel",
      date_label = date_label,
      fichier    = file.path(dir_out, glue("chapy_generations_{date_j}.png"))
    )
  } else {
    catg("  Carte generations : aucune generation complete (normal avant l'ete).")
  }
  
  catg("  Cartes exportees dans {dir_out}/\n")
  invisible(TRUE)
}


# =============================================================================
# 13. VISUALISATIONS COMPLEMENTAIRES (specifiques CHAPY)
# =============================================================================
# Cartes mensuelles de generations et courbes temporelles nationales.
# Ces sorties n'ont pas d'equivalent PHENIPS et sont conservees telles quelles.


fond_france <- function() {
  tryCatch({
    if (requireNamespace("rnaturalearth", quietly = TRUE)) {
      ne <- rnaturalearth::ne_countries(country = "France",
                                        scale = "medium", returnclass = "sf")
      return(geom_sf(data = ne, fill = NA, color = "grey50", linewidth = 0.3))
    }
    NULL
  }, error = function(e) NULL)
}


cartographier_generations <- function(pheno,
                                      date_carte = DATE_FIN,
                                      titre      = NULL,
                                      fichier    = NULL) {
  date_carte <- as.Date(date_carte)
  xlim <- c(-5.2, 9.7)
  ylim <- c(41.2, 51.2)
  
  rst_gen  <- tryCatch(get_generations_rst(pheno, as.character(date_carte)),
                       error = function(e) NULL)
  gen_vide <- raster_est_vide(rst_gen)
  if (!gen_vide) {
    # CORRECTION v3 [2] : on teste uniquement si les valeurs discretes sont presentes.
    # L'ancien test `|| is.numeric(vals)` etait toujours TRUE (barrks stocke les
    # generations en entiers numeriques) et bloquait l'affichage des generations.
    vals_uniques <- unique(na.omit(as.integer(values(rst_gen))))
    if (length(vals_uniques) == 0) gen_vide <- TRUE
  }
  
  if (gen_vide) {
    catg("  [{date_carte}] Pas encore de generation — affichage onset")
    dates_onset <- as.Date(time(pheno$onset))
    rst_onset   <- pheno$onset[[which.min(abs(dates_onset - date_carte))]]
    if (raster_est_vide(rst_onset)) {
      warning("Pas d'onset pour ", date_carte); return(invisible(NULL))
    }
    titre_auto <- titre %||% glue("Vol printanier P. chalcographus — {date_carte}")
    p <- ggplot() +
      geom_spatraster(data = project(rst_onset, CRS_WGS)) +
      fond_france() +
      scale_fill_gradient(low = "lightyellow", high = "darkorange3",
                          na.value = "grey90", name = "Essaimage\n(0=non, 1=oui)",
                          limits = c(0, 1)) +
      coord_sf(xlim = xlim, ylim = ylim) +
      labs(title    = titre_auto,
           subtitle = glue("Seuil vol : Tmax >= {SEUIL_VOL_TMAX} C | CHAPY | SAFRAN GeoSAS"),
           caption  = "Source : Meteo France — CNPF/CRPF AuRA", x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(plot.title   = element_text(face = "bold"),
            plot.caption = element_text(size = 8, color = "grey50"))
  } else {
    catg("  [{date_carte}] Generations presentes")
    titre_auto <- titre %||% glue("P. chalcographus — France — {date_carte}")
    p <- ggplot() +
      geom_spatraster(data = project(rst_gen, CRS_WGS)) +
      fond_france() +
      scale_fill_manual(values = barrks_colors(), labels = barrks_labels(),
                        na.value = "grey90", name = "Generation") +
      coord_sf(xlim = xlim, ylim = ylim) +
      labs(title    = titre_auto,
           subtitle = "CHAPY (Ogris et al. 2020) | SAFRAN 8 km — GeoSAS INRAE",
           caption  = "Source : Meteo France — CNPF/CRPF AuRA", x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(plot.title   = element_text(face = "bold"),
            plot.caption = element_text(size = 8, color = "grey50"))
  }
  
  if (!is.null(fichier)) {
    ggsave(fichier, plot = p, width = 10, height = 8, dpi = 150)
    catg("  Carte : {fichier}")
  }
  p
}


courbe_onset <- function(pheno, fichier = NULL) {
  df <- data.frame(date    = as.Date(time(pheno$onset)),
                   pct_vol = global(pheno$onset, fun = "mean", na.rm = TRUE)$mean * 100)
  df <- df[df$pct_vol > 0 & !is.nan(df$pct_vol), ]
  if (nrow(df) == 0) { cat("  Vol non encore declenche\n"); return(invisible(NULL)) }
  
  p <- ggplot(df, aes(date, pct_vol)) +
    geom_area(fill = "darkorange3", alpha = 0.15) +
    geom_line(color = "darkorange3", linewidth = 1.2) +
    geom_vline(xintercept = as.numeric(DATE_FIN), linetype = "dashed", color = "grey40") +
    annotate("text", x = DATE_FIN, y = max(df$pct_vol) * .85,
             label = "Fin donnees", hjust = -0.1, size = 3, color = "grey40") +
    scale_y_continuous(limits = c(0, 100), labels = percent_format(scale = 1)) +
    labs(title    = "Progression du vol printanier — P. chalcographus — France",
         subtitle = glue("% mailles ayant atteint Tmax >= {SEUIL_VOL_TMAX} C | {ANNEE}"),
         x = NULL, y = "% mailles avec vol printanier",
         caption  = "CHAPY | SAFRAN GeoSAS — Source : Meteo France — CNPF/CRPF AuRA") +
    theme_minimal(base_size = 12) +
    theme(plot.title   = element_text(face = "bold"),
          plot.caption = element_text(size = 8, color = "grey50"))
  
  if (!is.null(fichier)) ggsave(fichier, plot = p, width = 12, height = 5, dpi = 150)
  p
}


courbe_tmax_zone <- function(meteo, fichier = NULL) {
  df <- data.frame(
    date  = as.Date(time(meteo$tmax)),
    tmax  = global(meteo$tmax,  fun = "mean", na.rm = TRUE)$mean,
    tmean = global(meteo$tmean, fun = "mean", na.rm = TRUE)$mean
  )
  p <- ggplot(df, aes(date)) +
    geom_line(aes(y = tmean, color = "Tmean"), linewidth = 0.7, alpha = 0.8) +
    geom_line(aes(y = tmax,  color = "Tmax"),  linewidth = 1.0) +
    geom_hline(yintercept = SEUIL_VOL_TMAX, linetype = "dashed",
               color = "darkorange3", linewidth = 0.8) +
    annotate("text", x = min(df$date) + 5, y = SEUIL_VOL_TMAX + 0.8,
             label = glue("Seuil vol ({SEUIL_VOL_TMAX} C)"),
             color = "darkorange3", size = 3, hjust = 0) +
    geom_vline(xintercept = as.numeric(DATE_FIN), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Tmax" = "firebrick", "Tmean" = "steelblue"),
                       name = "Temperature") +
    labs(title    = "Temperatures moyennes nationales — France",
         subtitle = glue("Moyenne spatiale SAFRAN 8 km | {ANNEE}"),
         x = NULL, y = "Temperature (C)",
         caption  = "CHAPY | SAFRAN GeoSAS — Source : Meteo France — CNPF/CRPF AuRA") +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(face = "bold"),
          plot.caption    = element_text(size = 8, color = "grey50"),
          legend.position = "top")
  
  if (!is.null(fichier)) ggsave(fichier, plot = p, width = 12, height = 5, dpi = 150)
  p
}


courbe_photoperiode <- function(meteo, fichier = NULL) {
  df  <- data.frame(date      = as.Date(time(meteo$daylength)),
                    daylength = global(meteo$daylength, fun = "mean", na.rm = TRUE)$mean)
  d_d <- df$date[df$daylength < SEUIL_DIAPAUSE &
                   df$date > as.Date(glue("{ANNEE}-06-01"))][1]
  
  p <- ggplot(df, aes(date, daylength)) +
    geom_line(color = "royalblue", linewidth = 1.2) +
    geom_hline(yintercept = SEUIL_DIAPAUSE, linetype = "dashed",
               color = "purple", linewidth = 0.8) +
    annotate("text", x = min(df$date) + 5, y = SEUIL_DIAPAUSE + 0.3,
             label = glue("Seuil diapause ({SEUIL_DIAPAUSE} h)"),
             color = "purple", size = 3, hjust = 0) +
    {if (!is.na(d_d)) list(
      geom_vline(xintercept = as.numeric(d_d), linetype = "dotted",
                 color = "purple", linewidth = 0.8),
      annotate("text", x = d_d, y = 10, label = format(d_d, "%d %b"),
               hjust = -0.1, size = 3, color = "purple")
    )} +
    labs(title    = "Photopériode et seuil de diapause — France (moyenne nationale)",
         subtitle = glue("Diapause initiee si photopériode < {SEUIL_DIAPAUSE} h | {ANNEE}"),
         x = NULL, y = "Duree du jour (h)",
         caption  = "CHAPY (Ogris et al. 2020) — CNPF/CRPF AuRA") +
    theme_minimal(base_size = 12) +
    theme(plot.title   = element_text(face = "bold"),
          plot.caption = element_text(size = 8, color = "grey50"))
  
  if (!is.null(fichier)) ggsave(fichier, plot = p, width = 12, height = 5, dpi = 150)
  p
}


courbe_gen1 <- function(pheno, fichier = NULL) {
  df <- data.frame(
    date    = as.Date(time(pheno$development$gen_1)),
    dev_pct = global(pheno$development$gen_1, fun = "mean", na.rm = TRUE)$mean * 100
  )
  df <- df[!is.nan(df$dev_pct) & !is.na(df$dev_pct), ]
  if (nrow(df) == 0) { cat("  Gen 1 pas encore demarree\n"); return(invisible(NULL)) }
  
  p <- ggplot(df, aes(date, dev_pct)) +
    geom_area(fill = "darkorange3", alpha = 0.15) +
    geom_line(color = "darkorange3", linewidth = 1.2) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    scale_y_continuous(limits = c(0, 110)) +
    labs(title    = "Developpement generation 1 — P. chalcographus — France",
         subtitle = glue("% de completion moyen (France entiere) | {ANNEE}"),
         x = NULL, y = "Developpement gen 1 (%)",
         caption  = "CHAPY (Ogris et al. 2020) | SAFRAN — CNPF/CRPF AuRA") +
    theme_minimal(base_size = 12) +
    theme(plot.title   = element_text(face = "bold"),
          plot.caption = element_text(size = 8, color = "grey50"))
  
  if (!is.null(fichier)) ggsave(fichier, plot = p, width = 12, height = 5, dpi = 150)
  p
}


#' Diagramme de developpement pour une station ponctuelle (barrks natif)
diagramme_station <- function(pheno, lon, lat, nom = "Station", fichier = NULL) {
  catg(">> Diagramme {nom}...")
  pt  <- vect(data.frame(x = lon, y = lat), geom = c("x","y"), crs = CRS_WGS) |>
    project(CRS_L93)
  idx <- cellFromXY(pheno$onset[[1]], crds(pt))
  if (is.na(idx)) { warning("Station hors emprise."); return(invisible(NULL)) }
  p <- plot_development_diagram(pheno, stations_create(nom, idx),
                                xlim = range(as.Date(time(pheno$onset))))
  if (!is.null(fichier)) {
    png(fichier, width = 1600, height = 600, res = 150)
    print(p); dev.off()
    catg("  Diagramme : {fichier}")
  }
  invisible(p)
}


# =============================================================================
# 14. SYNTHESE CONSOLE
# =============================================================================

#' Rapport statistique national, structure comme celui de phenips_v3.
#' @param meteo optionnel : active le bloc temperature maximale
synthese_console <- function(pheno, meteo = NULL) {
  cat("\n", strrep("=", 68), "\n")
  cat("  SYNTHESE CHAPY — Pityogenes chalcographus — France metropolitaine\n")
  cat(strrep("=", 68), "\n")
  
  date_ref <- date_fin_chapy(pheno)
  
  # [3] v3 : get_onset_rst() plutot que pheno$onset[[nlyr]]
  onset_final <- tryCatch(
    get_onset_rst(pheno, as.character(date_ref)),
    error = function(e) pheno$onset[[nlyr(pheno$onset)]]
  )
  
  v        <- values(onset_final)
  n_total  <- ncell(onset_final)
  n_valide <- sum(!is.na(v))
  pct      <- function(n) if (n_valide > 0) round(100 * n / n_valide, 1) else NA_real_
  
  catg("  Date fin donnees : {format(date_ref, '%d/%m/%Y')}")
  catg("  Cellules         : {n_valide} valides / {n_total} total\n")
  
  # --- Vol printanier ---------------------------------------------------------
  # [v3.1] v > 0 fonctionne quel que soit le type de stockage barrks :
  # une valeur logique se coerce en 0/1, donc TRUE > 0 vaut TRUE.
  n_vol <- sum(v > 0, na.rm = TRUE)
  cat("  --- Vol printanier ---\n")
  catg("  Seuil            : Tmax >= {SEUIL_VOL_TMAX} C")
  catg("  Mailles avec vol : {n_vol} ({pct(n_vol)} %)")
  
  envol <- raster_premier_envol_doy(pheno)
  doy   <- na.omit(values(envol$rst))
  if (length(doy) > 0) {
    origine <- as.Date(glue("{ANNEE}-01-01"))
    q       <- quantile(doy, c(0, 0.5, 1), names = FALSE)
    d_txt   <- format(origine + q - 1L, "%d/%m/%Y")
    catg("  Envol le + precoce : {d_txt[1]} (DOY {round(q[1])})")
    catg("  Envol median       : {d_txt[2]} (DOY {round(q[2])})")
    catg("  Envol le + tardif  : {d_txt[3]} (DOY {round(q[3])})")
  }
  cat("\n")
  
  # --- Temperature maximale ---------------------------------------------------
  if (!is.null(meteo) && !is.null(meteo$tmax)) {
    t_vals <- values(meteo$tmax[[nlyr(meteo$tmax)]])
    n_seuil <- sum(t_vals >= SEUIL_VOL_TMAX, na.rm = TRUE)
    cat("  --- Temperature maximale (dernier jour) ---\n")
    catg("  Moyenne : {round(mean(t_vals, na.rm = TRUE), 1)} C")
    catg("  Min     : {round(min(t_vals,  na.rm = TRUE), 1)} C")
    catg("  Max     : {round(max(t_vals,  na.rm = TRUE), 1)} C")
    catg("  >= seuil vol ({SEUIL_VOL_TMAX} C) : {n_seuil} mailles ({pct(n_seuil)} %)\n")
  }
  
  # --- Developpement generation 1 ---------------------------------------------
  gen1 <- get_gen1_chapy(pheno)
  if (!is.null(gen1)) {
    dev_j <- values(gen1[[nlyr(gen1)]]) * 100
    n_c   <- sum(dev_j >= 100, na.rm = TRUE)
    cat("  --- Developpement generation 1 ---\n")
    catg("  Moyen          : {round(mean(dev_j, na.rm = TRUE), 1)} %")
    catg("  Max            : {round(max(dev_j,  na.rm = TRUE), 1)} %")
    catg("  Gen1 complete  : {n_c} mailles ({pct(n_c)} %)\n")
  } else {
    cat("  [INFO] Gen1 non disponible (normal avant mi-saison).\n\n")
  }
  
  # --- Generations completes --------------------------------------------------
  rst_gen <- get_gen_rst_chapy(pheno, date_ref)
  if (!rst_vide(rst_gen)) {
    # [E] normalisation des deux conventions barrks avant comptage
    gen_vals <- normaliser_generation(na.omit(values(rst_gen)))
    gen_vals <- gen_vals[!is.na(gen_vals)]
    cat("  --- Generations presentes ---\n")
    for (g in sort(unique(gen_vals))) {
      n_g <- sum(gen_vals == g)
      catg("  {libelle_generation(g)} : {n_g} mailles ({pct(n_g)} %)")
    }
    cat("\n")
  }
  
  cat(strrep("=", 68), "\n\n")
  invisible(TRUE)
}


# =============================================================================
# 15. PIPELINE PRINCIPAL
# =============================================================================

#' Pipeline complet France entiere — mise a jour mensuelle.
#' Etapes : cache SAFRAN -> assemblage -> daylength -> CHAPY -> sorties.
#'
#' @param mode         "max" (defaut) ou "mean" — scenario thermique CHAPY
#' @param dates_cartes dates des cartes mensuelles ; NULL = fin de chaque mois
#' @param force_dl     TRUE = purge le cache et retelecharge tout
#' @return liste(pheno, cartes)
pipeline_chapy_france <- function(mode         = "max",
                                  dates_cartes = NULL,
                                  force_dl     = FALSE) {
  cat("\n", strrep("=", 68), "\n")
  cat("  CHAPY v4 — Pityogenes chalcographus — FRANCE ENTIERE\n")
  cat("  Source : API GeoSAS INRAE / SAFRAN — Source : Meteo France\n")
  cat(strrep("=", 68), "\n\n")
  t0 <- proc.time()
  
  # --- 1. Telecharger / mettre a jour le cache par dalle ---
  if (force_dl) {
    fichiers_cache <- list.files(DIR_CACHE, pattern = "\\.csv$", full.names = TRUE)
    if (length(fichiers_cache) > 0) {
      file.remove(fichiers_cache)
      catg("  {length(fichiers_cache)} fichiers cache supprimes (force_dl=TRUE)")
    }
  }
  recuperer_safran()
  
  # --- 2. Lire et assembler ---
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  # [5] extraire les dates depuis l'attribut avant rm(dt) : evite le crash si
  # time() retourne NULL sur un raster sans axe temporel
  dates_serie <- attr(meteo, "dates")
  if (is.null(dates_serie)) dates_serie <- as.Date(time(meteo$tmean))
  rm(dt); gc()
  
  # --- 3. Daylength via create_daylength_rst() barrks ([1]) ---
  meteo$daylength <- calculer_daylength(RST_FRANCE, dates_serie)
  
  # --- 4. Modelisation CHAPY ---
  pheno <- modeliser_chapy(meteo, mode = mode)
  
  # --- 5. Sorties format PHENIPS : GeoTIFF + GeoJSON + cartes nationales ---
  exporter_rasters_chapy(pheno, meteo)
  produire_geojson_generations_chapy(pheno)
  produire_geojson_envol_chapy(pheno)
  generer_cartes_nationales_chapy(pheno, meteo)
  
  # --- 6. Cartes mensuelles de generations (specifique CHAPY) ---
  if (is.null(dates_cartes)) {
    mois_dispo   <- seq(DATE_DEB, DATE_FIN, by = "month")
    dates_cartes <- unique(c(
      as.Date(sapply(mois_dispo,
                     function(m) min(ceiling_date(m, "month") - 1, DATE_FIN))),
      DATE_FIN
    ))
  }
  catg("\n>> Cartes mensuelles de generations...")
  cartes <- lapply(as.Date(dates_cartes), function(d) {
    cartographier_generations(
      pheno, date_carte = d,
      fichier = file.path(DIR_OUT, glue("chapy_france_{format(d, '%Y%m%d')}.png"))
    )
  })
  
  # --- 7. Courbes temporelles nationales (specifique CHAPY) ---
  catg("\n>> Courbes temporelles...")
  courbe_onset(pheno,        fichier = file.path(DIR_OUT, "chapy_france_onset.png"))
  courbe_tmax_zone(meteo,    fichier = file.path(DIR_OUT, "chapy_france_tmax.png"))
  courbe_photoperiode(meteo, fichier = file.path(DIR_OUT, "chapy_france_photoperiode.png"))
  courbe_gen1(pheno,         fichier = file.path(DIR_OUT, "chapy_france_gen1.png"))
  
  # --- 8. Synthese console (meteo encore en memoire pour le bloc Tmax) ---
  synthese_console(pheno, meteo)
  rm(meteo); gc()
  
  duree         <- (proc.time() - t0)[["elapsed"]]
  prochaine_maj <- format(floor_date(DATE_FIN + 40, "month") + 5, "%d/%m/%Y")
  catg("  Duree totale : {round(duree / 60, 1)} min")
  catg("  Prochaine mise a jour recommandee : apres le {prochaine_maj}")
  catg("  Sorties : {DIR_OUT}/\n")
  
  invisible(list(pheno = pheno, cartes = cartes))
}


# =============================================================================
# 16. LANCEMENT
# =============================================================================

# --- Test connexion API sur un point unique ---
if (FALSE) {
  url_test <- paste0(
    "https://api.geosas.fr/edr/collections/safran-isba/position",
    "?parameter-name=T_Q",
    "&coords=POINT(700000%206440000)",
    "&crs=EPSG:2154",
    "&datetime=", ANNEE, "-01-01/", ANNEE, "-01-07",
    "&f=CSV"
  )
  resp <- request(url_test) |> req_perform()
  cat(resp_body_string(resp))
}

# --- Test telechargement d'une dalle ---
if (FALSE) {
  d_test <- DALLES[1, ]
  telecharger_dalle(d_test)
  fread(chemin_cache(d_test$id)) |> head()
}

# --- Test lecture + rasterisation ---
if (FALSE) {
  recuperer_safran()
  dt_test <- lire_cache_complet()
  print(dt_test[date == min(date)][1:5])
}

# --- Verification des sorties apres un run ---
# Liste les fichiers produits, tries par type, avec leur taille.
if (FALSE) {
  fichiers <- list.files(DIR_OUT, pattern = "\\.(tif|geojson|png)$", full.names = TRUE)
  info     <- data.frame(
    fichier = basename(fichiers),
    type    = toupper(tools::file_ext(fichiers)),
    ko      = round(file.size(fichiers) / 1024, 0),
    row.names = NULL
  )
  print(info[order(info$type, info$fichier), ])
}

# --- Regeneration des sorties sans relancer la modelisation ---
# Utile pour ajuster une palette ou un libelle : reutilise l'objet pheno.
if (FALSE) {
  exporter_rasters_chapy(resultats$pheno)          # meteo optionnel
  produire_geojson_generations_chapy(resultats$pheno)
  produire_geojson_envol_chapy(resultats$pheno)
  generer_cartes_nationales_chapy(resultats$pheno)
  synthese_console(resultats$pheno)
}

# --- Diagrammes de developpement par station ---
if (FALSE) {
  diagramme_station(resultats$pheno, 3.087, 45.777, "Clermont-Ferrand",
                    file.path(DIR_OUT, "diagramme_clermont.png"))
  diagramme_station(resultats$pheno, 2.347, 48.859, "Paris",
                    file.path(DIR_OUT, "diagramme_paris.png"))
  diagramme_station(resultats$pheno, 7.268, 47.740, "Strasbourg",
                    file.path(DIR_OUT, "diagramme_strasbourg.png"))
  diagramme_station(resultats$pheno, 5.921, 45.173, "Grenoble",
                    file.path(DIR_OUT, "diagramme_grenoble.png"))
}

# --- PIPELINE COMPLET FRANCE ENTIERE ---
# DESACTIVE : la version production passe par R/pipeline_chapy.R (charge le
# jeu SAFRAN partage, appelle modeliser_chapy() directement).
if (FALSE) {
  resultats <- pipeline_chapy_france(mode = "max")
}