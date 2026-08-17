# =============================================================================
# MONOCHAMUS galloprovincialis — France metropolitaine
# v1 — Modele DJ10 + bivoltinisme | Source : SAFRAN via API GéoSAS INRAE
# -----------------------------------------------------------------------------
# Contexte : crise nematode du pin (Bursaphelenchus xylophilus) en France 2024-
#
# Source meteo : API OGC EDR GéoSAS INRAE
#   URL        : https://api.geosas.fr/edr/collections/safran-isba/cube
#   Variables  : T_Q (Tmoy), TSUP_H_Q (Tmax), TINF_H_Q (Tmin)
#   Resolution : grille SAFRAN 8x8 km, France metropolitaine, 1958-2026
#   Avantage   : requete bbox en L93 directe, pas de point par point
#   Limite     : donnees disponibles jusqu'a J-5 environ (pas de temps reel)
#
# Strategie France entiere :
#   - Decoupage de la metropole en dalles de 200x200 km (L93)
#   - 1 requete API par dalle et par variable -> CSV
#   - Assemblage en raster national SAFRAN natif (8 km, L93)
#   - Cache local : ne re-telecharge que les jours manquants (relance hebdo)
#
# Modele pheno : DJ10 (Degres-Jours base 10C) depuis le 1er janvier
#   Seuils     : emergence=320, pic=550, fin=950 DJ10
#   Bivolt [A] : DJ10 annuel projete >= 1800
#   Bivolt [B] : +1000 DJ depuis pic gen1 avant le 1er octobre
#
# Relance      : hebdomadaire (cron ou manuel), idempotent (cache)
# Sorties      : GeoTIFF nationaux L93 + PNGs de synthese
#
# References   : Sousa et al. 2011, Naves et al. 2006, EPPO PM 9/5(1)
#                GéoSAS INRAE : https://geosas.fr/web/?page_id=6345
# Auteur       : CNPF / CRPF AuRA — avril 2026
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

# install.packages(c("terra","httr2","ggplot2","tidyterra",
#                    "lubridate","glue","sf","scales","data.table"))

library(terra)
library(httr2)
library(ggplot2)
library(tidyterra)
library(lubridate)
library(glue)
library(sf)
library(scales)
library(data.table)   # lecture CSV rapide


# =============================================================================
# 1. PARAMETRES GLOBAUX
# =============================================================================

# --- Periode ---
ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))

# SAFRAN : donnees disponibles jusqu'a ~J-5
# On interroge jusqu'a hier pour etre sur d'avoir des donnees completes
DATE_FIN_SAFRAN <- Sys.Date() - 5
if (DATE_FIN_SAFRAN < DATE_DEB) DATE_FIN_SAFRAN <- DATE_DEB

# --- API GéoSAS ---
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
GEOSAS_VARS <- "T_Q,TSUP_H_Q,TINF_H_Q"  # Tmoy, Tmax, Tmin

# --- Grille SAFRAN France metropolitaine (L93) ---
# Emprise approx. metropole en L93
FRANCE_XMIN <- 100000
FRANCE_XMAX <- 1200000
FRANCE_YMIN <- 6050000
FRANCE_YMAX <- 7150000
SAFRAN_RES  <- 8000        # 8 km

# Decoupage en dalles pour eviter les timeouts API
# 200 km = 25 tuiles SAFRAN par cote -> ~625 tuiles par dalle
DALLE_KM    <- 200
DALLE_M     <- DALLE_KM * 1000

CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

# --- Parametres biologiques ---
PARAM_MONO <- list(
  base_temp        = 10.0,
  seuil_vol_tmax   = 13.0,
  dj_emergence     = 320,
  dj_pic           = 550,
  dj_fin_vol       = 950,
  dj_risque_bivolt = 1800,
  dj_gen2_seuil    = 1000,
  date_limite_gen2 = "10-01"
)

# --- Repertoires ---
DIR_OUT   <- "data"
DIR_CACHE <- file.path(DIR_OUT, "cache_safran")
dir.create(DIR_OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

cat("=== Monochamus galloprovincialis — France metropolitaine (SAFRAN) ===\n")
cat(glue("  Periode     : {DATE_DEB} -> {DATE_FIN_SAFRAN}\n"))
cat(glue("  Source      : GéoSAS INRAE | safran-isba | {GEOSAS_VARS}\n"))
cat(glue("  Resolution  : {SAFRAN_RES/1000} km | Dalles : {DALLE_KM}x{DALLE_KM} km\n"))
cat(glue("  Cache       : {DIR_CACHE}/\n\n"))


# =============================================================================
# 2. GRILLE DE DALLES FRANCE
# =============================================================================

# =============================================================================
# 3. RASTER DE REFERENCE FRANCE (grille SAFRAN native)
# =============================================================================

# =============================================================================
# 4. TELECHARGEMENT SAFRAN PAR DALLE (avec cache)
# =============================================================================
#
# Format CSV retourne par GéoSAS :
#   colonnes : date, x, y, T_Q, TSUP_H_Q, TINF_H_Q
#   x, y     : coordonnees centroide tuile SAFRAN en L93
#   date     : YYYY-MM-DD
#
# Cache : un fichier CSV par dalle et par annee
#   -> si le fichier existe et couvre DATE_FIN_SAFRAN, on ne retelecharge pas
#   -> sinon on retelecharge depuis DATE_DEB (simpler que le diff partiel)
# =============================================================================

# =============================================================================
# 5. ASSEMBLAGE RASTER NATIONAL DEPUIS LES CSV
# =============================================================================
#
# Pour chaque jour, on lit tous les CSV de dalle, on extrait les valeurs
# du jour, on rasterise sur RST_FRANCE.
# Optimisation : on lit tous les CSV une seule fois en memoire (data.table),
# puis on pivote par date.
# =============================================================================

# =============================================================================
# 6. MODELE DJ10 + BIVOLTINISME
# (identique au script Massif Central mais adapte France entiere)
# =============================================================================

calculer_dj <- function(meteo, params = PARAM_MONO) {
  cat(">> Calcul DJ10 + module bivoltinisme...\n")
  
  tmean <- meteo$tmean
  tmax  <- meteo$tmax
  dates <- as.Date(time(tmean))
  n_j   <- nlyr(tmean)
  n_c   <- ncell(tmean)
  
  date_limite <- as.Date(glue("{year(dates[1])}-{params$date_limite_gen2}"))
  
  # Matrices [ncell x nlyr] via values()
  mat_tmean <- values(tmean)
  mat_tmax  <- values(tmax)
  cat(glue("  Matrices : {nrow(mat_tmean)} cellules x {ncol(mat_tmean)} jours\n"))
  
  # Masque mer/hors-France : cellules avec NA sur TOUTE la serie
  masque_na <- rowSums(!is.na(mat_tmean)) == 0
  cat(glue("  Cellules hors domaine (mer/NA) : {sum(masque_na)}/{n_c}\n"))
  
  # Cellules avec NAs partiels (donnees manquantes sur certains jours)
  n_partial <- sum(rowSums(is.na(mat_tmean)) > 0 & !masque_na)
  if (n_partial > 0)
    cat(glue("  Cellules avec NAs partiels : {n_partial} -> traites comme 0 DJ\n"))
  
  # DJ journalier : NA de temperature = 0 contribution (pas de chaleur accumulee)
  # IMPORTANT : on ne met a 0 que les NAs des cellules qui ont AU MOINS une valeur
  # Les cellules 100% NA (mer) restent NA et sont filtrees par masque_na
  dj_j <- pmax(mat_tmean - params$base_temp, 0, na.rm = FALSE)
  dj_j[is.na(dj_j) & !masque_na] <- 0  # NAs partiels -> 0 DJ (conservateur)
  
  # Cumul DJ10 par cellule [ncell x nlyr]
  # apply(, 1, cumsum) sur [ncell x nlyr] retourne [nlyr x ncell] -> t()
  dj_cum <- t(apply(dj_j, 1, cumsum))
  if (!is.matrix(dj_cum)) dj_cum <- matrix(dj_cum, nrow = n_c, ncol = n_j)
  # Remettre NA pour les cellules hors domaine
  dj_cum[masque_na, ] <- NA_real_
  
  # --- Diagnostic de coherence sur un point de reference ---
  # Bordeaux L93 ~ x=368000, y=6421000
  cell_ref <- cellFromXY(tmean[[1]], matrix(c(368000, 6421000), ncol = 2))
  if (!is.na(cell_ref) && !masque_na[cell_ref]) {
    tmoy_ref <- mat_tmean[cell_ref, ]
    dj_ref   <- dj_cum[cell_ref, ]
    cat(glue("  Verification Bordeaux (cellule {cell_ref}) :\n"))
    cat(glue("    Tmoy janvier (moy) : {round(mean(tmoy_ref[1:31],  na.rm=TRUE), 1)} C\n"))
    cat(glue("    Tmoy mars   (moy)  : {round(mean(tmoy_ref[60:90], na.rm=TRUE), 1)} C\n"))
    cat(glue("    DJ10 cumule J-fin   : {round(dj_ref[n_j], 0)} DJ10\n"))
  }
  
  vol_mat <- ifelse(mat_tmax >= params$seuil_vol_tmax, 1L, 0L)
  vol_mat[masque_na, ] <- NA_integer_
  
  # --- [A] Risque bivolt annuel ---
  cat("  [A] Risque bivolt annuel...\n")
  dj_fin_serie   <- dj_cum[, n_j]       # DJ cumule a la derniere date [ncell]
  date_fin_serie <- max(dates)
  jours_restants <- as.integer(as.Date(glue("{year(dates[1])}-12-31")) - date_fin_serie)
  
  if (jours_restants > 0 && jours_restants < 365) {
    # Rythme moyen des 30 derniers jours (dj_j est deja nettoye des NAs partiels)
    fenetre   <- min(30, n_j)
    rythme_dj <- rowMeans(dj_j[, (n_j - fenetre + 1):n_j, drop = FALSE], na.rm = TRUE)
    rythme_dj[masque_na] <- NA_real_
    mois_fin  <- month(date_fin_serie)
    coeff     <- ifelse(mois_fin >= 9, 0.4, ifelse(mois_fin >= 7, 0.7, 1.0))
    dj_projete <- dj_fin_serie + rythme_dj * jours_restants * coeff
    cat(glue("    Projection +{jours_restants}j (coeff={coeff}) vers 31-12\n"))
    cat(glue("    DJ10 projete moyen (hors mer) : ",
             "{round(mean(dj_projete[!masque_na], na.rm=TRUE), 0)} DJ10\n"))
    cat(glue("    DJ10 projete max   (hors mer) : ",
             "{round(max(dj_projete[!masque_na],  na.rm=TRUE), 0)} DJ10\n"))
  } else {
    dj_projete <- dj_fin_serie
  }
  
  risque_bivolt_vec <- as.integer(dj_projete >= params$dj_risque_bivolt)
  risque_bivolt_vec[masque_na] <- NA_integer_
  n_risque <- sum(risque_bivolt_vec, na.rm = TRUE)
  cat(glue("    {n_risque} cellules a risque bivolt [A]\n"))
  
  # --- [B] Module bivolt complet ---
  cat("  [B] Module bivolt gen2...\n")
  idx_pic_gen1 <- apply(dj_cum, 1, function(row) {
    idx <- which(row >= params$dj_pic)
    if (length(idx) == 0) return(NA_integer_)
    idx[1]
  })
  
  # Date butoir : verifier si elle est dans la serie disponible
  # Si non (donnees s'arretent avant le 1er octobre), idx_limite = NA
  # -> les cellules seront "en cours" et non "avortees" prematurement
  dates_apres_butoir <- which(dates >= date_limite)
  if (length(dates_apres_butoir) > 0) {
    idx_limite <- dates_apres_butoir[1]
    cat(glue("    Date butoir ({date_limite}) presente dans la serie (jour {idx_limite})\n"))
  } else {
    idx_limite <- NA_integer_
    cat(glue("    Date butoir ({date_limite}) HORS serie (donnees jusqu'au {max(dates)})\n"))
    cat(glue("    -> Les cellules seront classees 'en cours' en fin de serie\n"))
  }
  
  gen2_etat_mat <- matrix(0L, nrow = n_c, ncol = n_j)
  
  cells_actives <- which(!is.na(idx_pic_gen1) & !masque_na)
  cat(glue("    {length(cells_actives)} cellules avec gen1 atteinte\n"))
  
  for (c in cells_actives) {
    ip      <- idx_pic_gen1[c]
    dj_gen2 <- 0.0
    
    for (j in seq(ip, n_j)) {
      dj_gen2 <- dj_gen2 + dj_j[c, j]
      if (is.na(dj_gen2)) break
      
      if (dj_gen2 >= params$dj_gen2_seuil) {
        # Seuil atteint avant le butoir -> gen2 viable
        gen2_etat_mat[c, j:n_j] <- 2L
        break
      } else if (!is.na(idx_limite) && j >= idx_limite) {
        # Butoir depasse sans atteindre le seuil -> gen2 avortee
        gen2_etat_mat[c, j:n_j] <- 3L
        break
      } else {
        # Accumulation en cours (butoir pas encore atteint ou hors serie)
        gen2_etat_mat[c, j] <- 1L
      }
    }
  }
  
  gen2_final_etat              <- gen2_etat_mat[, n_j]
  gen2_final_etat[masque_na]   <- NA_integer_
  gen2_viable_vec              <- as.integer(gen2_final_etat == 2L)
  
  cat(glue("    Gen2 viable : {sum(gen2_viable_vec, na.rm=TRUE)} cellules\n"))
  cat(glue("    Gen2 en cours : {sum(gen2_final_etat==1L, na.rm=TRUE)} cellules\n"))
  cat(glue("    Gen2 avortee  : {sum(gen2_final_etat==3L, na.rm=TRUE)} cellules\n"))
  
  # Etat gen1
  etat_mat <- matrix(0L, nrow = n_c, ncol = n_j)
  etat_mat[dj_cum >= params$dj_emergence & dj_cum < params$dj_pic]    <- 1L
  etat_mat[dj_cum >= params$dj_pic       & dj_cum < params$dj_fin_vol] <- 2L
  etat_mat[dj_cum >= params$dj_fin_vol]                                 <- 3L
  etat_mat[masque_na, ] <- NA_integer_
  
  cat(glue("  DJ10 moyen zone (hors mer) J-fin : ",
           "{round(mean(dj_cum[!masque_na, n_j], na.rm=TRUE), 0)} DJ10\n\n"))
  
  # Reconstruction SpatRasters
  to_rst <- function(mat, ref = tmean) {
    couches <- lapply(seq_len(n_j), function(i) {
      r <- ref[[1]]; values(r) <- mat[, i]; names(r) <- as.character(dates[i]); r
    })
    rst <- rast(couches); time(rst) <- dates; rst
  }
  to_rst1 <- function(vec, ref = tmean[[1]]) {
    r <- ref; values(r) <- vec; r
  }
  
  list(
    dj_cum        = to_rst(dj_cum),
    vol_actif     = to_rst(vol_mat),
    etat          = to_rst(etat_mat),
    dates         = dates,
    risque_bivolt = to_rst1(risque_bivolt_vec),
    dj_projete    = to_rst1(dj_projete),
    gen2_etat     = to_rst(gen2_etat_mat),
    gen2_viable   = to_rst1(gen2_viable_vec),
    gen2_final    = to_rst1(gen2_final_etat)
  )
}


# =============================================================================
# 7. EXPORT RASTERS
# =============================================================================

exporter_rasters <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Export GeoTIFF (L93, LZW)...\n")
  date_j  <- format(max(pheno$dates), "%Y%m%d")
  
  ecrire <- function(rst, nom) {
    ch <- file.path(dir_out, nom)
    writeRaster(rst, ch, overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    cat(glue("   {nom}\n"))
  }
  
  cat("  [1] DJ10 accumule :\n")
  # Serie complete (multi-bandes) — peut etre volumineuse (~600x750x122 float32)
  rst_serie <- pheno$dj_cum
  names(rst_serie) <- as.character(pheno$dates)
  ecrire(rst_serie, glue("mono_DJ10_serie_{ANNEE}_L93.tif"))
  
  # Couche J final
  rst_j <- pheno$dj_cum[[ nlyr(pheno$dj_cum) ]]
  names(rst_j) <- glue("DJ10_{date_j}")
  ecrire(rst_j, glue("mono_DJ10_jourJ_{date_j}_L93.tif"))
  
  cat("  [2] Bivolt [A] seuil annuel :\n")
  rst_proj <- pheno$dj_projete
  names(rst_proj) <- glue("DJ10proj_{ANNEE}")
  ecrire(rst_proj, glue("mono_bivoltA_DJ10projete_{ANNEE}_L93.tif"))
  
  rst_ra <- classify(pheno$risque_bivolt, cbind(NA, NA))
  names(rst_ra) <- glue("RisqueA_{ANNEE}")
  ecrire(rst_ra, glue("mono_bivoltA_risque01_{ANNEE}_L93.tif"))
  
  cat("  [3] Bivolt [B] gen2 :\n")
  rst_ge <- pheno$gen2_final
  names(rst_ge) <- glue("EtatGen2_{date_j}")
  ecrire(rst_ge, glue("mono_bivoltB_etatGen2_{date_j}_L93.tif"))
  
  rst_gv <- pheno$gen2_viable
  names(rst_gv) <- glue("Gen2viable_{date_j}")
  ecrire(rst_gv, glue("mono_bivoltB_gen2viable01_{date_j}_L93.tif"))
  
  cat(glue("\n  6 rasters exportes dans {dir_out}/\n"))
  cat("  Legende etat gen2 : NA=mer | 0=inactif | 1=en cours | 2=viable | 3=avortee\n\n")
}


# =============================================================================
# 8. EXPORT GEOJSON
# =============================================================================
#
# Deux GeoJSON en WGS84 (EPSG:4326), polygones dissous par valeur :
#
# [1] mono_stade_YYYYMMDD_WGS84.geojson
#     Champs : stade (int 0-3), stade_label (str), dj_seuil_min, dj_seuil_max
#
# [2] mono_bivoltB_gen2_YYYYMMDD_WGS84.geojson
#     Champs : statut_gen2 (int 0-3), statut_label (str)
#
# Note : as.polygons() avec dissolve=TRUE peut prendre 2-5 min sur la France
#        entiere. Les polygones sont simplifies (snap) pour alleger le fichier.
# =============================================================================

exporter_geojson <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Export GeoJSON (WGS84)...\n")
  date_j  <- format(max(pheno$dates), "%Y%m%d")
  idx_j   <- nlyr(pheno$dj_cum)
  
  # Labels pour les attributs
  stade_labels  <- c("0" = "Inactif",
                     "1" = "Pre-emergence",
                     "2" = "Vol actif",
                     "3" = "Fin de vol")
  stade_dj_min  <- c("0" = 0,   "1" = 320, "2" = 550, "3" = 950)
  stade_dj_max  <- c("0" = 320, "1" = 550, "2" = 950, "3" = 9999)
  # Intensite croissante (activite) : gris -> jaune -> orange -> rouge
  stade_couleurs <- c("0" = "#E8E8E8", "1" = "#FDE725FF",
                      "2" = "#F39C12", "3" = "#C0392B")

  gen2_labels   <- c("0" = "Inactif (gen1 non atteinte)",
                     "1" = "Gen2 en cours",
                     "2" = "Gen2 viable",
                     "3" = "Gen2 avortee")
  # Statuts non ordonnes (une gen2 "avortee" n'est pas "moins grave" qu'une
  # gen2 "en cours") : palette qualitative plutot que sequentielle.
  #   gris = inactif | jaune = en cours | rouge = viable (risque confirme)
  #   gris-bleu = avortee (fin de cycle sans descendance)
  gen2_couleurs <- c("0" = "#E8E8E8", "1" = "#FDE725FF",
                     "2" = "#C0392B",  "3" = "#7F8C8D")
  
  # Helper : raster categoriel -> polygones dissous -> GeoJSON WGS84
  raster_vers_geojson <- function(rst_categoriel, nom_champ, labels,
                                  couleurs = NULL,
                                  champs_extra = NULL, raster_extra = NULL,
                                  fichier = NULL) {
    
    cat(glue("  Conversion polygones {fichier}...\n"))
    
    # Masquer les NA (mer) explicitement avec une valeur sentinelle hors plage
    rst_clean <- rst_categoriel
    rst_clean[is.na(rst_clean)] <- -9999L
    
    # Polygonisation avec dissolution par valeur
    # snap = TRUE : agrège les cellules adjacentes de même valeur
    poly <- as.polygons(rst_clean, dissolve = TRUE, na.rm = FALSE)
    
    # Filtrer la valeur sentinelle (mer)
    poly <- poly[values(poly)[, 1] != -9999L, ]
    
    # Renommer le champ valeur
    names(poly)[1] <- nom_champ
    
    # Ajouter le label textuel
    valeurs <- as.character(as.integer(values(poly)[, 1]))
    poly$label <- unname(labels[valeurs])
    
    # Ajouter la couleur (calculee ici, pas cote app : coherence garantie
    # entre QGIS, l'app Shiny, ou tout autre client qui lit ce GeoJSON)
    if (!is.null(couleurs)) {
      poly$couleur <- unname(couleurs[valeurs])
    }
    
    # Champs extra (DJ seuils pour le stade)
    if (!is.null(champs_extra)) {
      for (nm in names(champs_extra)) {
        poly[[nm]] <- unname(champs_extra[[nm]][valeurs])
      }
    }
    
    # Bivoltisme [A] : rasters mono-couche independants de la classification
    # du stade (risque annuel 0/1 + DJ10 projete). Extraits par centroide de
    # chaque polygone dissous, en L93 (avant reprojection WGS84).
    if (!is.null(raster_extra)) {
      cent_l93 <- centroids(poly)
      for (nm in names(raster_extra)) {
        poly[[nm]] <- as.numeric(extract(raster_extra[[nm]], cent_l93)[, 2])
      }
    }
    
    # Ajouter la date
    poly$date_mise_a_jour <- as.character(max(pheno$dates))
    
    # Reprojection WGS84
    poly_wgs <- project(poly, CRS_WGS)
    
    # Export GeoJSON
    chemin <- file.path(dir_out, fichier)
    writeVector(poly_wgs, chemin,
                filetype = "GeoJSON", overwrite = TRUE)
    
    # Taille du fichier
    taille_ko <- round(file.size(chemin) / 1024)
    cat(glue("   {fichier} ({taille_ko} Ko, {nrow(poly_wgs)} polygones)\n"))
    invisible(poly_wgs)
  }
  
  # --- [1] Stade phenologique + bivoltisme [A] (risque annuel + DJ10 projete) ---
  cat("  [1] Stade phenologique gen1 + bivoltisme [A]...\n")
  rst_stade <- pheno$etat[[idx_j]]
  raster_vers_geojson(
    rst_categoriel = rst_stade,
    nom_champ      = "stade",
    labels         = stade_labels,
    couleurs       = stade_couleurs,
    champs_extra   = list(dj_min = stade_dj_min, dj_max = stade_dj_max),
    raster_extra   = list(
      risque_bivoltA = pheno$risque_bivolt,   # 0/1
      dj10_projete   = pheno$dj_projete       # DJ10 projetes au 31/12
    ),
    fichier        = glue("mono_stade_{date_j}_WGS84.geojson")
  )
  
  # --- [2] Statut gen2 bivolt B ---
  cat("  [2] Statut generation 2 [B]...\n")
  rst_gen2 <- pheno$gen2_final
  raster_vers_geojson(
    rst_categoriel = rst_gen2,
    nom_champ      = "statut_gen2",
    labels         = gen2_labels,
    couleurs       = gen2_couleurs,
    fichier        = glue("mono_bivoltB_gen2_{date_j}_WGS84.geojson")
  )
  
  cat(glue("\n  2 GeoJSON exportes dans {dir_out}/\n"))
  cat("  CRS : WGS84 (EPSG:4326) — ouvrable directement dans QGIS, Leaflet, MapLibre\n\n")
}


# =============================================================================
# 8. VISUALISATIONS NATIONALES (PNG)
# =============================================================================

STADE_LABELS <- c(
  "0" = "Inactif (< 320 DJ10)",
  "1" = "Pre-emergence (320-550)",
  "2" = "Vol actif (550-950)",
  "3" = "Fin de vol (> 950)"
)
STADE_COULEURS <- c(
  "0" = "grey92", "1" = "#FDE725FF",
  "2" = "#35B779FF", "3" = "#31688EFF"
)
GEN2_LABELS <- c(
  "0" = "Inactif", "1" = "Gen2 en cours",
  "2" = "Gen2 viable", "3" = "Gen2 avortee"
)
GEN2_COULEURS <- c(
  "0" = "grey92", "1" = "#FFD580",
  "2" = "#C0392B", "3" = "#85929E"
)

carte_france <- function(rst_l93, titre, legende, palette, breaks = NULL,
                         fichier = NULL, type = c("gradient","categoriel"),
                         date_label = NULL) {
  type <- match.arg(type)
  if (is.null(date_label)) date_label <- format(Sys.Date(), "%d/%m/%Y")
  rst_wgs <- project(rst_l93, CRS_WGS)
  
  p <- ggplot() +
    geom_spatraster(data = rst_wgs) +
    {
      if (type == "gradient") {
        scale_fill_gradientn(
          colours  = palette,
          values   = if (!is.null(breaks)) scales::rescale(breaks) else NULL,
          na.value = "lightblue",
          name     = legende,
          guide    = guide_colorbar(barwidth = 14, barheight = 0.8,
                                    title.position = "top")
        )
      } else {
        scale_fill_manual(
          values   = palette,
          na.value = "lightblue",
          name     = legende,
          drop     = FALSE
        )
      }
    } +
    labs(
      title    = titre,
      subtitle = glue("SAFRAN INRAE | {date_label} | {ANNEE}"),
      caption  = "Modele DJ10 | GéoSAS INRAE — CNPF/CRPF AuRA",
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
    cat(glue("  Carte : {fichier}\n"))
  }
  invisible(p)
}

generer_cartes <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Generation des cartes nationales...\n")
  date_j     <- format(max(pheno$dates), "%Y%m%d")
  date_label <- format(max(pheno$dates), "%d/%m/%Y")
  idx_j      <- nlyr(pheno$dj_cum)
  
  # 1. Cumul DJ10
  carte_france(
    pheno$dj_cum[[idx_j]],
    titre       = "Cumul DJ10 — Monochamus galloprovincialis — France",
    legende     = "DJ10 cumules",
    palette     = c("grey95","#FFF3B0","#FDE725FF","#35B779FF","#31688EFF"),
    breaks      = c(0, 200, 320, 550, 950),
    type        = "gradient",
    date_label  = date_label,
    fichier     = file.path(dir_out, glue("mono_carte_DJ10_{date_j}.png"))
  )
  
  # 2. Stade phenologique gen1
  rst_etat <- pheno$etat[[idx_j]]
  values(rst_etat) <- factor(as.integer(values(rst_etat)),
                             levels = 0:3, labels = STADE_LABELS)
  carte_france(
    rst_etat,
    titre       = "Stade phenologique — Monochamus galloprovincialis — France",
    legende     = "Stade gen1",
    palette     = setNames(STADE_COULEURS, STADE_LABELS),
    type        = "categoriel",
    date_label  = date_label,
    fichier     = file.path(dir_out, glue("mono_carte_stade_{date_j}.png"))
  )
  
  # 3. Risque bivolt [A]
  carte_france(
    pheno$dj_projete,
    titre       = "Risque bivoltinisme [A] — DJ10 projete au 31-dec — France",
    legende     = glue("DJ10 annuel projete\n(seuil risque = {PARAM_MONO$dj_risque_bivolt})"),
    palette     = c("grey95","#F9E79F","#F39C12","#C0392B","#7B241C"),
    breaks      = c(0, 900, 1400, 1800, 2200),
    type        = "gradient",
    date_label  = date_label,
    fichier     = file.path(dir_out, glue("mono_carte_bivoltA_{ANNEE}.png"))
  )
  
  # 4. Statut gen2 [B]
  rst_gen2 <- pheno$gen2_final
  values(rst_gen2) <- factor(as.integer(values(rst_gen2)),
                             levels = 0:3, labels = GEN2_LABELS)
  carte_france(
    rst_gen2,
    titre       = "Bivoltinisme [B] — Statut generation 2 — France",
    legende     = "Statut gen2",
    palette     = setNames(GEN2_COULEURS, GEN2_LABELS),
    type        = "categoriel",
    date_label  = date_label,
    fichier     = file.path(dir_out, glue("mono_carte_bivoltB_{date_j}.png"))
  )
  
  cat("  4 cartes nationales generees.\n\n")
}


# =============================================================================
# 9. SYNTHESE CONSOLE
# =============================================================================

synthese_console <- function(pheno) {
  cat("\n", strrep("=", 60), "\n")
  cat("  SYNTHESE — France metropolitaine\n")
  cat(strrep("=", 60), "\n")
  
  n_total <- ncell(pheno$risque_bivolt)
  n_valide <- sum(!is.na(values(pheno$risque_bivolt)))
  idx_j    <- nlyr(pheno$dj_cum)
  
  dj_mean <- mean(values(pheno$dj_cum[[idx_j]]), na.rm = TRUE)
  dj_min  <- min(values(pheno$dj_cum[[idx_j]]),  na.rm = TRUE)
  dj_max  <- max(values(pheno$dj_cum[[idx_j]]),  na.rm = TRUE)
  
  etat_j  <- values(pheno$etat[[idx_j]])
  n_inact <- sum(etat_j == 0L, na.rm = TRUE)
  n_pre   <- sum(etat_j == 1L, na.rm = TRUE)
  n_vol   <- sum(etat_j == 2L, na.rm = TRUE)
  n_fin   <- sum(etat_j == 3L, na.rm = TRUE)
  
  n_risqueA <- sum(values(pheno$risque_bivolt) == 1L, na.rm = TRUE)
  gen2_f    <- values(pheno$gen2_final)
  n_viable  <- sum(gen2_f == 2L, na.rm = TRUE)
  n_avorte  <- sum(gen2_f == 3L, na.rm = TRUE)
  
  pct <- function(n) round(100 * n / n_valide, 1)
  
  cat(glue("  Date     : {format(max(pheno$dates), '%d/%m/%Y')}\n"))
  cat(glue("  Cellules : {n_valide} valides / {n_total} total\n\n"))
  
  cat("  --- DJ10 cumule ---\n")
  cat(glue("  Moyen : {round(dj_mean,0)} | Min : {round(dj_min,0)} | Max : {round(dj_max,0)}\n\n"))
  
  cat("  --- Stades phenologiques gen1 ---\n")
  cat(glue("  Inactif       : {n_inact} cel. ({pct(n_inact)} %)\n"))
  cat(glue("  Pre-emergence : {n_pre}  cel. ({pct(n_pre)} %)\n"))
  cat(glue("  Vol actif     : {n_vol}  cel. ({pct(n_vol)} %)\n"))
  cat(glue("  Fin de vol    : {n_fin}  cel. ({pct(n_fin)} %)\n\n"))
  
  cat("  --- Bivoltinisme ---\n")
  cat(glue("  [A] Risque annuel (>={PARAM_MONO$dj_risque_bivolt} DJ) : ",
           "{n_risqueA} cel. ({pct(n_risqueA)} %)\n"))
  cat(glue("  [B] Gen2 viable  : {n_viable} cel. ({pct(n_viable)} %)\n"))
  cat(glue("  [B] Gen2 avortee : {n_avorte} cel. ({pct(n_avorte)} %)\n"))
  
  if (pct(n_viable) > 5)
    cat("\n  ** ALERTE : > 5% du territoire avec gen2 viable.\n")
  
  cat(strrep("=", 60), "\n\n")
}


# =============================================================================
# 10. PIPELINE PRINCIPAL
# =============================================================================

pipeline_france <- function() {
  cat("\n", strrep("=", 60), "\n")
  cat("  MONOCHAMUS — France metropolitaine — SAFRAN\n")
  cat(strrep("=", 60), "\n\n")
  t0 <- proc.time()
  
  # 1. Telecharger / mettre a jour le cache
  recuperer_safran()
  
  # 2. Lire et assembler
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  rm(dt); gc()  # liberer la RAM
  
  # 3. Modele DJ10 + bivolt
  pheno <- calculer_dj(meteo)
  rm(meteo); gc()
  
  # 4. Export rasters
  exporter_rasters(pheno)
  
  # 5. Export GeoJSON
  exporter_geojson(pheno)
  
  # 5. Cartes
  generer_cartes(pheno)
  
  # 6. Synthese console
  synthese_console(pheno)
  
  duree <- (proc.time() - t0)[["elapsed"]]
  cat(glue("  Duree totale : {round(duree/60, 1)} min\n"))
  cat(glue("  Sorties      : {DIR_OUT}/\n\n"))
  
  invisible(pheno)
}


# =============================================================================
# 11. LANCEMENT
# =============================================================================

# Test connexion API sur un point unique
if (FALSE) {
  url_test <- paste0(
    "https://api.geosas.fr/edr/collections/safran-isba/position",
    "?parameter-name=T_Q",
    "&coords=POINT(700000%206440000)",
    "&crs=EPSG:2154",
    "&datetime=2026-01-01/2026-01-10",
    "&f=CSV"
  )
  resp <- request(url_test) |> req_perform()
  cat(resp_body_string(resp))
}

# Lancement complet DESACTIVE : la version production passe par
# R/pipeline_monochamus.R (charge le jeu SAFRAN partage, appelle calculer_dj()
# directement). pipeline_france() ci-dessus reste utilisable en local/debug
# UNIQUEMENT si tu resources aussi R/safran_common.R avant (sinon
# recuperer_safran()/lire_cache_complet() sont introuvables ici).
# pheno <- pipeline_france()