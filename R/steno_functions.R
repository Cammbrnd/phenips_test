# =============================================================================
# STENO-DJ — Phenologie du stenographe (Ips sexdentatus) — France metropolitaine
# v1 — SAFRAN via API GeoSAS INRAE | Architecture dalles + cache CSV
# -----------------------------------------------------------------------------
# POURQUOI CE SCRIPT N'UTILISE PAS barrks
#   barrks n'implemente que des modeles calibres pour Ips typographus
#   (phenips, phenips-clim, rity, chapy, lange, bso). Il n'existe aucun modele
#   publie "cle en main" pour Ips sexdentatus. On implemente donc un modele
#   degres-jours explicite, calibre sur les parametres thermiques mesures en
#   laboratoire sur pin maritime par Pineau et al. (2017).
#
# PARAMETRES BIOLOGIQUES (sources)
#   Seuil inferieur de developpement  : 10.9 C   (Pineau et al. 2017)
#   Seuil superieur de developpement  : 36.0 C   (Pineau et al. 2017)
#   Optimum thermique                 : 29.0 C   (Pineau et al. 2017)
#   Constante thermique 1 generation  : 517 DJ   (Pineau et al. 2017)
#   Declenchement du 1er essaimage    : Tair ~20 C plusieurs jours consecutifs,
#                                       sans pluie ni gelee nocturne
#                                       (Ephytia/DSF, fiche Stenographe)
#   Nb de generations attendu         : 2/an en zone temperate, 3 en annee
#                                       chaude ou en zone meridionale
#                                       (Ephytia/DSF ; Pineau et al. 2017 pour
#                                        le Sud-Ouest)
#
#   /!\ A la difference de PHENIPS (typographe), il n'existe PAS de modele
#       valide de terrain pour le stenographe. Ce script produit un INDICATEUR
#       d'aide a la decision (dates de vol probables, nb de generations
#       potentielles), pas une prediction validee. Toute utilisation en appui
#       technique doit s'accompagner d'un piegeage de controle.
#
# -----------------------------------------------------------------------------
# Source meteo : API OGC EDR GeoSAS INRAE  (identique au script PHENIPS)
#   URL        : https://api.geosas.fr/edr/collections/safran-isba/cube
#   Variables  : T_Q (Tmoy), TINF_H_Q (Tmin), TSUP_H_Q (Tmax), SSI_Q (rayont.)
#   Resolution : grille SAFRAN 8x8 km | mise a jour mensuelle, ~J-5
#
#   >> Le cache CSV est PARTAGE avec le script PHENIPS (meme dossier par defaut,
#      meme requete, memes variables) : aucun retelechargement si tu as deja
#      lance phenips_france_safran_v3.R sur la meme annee.
#
# Auteur      : CNPF / CRPF AuRA — aout 2026
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

# install.packages(c("terra", "httr2", "ggplot2", "tidyterra",
#                    "lubridate", "glue", "sf", "scales", "data.table"))

library(terra)
library(httr2)
library(ggplot2)
library(tidyterra)
library(lubridate)
library(glue)
library(sf)
library(scales)
library(data.table)


# =============================================================================
# 1. PARAMETRES GLOBAUX
# =============================================================================

# --- Periode ---
ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))

DATE_FIN_SAFRAN <- Sys.Date() - 5
if (DATE_FIN_SAFRAN < DATE_DEB) DATE_FIN_SAFRAN <- DATE_DEB

# --- API GeoSAS ---
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
GEOSAS_VARS <- "T_Q,TINF_H_Q,TSUP_H_Q,SSI_Q"

# --- Grille SAFRAN France metropolitaine (L93) ---
FRANCE_XMIN <- 100000
FRANCE_XMAX <- 1200000
FRANCE_YMIN <- 6050000
FRANCE_YMAX <- 7150000
SAFRAN_RES  <- 8000

DALLE_KM <- 200
DALLE_M  <- DALLE_KM * 1000

CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

# -----------------------------------------------------------------------------
# RAYONNEMENT : detection automatique de l'unite  [CORRECTION vs PHENIPS v3]
# -----------------------------------------------------------------------------
# Doc SAFRAN : SSI_Q = "Rayonnement visible (cumul quotidien)" en J/cm2.
# Ordres de grandeur ATTENDUS en J/cm2 pour un cumul journalier :
#     hiver ~ 200-400 J/cm2   |   ete ~ 2000-2600 J/cm2
# soit en Wh/m2 (unite attendue par l'equation de temperature d'ecorce) :
#     hiver ~ 550-1100        |   ete ~ 5500-7200
#
# /!\ Le script PHENIPS v3 supposait une mediane brute de ~15-25 "J/cm2" et un
#     garde-fou SSI_CHECK_MAX = 50. Si l'API renvoie effectivement ~15-25, alors
#     la valeur est en MJ/m2 (25 MJ/m2 = 2500 J/cm2) et le facteur 2.7778 est
#     100x trop petit : l'effet du rayonnement sur la temperature sous ecorce
#     devient ~0.06 C au lieu de ~5 C, donc PHENIPS tourne en pratique en mode
#     "shaded". >> A verifier sur ton cache, voir diagnostiquer_rayonnement().
#
# On evite le probleme : l'unite est deduite de la magnitude observee.
detecter_facteur_rad <- function(mediane_brute) {
  # Retourne le facteur de conversion -> Wh/m2 et l'unite supposee
  if (is.na(mediane_brute) || mediane_brute <= 0) {
    return(list(facteur = 10000 / 3600, unite = "J/cm2 (defaut, non verifiable)"))
  }
  if (mediane_brute < 60) {
    # 0-60 : compatible MJ/m2 (mediane annuelle France ~ 11 MJ/m2)
    list(facteur = 1e6 / 3600,        unite = "MJ/m2")
  } else if (mediane_brute < 600) {
    # 60-600 : compatible W/m2 moyen journalier (mediane ~ 130 W/m2)
    list(facteur = 24,                unite = "W/m2 (moyenne journaliere)")
  } else {
    # > 600 : compatible J/cm2 (mediane annuelle France ~ 1100 J/cm2)
    list(facteur = 10000 / 3600,      unite = "J/cm2")
  }
}

# --- Repertoires ---
DIR_OUT   <- "data"
# Cache PARTAGE avec les 3 autres modeles (voir safran_common.R)
DIR_CACHE <- file.path(DIR_OUT, "cache_safran")
dir.create(DIR_OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 1bis. PARAMETRES BIOLOGIQUES DU STENOGRAPHE
# =============================================================================
# Tout est regroupe ici : c'est le seul bloc a toucher pour recalibrer.
# Les valeurs marquees [PUBLIE] sont issues de la litterature ;
# les valeurs marquees [EXPERT] sont des hypotheses de travail a calibrer
# sur piegeage local (pheromone : ipsdienol / cis-verbenol).

STENO <- list(
  # --- Developpement (Pineau et al. 2017, AFE, pin maritime) ---
  T_MIN_DEV    = 10.9,   # [PUBLIE] seuil inferieur de developpement (C)
  T_MAX_DEV    = 36.0,   # [PUBLIE] seuil superieur (cutoff horizontal) (C)
  T_OPT        = 29.0,   # [PUBLIE] optimum (informatif, non utilise en DJ lineaire)
  DJ_GEN       = 517,    # [PUBLIE] DJ pour un developpement complet oeuf -> adulte
  
  # --- Maturation / alimentation complementaire avant nouvel essaimage ---
  # Aucune valeur publiee pour I. sexdentatus. Le cycle effectif retenu est donc
  # DJ_GEN + DJ_MATURATION : les 517 DJ couvrent oeuf -> adulte, la maturation
  # couvre l'alimentation complementaire, la dispersion, le forage et
  # l'accouplement avant la ponte suivante. A 0, le modele enchaine les
  # generations instantanement et surestime le voltinisme.
  DJ_MATURATION = 100,   # [EXPERT] cycle effectif = 617 DJ
  
  # --- Essaimage : somme thermique prealable (structure PHENIPS) ---
  # Le declenchement combine DEUX conditions, comme dans Baier et al. (2007) :
  #   (a) une somme thermique sous ecorce accumulee depuis DOY_DEB_CUMUL
  #       doit atteindre DJ_VOL  -> donne le gradient spatial (latitude,
  #       altitude, continentalite)
  #   (b) la journee elle-meme doit etre volable : Tair max >= SEUIL_VOL,
  #       sans gelee -> donne la date exacte de declenchement
  # PHENIPS utilise 140 DJ base 8.3 C depuis le 1er avril + Tair max 16.5 C.
  # Aucun equivalent publie pour le stenographe : DJ_VOL est le parametre a
  # calibrer en priorite (voir calibrer_seuil_vol()).
  DOY_DEB_CUMUL = 1,     # [EXPERT] debut du cumul prevol (1 = 1er janvier)
  DJ_VOL       = 120,    # [EXPERT] somme thermique prevol (base T_MIN_DEV)
  SEUIL_VOL    = 20.0,   # [PUBLIE/qualitatif] Tair max declenchant le vol (C)
  SEUIL_GEL    = 0.0,    # [EXPERT] Tair min en dessous de laquelle le vol est nul
  DOY_MIN_VOL  = 46,     # [EXPERT] garde-fou : pas d'essaimage avant le 15/02
  
  # --- Arret d'initiation des generations (photoperiode) ---
  # Ephytia : vols en continu jusqu'a octobre environ. 12 h de jour ~ 1er octobre
  # a 45 N. Au-dela, la descendance hiverne (larve/nymphe/jeune adulte).
  DL_INIT_MIN  = 12.0,   # [EXPERT] duree du jour mini pour initier une generation
  
  # --- Divers ---
  N_GEN_MAX    = 4       # nb max de generations suivies dans les sorties
)

# Coefficient d'amortissement du terme radiatif de l'equation de Baier (2007),
# pour tenir compte de l'ecorce epaisse et fissuree du bas de tronc des pins,
# ou se developpe le stenographe. Voir le commentaire detaille de temp_ecorce().
#   1.00 = equation de Baier telle quelle (epicea, ecorce fine)  -> borne haute
#   0.50 = hypothese de travail par defaut                       -> ~+2 a +4 C
#   0.20 = proche des mesures de phloeme de pin (+1 a +2 C)      -> borne basse
#   0.00 = equivalent a exposition="shaded"
K_ECORCE <- 0.2   # [EXPERT] a calibrer sur sondes sous ecorce

cat("=== STENO-DJ — Ips sexdentatus — France metropolitaine (SAFRAN) v1 ===\n")
cat(glue("  Periode    : {DATE_DEB} -> {DATE_FIN_SAFRAN}\n"))
cat(glue("  Modele     : degres-jours | seuil {STENO$T_MIN_DEV} C | ",
         "{STENO$DJ_GEN} DJ/generation\n"))
cat(glue("  Resolution : {SAFRAN_RES/1000} km | Dalles : {DALLE_KM}x{DALLE_KM} km\n"))
cat(glue("  Cache      : {DIR_CACHE}/\n\n"))


# =============================================================================
# 2. GRILLE DE DALLES FRANCE   (identique PHENIPS)
# =============================================================================

# =============================================================================
# 3. RASTER DE REFERENCE FRANCE (grille SAFRAN native, L93)
# =============================================================================

RST_FRANCE <- rast(
  xmin = FRANCE_XMIN, xmax = FRANCE_XMAX,
  ymin = FRANCE_YMIN, ymax = FRANCE_YMAX,
  resolution = SAFRAN_RES, crs = CRS_L93
)

cat(glue("  Raster France : {nrow(RST_FRANCE)}L x {ncol(RST_FRANCE)}C",
         " = {ncell(RST_FRANCE)} cellules\n"))
cat(glue("  Dalles France : {nrow(DALLES)}\n\n"))


# =============================================================================
# 4. TELECHARGEMENT SAFRAN PAR DALLE (avec cache)   (identique PHENIPS)
# =============================================================================

verifier_dalles <- function(dalles = DALLES, reparer = FALSE) {
  cat(">> Controle qualite des caches...\n")
  suspects <- character(0)
  
  for (i in seq_len(nrow(dalles))) {
    f <- chemin_cache(dalles$id[i])
    if (!file.exists(f)) next
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    
    nx <- grep("^x$", tolower(names(dt)))
    ny <- grep("^y$", tolower(names(dt)))
    nt <- grep("^time$|^date$", tolower(names(dt)))
    if (length(nx) == 0 || length(ny) == 0 || length(nt) == 0) next
    
    setnames(dt, c(names(dt)[nx[1]], names(dt)[ny[1]], names(dt)[nt[1]]),
             c(".x", ".y", ".t"))
    par_pt <- dt[, .N, by = .(.x, .y)]
    n_dates <- uniqueN(dt$.t)
    
    if (nrow(par_pt) > 1 && min(par_pt$N) < n_dates) {
      n_incomplets <- sum(par_pt$N < n_dates)
      cat(glue("  [TRONQUE] dalle {dalles$id[i]} : {n_incomplets}/{nrow(par_pt)} ",
               "points incomplets ({min(par_pt$N)}/{n_dates} dates)\n"))
      suspects <- c(suspects, f)
    }
  }
  
  if (length(suspects) == 0) {
    cat("  Aucune dalle tronquee detectee.\n\n")
  } else if (reparer) {
    file.remove(suspects)
    cat(glue("  {length(suspects)} cache(s) supprime(s) -> relancer recuperer_safran()\n\n"))
  } else {
    cat(glue("  {length(suspects)} cache(s) suspect(s). Relancer avec reparer = TRUE.\n\n"))
  }
  invisible(suspects)
}

#' Detecte les trous internes d'un raster (colonnes/lignes manquantes au milieu
#' de l'emprise valide) — signature typique d'une dalle tronquee.
verifier_trous <- function(rst, seuil = 5) {
  m <- as.matrix(rst, wide = TRUE)
  trous <- 0L
  for (j in seq_len(ncol(m))) {
    ok <- which(!is.na(m[, j]))
    if (length(ok) < 2) next
    n <- sum(is.na(m[min(ok):max(ok), j]))
    if (n >= seuil) {
      trous <- trous + n
      cat(glue("  [TROU] colonne {j} (x~{round(xFromCol(rst, j))}) : {n} cellules\n"))
    }
  }
  if (trous == 0) cat("  Aucun trou interne detecte.\n")
  invisible(trous)
}


# =============================================================================
# 5. LECTURE ET ASSEMBLAGE RASTERS NATIONAUX
# =============================================================================

# =============================================================================
# 6. DUREE DU JOUR
# =============================================================================

# =============================================================================
# 7. TEMPERATURE SOUS ECORCE ET DEGRES-JOURS
# =============================================================================

#' Temperature sous ecorce — equation de Baier et al. (2007), utilisee par PHENIPS
#'
#' BT = -0.173 + K_ECORCE * 0.0008518 * I + 1.054 * Tair
#'   I    : rayonnement global journalier (Wh/m2)
#'   Tair : temperature de l'air (C)
#'
#' -------------------------------------------------------------------------
#' POURQUOI UN COEFFICIENT K_ECORCE POUR LE PIN
#' -------------------------------------------------------------------------
#' L'equation de Baier est calibree sur epicea (Picea abies), sur des zones de
#' tronc a ecorce relativement fine — l'habitat d'Ips typographus. Elle donne
#' un echauffement de +3 a +8 C en ete.
#'
#' Le stenographe occupe un compartiment tres different : il se developpe dans
#' la MOITIE INFERIEURE du tronc, dans les ecorces epaisses et fissurees des
#' deux premiers metres (Nageleisen/DSF ; Chararas). C'est precisement le
#' compartiment le mieux isole thermiquement — a l'oppose d'Ips acuminatus qui
#' colonise les ecorces fines du haut de tige, ou l'effet du rayonnement est
#' maximal.
#'
#' Ordres de grandeur mesures dans la litterature :
#'   - Temperature de SURFACE d'ecorce de pin exposee au sud : jusqu'a +28 a
#'     +37 C au-dessus de l'air (Stoutjesdijk, faible angle solaire hivernal).
#'     >> non pertinent : c'est la surface, pas le phloeme.
#'   - Temperature de PHLOEME (sous ecorce) vs air : +1 a +2 C seulement
#'     (Bolstad et al. 1997 ; Ungerer et al. 1999, pins nord-americains).
#'   - Baier et al. 2007 (epicea, sunny) : +3 a +8 C.
#'
#' Il n'existe aucune equation publiee pour le phloeme de pin sous ecorce
#' epaisse. K_ECORCE = 0.5 est donc une HYPOTHESE DE TRAVAIL qui place
#' l'echauffement estival vers +2 a +4 C, entre les mesures nord-americaines
#' et Baier. A calibrer avec des sondes sous ecorce (voir note en fin de
#' fichier).
#'
#' Nuance importante : le stenographe SELECTIONNE activement les faces chaudes
#' (face superieure des troncs abattus, faces SE plutot que NO — Jactel &
#' Lieutier). Un scenario entierement "shaded" sous-estime donc la realite.
#' D'ou l'encadrement systematique sunny / shaded.
temp_ecorce <- function(tair, rad_whm2,
                        exposition = c("sunny", "shaded"),
                        k_ecorce   = K_ECORCE) {
  exposition <- match.arg(exposition)
  if (exposition == "shaded") return(tair)
  rad_whm2[is.na(rad_whm2)] <- 0
  -0.173 + k_ecorce * 0.0008518 * rad_whm2 + 1.054 * tair
}

#' Degres-jours journaliers — methode sinus simple (Baskerville & Emin 1969)
#' avec seuil inferieur TL et cutoff horizontal superieur TU.
#' Vectorise sur des vecteurs tmin/tmax de meme longueur.
dj_sinus <- function(tmin, tmax, TL, TU) {
  n  <- length(tmin)
  dj <- rep(NA_real_, n)
  ok <- !is.na(tmin) & !is.na(tmax)
  if (!any(ok)) return(dj)
  
  a <- tmin[ok]; b <- tmax[ok]
  # securite : inversion eventuelle
  inv <- b < a
  if (any(inv)) { tmp <- a[inv]; a[inv] <- b[inv]; b[inv] <- tmp }
  
  # cutoff horizontal : tout ce qui depasse TU compte comme TU
  a <- pmin(a, TU); b <- pmin(b, TU)
  
  moy <- (a + b) / 2
  amp <- (b - a) / 2
  
  d <- numeric(length(a))
  
  # cas 1 : journee entierement sous le seuil
  c_bas <- b <= TL
  d[c_bas] <- 0
  
  # cas 2 : journee entierement au-dessus du seuil (ou amplitude nulle)
  c_haut <- !c_bas & (a >= TL | amp <= 0)
  d[c_haut] <- pmax(moy[c_haut] - TL, 0)
  
  # cas 3 : le seuil coupe la sinusoide
  c_mix <- !c_bas & !c_haut
  if (any(c_mix)) {
    theta <- asin(pmin(pmax((TL - moy[c_mix]) / amp[c_mix], -1), 1))
    d[c_mix] <- ((moy[c_mix] - TL) * (pi/2 - theta) +
                   amp[c_mix] * cos(theta)) / pi
  }
  
  dj[ok] <- pmax(d, 0)
  dj
}


# =============================================================================
# 8. MODELE STENO-DJ
# =============================================================================

#' Modele phenologique degres-jours du stenographe
#'
#' Automate a etats parcourant la saison jour par jour, vectorise sur toutes
#' les cellules SAFRAN simultanement.
#'
#' Regles :
#'   1. Essaimage initial : N_JOURS_VOL jours (non necessairement consecutifs)
#'      avec Tair_max >= SEUIL_VOL et Tair_min > SEUIL_GEL, a partir de DOY_MIN_VOL.
#'   2. Developpement : cumul de DJ (sinus simple, base T_MIN_DEV, cutoff
#'      T_MAX_DEV) calcules sur la temperature SOUS ECORCE.
#'   3. Generation achevee a DJ_GEN + DJ_MATURATION ; le reliquat de DJ est
#'      reporte sur la generation suivante.
#'   4. Une nouvelle generation n'est initiee que si la duree du jour est
#'      >= DL_INIT_MIN ; sinon la descendance est consideree comme hivernante
#'      et le cumul s'arrete pour l'annee.
#'
#' @param meteo liste de SpatRasters : tmean, tmin, tmax, rad, daylength
#' @param exposition "sunny" (defaut) ou "shaded"
#' @param par liste de parametres biologiques (defaut : STENO)
#' @return liste de SpatRasters + metadonnees
modeliser_steno <- function(meteo,
                            exposition = "sunny",
                            k_ecorce   = K_ECORCE,
                            par        = STENO) {
  cat(glue(">> Modelisation STENO-DJ [exposition={exposition}]...\n"))
  
  dates <- as.Date(time(meteo$tmean))
  n_j   <- length(dates)
  n_c   <- ncell(meteo$tmean)
  cat(glue("   {n_c} cellules x {n_j} jours\n"))
  
  # --- Extraction en matrices (cellules x jours) ---
  m_tmin <- values(meteo$tmin)
  m_tmax <- values(meteo$tmax)
  m_rad  <- values(meteo$rad)
  m_dl   <- values(meteo$daylength)
  
  # Temperature sous ecorce (min nocturne : pas de rayonnement -> tair)
  bt_min <- if (exposition == "sunny") 1.054 * m_tmin - 0.173 else m_tmin
  bt_max <- temp_ecorce(m_tmax, m_rad, exposition, k_ecorce)
  cat(glue("   K_ECORCE = {k_ecorce} | echauffement median : ",
           "+{round(k_ecorce * 0.0008518 * median(m_rad, na.rm = TRUE), 2)} C\n"))
  
  # Cellules exploitables (au moins une donnee valide)
  valide <- rowSums(!is.na(m_tmax)) > 0
  
  # --- Etats ---
  essaime   <- rep(FALSE, n_c)
  dj_prevol <- rep(0,     n_c)   # somme thermique accumulee avant l'essaimage
  doy_onset <- rep(NA_integer_, n_c)
  dev_cur   <- rep(0,     n_c)   # avancement (0-1) de la generation en cours
  dj_cum    <- rep(0,     n_c)   # DJ cumules depuis l'essaimage
  n_gen     <- rep(0L,    n_c)
  hiverne   <- rep(FALSE, n_c)   # plus d'initiation possible cette annee
  doy_gen   <- matrix(NA_integer_, nrow = n_c, ncol = par$N_GEN_MAX)
  
  DJ_CYCLE <- par$DJ_GEN + par$DJ_MATURATION
  
  # --- Boucle temporelle ---
  for (j in seq_len(n_j)) {
    if (j %% 30 == 0) cat(glue("\r   Jour {j}/{n_j}  "))
    doy <- yday(dates[j])
    
    tmx <- m_tmax[, j]; tmn <- m_tmin[, j]; dl <- m_dl[, j]
    
    # Degres-jours du jour (sous ecorce) — sert au cumul prevol ET au cycle
    dj <- dj_sinus(bt_min[, j], bt_max[, j], par$T_MIN_DEV, par$T_MAX_DEV)
    dj[is.na(dj)] <- 0
    
    # --- 1. Essaimage initial : somme thermique + journee volable ---
    if (doy >= par$DOY_DEB_CUMUL) {
      pre <- valide & !essaime
      dj_prevol[pre] <- dj_prevol[pre] + dj[pre]
    }
    
    vol_ok <- !is.na(tmx) & !is.na(tmn) &
      tmx >= par$SEUIL_VOL & tmn > par$SEUIL_GEL
    
    nouveau <- valide & !essaime & vol_ok &
      dj_prevol >= par$DJ_VOL & doy >= par$DOY_MIN_VOL
    if (any(nouveau)) {
      essaime[nouveau]   <- TRUE
      doy_onset[nouveau] <- doy
    }
    
    # --- 2. Developpement ---
    actif <- essaime & !hiverne
    if (!any(actif)) next
    
    dj_cum[actif]  <- dj_cum[actif]  + dj[actif]
    dev_cur[actif] <- dev_cur[actif] + dj[actif] / DJ_CYCLE
    
    # --- 3. Generation achevee ---
    fini <- actif & dev_cur >= 1
    if (any(fini)) {
      # une seule generation par jour au maximum (DJ journalier << 517)
      n_gen[fini] <- n_gen[fini] + 1L
      idx <- which(fini)
      rang <- pmin(n_gen[idx], par$N_GEN_MAX)
      doy_gen[cbind(idx, rang)] <- doy
      
      # 4. Nouvelle generation initiee seulement si photoperiode suffisante
      peut <- fini & !is.na(dl) & dl >= par$DL_INIT_MIN
      dev_cur[peut] <- dev_cur[peut] - 1
      bloq <- fini & !peut
      dev_cur[bloq] <- 1
      hiverne[bloq] <- TRUE
    }
  }
  cat("\n")
  
  # --- Rasterisation des sorties ---
  faire_rst <- function(v, nom) {
    r <- RST_FRANCE
    v[!valide] <- NA
    values(r) <- v
    names(r)  <- nom
    r
  }
  
  res <- list(
    dates      = dates,
    exposition = exposition,
    k_ecorce   = k_ecorce,
    par        = par,
    valide     = valide,
    onset      = faire_rst(as.numeric(doy_onset), "doy_essaimage"),
    dj_cum     = faire_rst(dj_cum,                "dj_cumules"),
    n_gen      = faire_rst(as.numeric(n_gen),     "generations_achevees"),
    dev_cur    = faire_rst(dev_cur * 100,         "dev_gen_en_cours_pct"),
    # indicateur continu : le plus lisible en carte
    gen_cont   = faire_rst(ifelse(essaime, n_gen + dev_cur, NA_real_),
                           "generations_cumulees"),
    hiverne    = faire_rst(as.numeric(hiverne),   "descendance_hivernante")
  )
  
  for (g in seq_len(par$N_GEN_MAX)) {
    res[[as.character(glue("doy_gen{g}"))]] <-
      faire_rst(as.numeric(doy_gen[, g]), as.character(glue("doy_fin_gen{g}")))
  }
  
  cat(glue("  Modelisation terminee ({n_j} jours).\n\n"))
  res
}


# =============================================================================
# 8bis. CALIBRATION DU SEUIL D'ESSAIMAGE
# =============================================================================

#' Explore la sensibilite de la date de premier envol au parametre DJ_VOL
#'
#' DJ_VOL est le parametre le moins contraint du modele. Cette fonction rejoue
#' uniquement le bloc "essaimage" (pas tout le cycle) pour plusieurs valeurs de
#' seuil et renvoie la distribution des dates obtenues. A confronter aux dates
#' de premieres captures en pieges a pheromone (ipsdienol / cis-verbenol).
#'
#' Un bon seuil doit produire un ETALEMENT spatial : si les quartiles D1 et D9
#' sont separes de moins de trois semaines, le critere se comporte encore en
#' fonction en escalier et ne discrimine rien.
#'
#' @param meteo  liste de SpatRasters (sortie d'assembler_rasters_nationaux + daylength)
#' @param seuils vecteur de valeurs de DJ_VOL a tester
calibrer_seuil_vol <- function(meteo,
                               seuils     = c(40, 60, 80, 100, 120, 160, 200, 250),
                               exposition = "sunny",
                               k_ecorce   = K_ECORCE,
                               par        = STENO) {
  dates <- as.Date(time(meteo$tmean))
  n_j   <- length(dates)
  n_c   <- ncell(meteo$tmean)
  
  m_tmin <- values(meteo$tmin); m_tmax <- values(meteo$tmax)
  m_rad  <- values(meteo$rad)
  bt_min <- if (exposition == "sunny") 1.054 * m_tmin - 0.173 else m_tmin
  bt_max <- temp_ecorce(m_tmax, m_rad, exposition, k_ecorce)
  valide <- rowSums(!is.na(m_tmax)) > 0
  
  # Cumul thermique et jours volables, calcules une seule fois
  dj_pre <- matrix(0, nrow = n_c, ncol = n_j)
  vol    <- matrix(FALSE, nrow = n_c, ncol = n_j)
  acc    <- rep(0, n_c)
  for (j in seq_len(n_j)) {
    d <- dj_sinus(bt_min[, j], bt_max[, j], par$T_MIN_DEV, par$T_MAX_DEV)
    d[is.na(d)] <- 0
    if (yday(dates[j]) >= par$DOY_DEB_CUMUL) acc <- acc + d
    dj_pre[, j] <- acc
    vol[, j] <- !is.na(m_tmax[, j]) & !is.na(m_tmin[, j]) &
      m_tmax[, j] >= par$SEUIL_VOL & m_tmin[, j] > par$SEUIL_GEL &
      yday(dates[j]) >= par$DOY_MIN_VOL
  }
  
  res <- lapply(seuils, function(sv) {
    ok  <- vol & (dj_pre >= sv)
    idx <- max.col(ok, ties.method = "first")
    idx[rowSums(ok) == 0 | !valide] <- NA
    doy <- yday(dates)[idx]
    q <- quantile(doy, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
    data.frame(
      dj_vol     = sv,
      pct_essaim = round(100 * sum(!is.na(doy)) / sum(valide), 1),
      d1         = format(doy_vers_date(q[1]), "%d/%m"),
      med        = format(doy_vers_date(q[3]), "%d/%m"),
      d9         = format(doy_vers_date(q[5]), "%d/%m"),
      etalement_j = as.integer(q[5] - q[1])
    )
  })
  out <- do.call(rbind, res)
  print(out, row.names = FALSE)
  invisible(out)
}


# =============================================================================
# 9. EXPORT RASTERS
# =============================================================================

exporter_rasters <- function(steno, dir_out = DIR_OUT) {
  cat(">> Export GeoTIFF (L93, LZW)...\n")
  date_j <- format(max(steno$dates), "%Y%m%d")
  
  ecrire <- function(rst, nom) {
    ch <- file.path(dir_out, nom)
    writeRaster(rst, ch, overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES"))
    cat(glue("  {nom}\n"))
  }
  
  ecrire(steno$onset,    glue("steno_essaimage_{ANNEE}_L93.tif"))
  ecrire(steno$dj_cum,   glue("steno_dj_{date_j}_L93.tif"))
  ecrire(steno$gen_cont, glue("steno_generations_{date_j}_L93.tif"))
  ecrire(steno$n_gen,    glue("steno_ngen_entier_{date_j}_L93.tif"))
  ecrire(steno$dev_cur,  glue("steno_dev_encours_{date_j}_L93.tif"))
  
  cat(glue("\n  Rasters exportes dans {dir_out}/\n\n"))
}


# =============================================================================
# 9bis. EXPORT GEOJSON (nb de generations + premier envol)
# =============================================================================
# GeoJSON = WGS84 obligatoire (RFC 7946). On reprojette depuis le L93.
# Deux formes de sortie, complementaires :
#   - "zones"   : polygones fusionnes par classe (leger, lisible en carto)
#   - "grille"  : une entite par maille SAFRAN 8 km, tous attributs
#                 (plus lourd mais permet le requetage et les jointures)

#' Vectorise un raster classe en polygones fusionnes
#' @param rst raster de classes entieres (NA = hors emprise)
#' @param nom_champ nom de l'attribut de classe dans la sortie
polygoniser_classes <- function(rst, nom_champ) {
  v <- as.polygons(rst, dissolve = TRUE, na.rm = TRUE)
  if (nrow(v) == 0) return(NULL)
  names(v) <- nom_champ
  sf_obj <- st_as_sf(v)
  sf_obj[[nom_champ]] <- as.integer(sf_obj[[nom_champ]])
  # Surface calculee en L93 (equivalent-surface) AVANT reprojection
  sf_obj$surface_km2 <- round(as.numeric(st_area(st_geometry(sf_obj))) / 1e6, 1)
  st_transform(sf_obj, 4326)
}

#' Ecrit un objet sf en GeoJSON (WGS84, precision 5 decimales ~1 m)
ecrire_geojson <- function(sf_obj, chemin) {
  if (is.null(sf_obj) || nrow(sf_obj) == 0) {
    cat(glue("  [SKIP] {basename(chemin)} : aucune entite\n")); return(invisible(NULL))
  }
  if (file.exists(chemin)) file.remove(chemin)   # GDAL ne gere pas l'overwrite ici
  st_write(sf_obj, chemin, driver = "GeoJSON", quiet = TRUE,
           layer_options = c("RFC7946=YES", "COORDINATE_PRECISION=5"))
  cat(glue("  {basename(chemin)} ({nrow(sf_obj)} entites)\n"))
  invisible(chemin)
}

#' Convertit un jour de l'annee en date lisible
doy_vers_date <- function(doy, annee = ANNEE) {
  as.Date(doy - 1, origin = as.Date(glue("{annee}-01-01")))
}

#' Export GeoJSON : nb de generations achevees + date du premier envol
#'
#' @param steno    sortie de modeliser_steno()
#' @param pas_essaimage largeur des classes de date d'essaimage, en jours
#' @param grille   TRUE = exporte aussi la grille SAFRAN complete (1 entite/maille)
exporter_geojson <- function(steno,
                             dir_out        = DIR_OUT,
                             pas_essaimage  = 15,
                             grille         = TRUE) {
  cat(">> Export GeoJSON (WGS84)...\n")
  date_j <- format(max(steno$dates), "%Y%m%d")
  
  # ---------------------------------------------------------------------------
  # 1. Nombre de generations achevees — zones fusionnees
  # ---------------------------------------------------------------------------
  zones_gen <- polygoniser_classes(steno$n_gen, "nb_generations")
  
  if (!is.null(zones_gen)) {
    zones_gen$libelle <- ifelse(
      zones_gen$nb_generations == 0,
      "Aucune generation achevee",
      as.character(glue("{zones_gen$nb_generations} generation(s) achevee(s)"))
    )
    # Palette sequentielle viridis, cohérente avec PHENIPS/CHAPY : 0 = gris
    # (aucune generation), puis jaune -> vert -> bleu -> violet a mesure que
    # le nombre de generations achevees augmente.
    palette_gen <- c("#E8E8E8", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF")
    idx_couleur <- pmin(zones_gen$nb_generations + 1L, length(palette_gen))
    zones_gen$couleur <- palette_gen[idx_couleur]
    zones_gen$date_situation <- format(max(steno$dates), "%Y-%m-%d")
    zones_gen$modele         <- "STENO-DJ v1"
    zones_gen$exposition     <- steno$exposition
    zones_gen$k_ecorce       <- steno$k_ecorce
    zones_gen <- zones_gen[order(zones_gen$nb_generations), ]
    
    ecrire_geojson(zones_gen,
                   file.path(dir_out, glue("steno_generations_{date_j}.geojson")))
  }
  
  # ---------------------------------------------------------------------------
  # 2. Date du premier envol — zones fusionnees par classe de N jours
  # ---------------------------------------------------------------------------
  ons <- steno$onset
  vals <- values(ons)[, 1]
  
  if (all(is.na(vals))) {
    cat("  [SKIP] essaimage : aucun envol modelise a ce stade\n")
  } else {
    dmin <- floor(min(vals, na.rm = TRUE) / pas_essaimage) * pas_essaimage
    dmax <- ceiling(max(vals, na.rm = TRUE) / pas_essaimage) * pas_essaimage
    brk  <- seq(dmin, dmax, by = pas_essaimage)
    if (length(brk) < 2) brk <- c(dmin, dmin + pas_essaimage)
    
    # classify : intervalles [from, to), classe = borne inferieure (doy)
    rcl <- cbind(brk[-length(brk)], brk[-1], brk[-length(brk)])
    ons_cl <- classify(ons, rcl, include.lowest = TRUE, right = FALSE)
    
    zones_ess <- polygoniser_classes(ons_cl, "doy_debut")
    
    if (!is.null(zones_ess)) {
      zones_ess$doy_fin    <- zones_ess$doy_debut + pas_essaimage - 1L
      zones_ess$date_debut <- format(doy_vers_date(zones_ess$doy_debut), "%Y-%m-%d")
      zones_ess$date_fin   <- format(doy_vers_date(zones_ess$doy_fin),   "%Y-%m-%d")
      zones_ess$libelle    <- as.character(glue(
        "Premier envol : {format(doy_vers_date(zones_ess$doy_debut), '%d/%m')} - ",
        "{format(doy_vers_date(zones_ess$doy_fin), '%d/%m')}"
      ))
      # Degrade jaune -> rouge fonce selon la precocite (meme logique que
      # l'onset PHENIPS) : precalcule ici plutot que cote app, pour que
      # QGIS/toute autre visu affiche la meme chose que l'app Shiny.
      pal_essaimage <- scales::col_numeric("YlOrRd", domain = range(zones_ess$doy_debut))
      zones_ess$couleur     <- pal_essaimage(zones_ess$doy_debut)
      zones_ess$annee       <- ANNEE
      zones_ess$modele      <- "STENO-DJ v1"
      zones_ess$seuil_vol_c <- steno$par$SEUIL_VOL
      zones_ess$dj_vol      <- steno$par$DJ_VOL
      zones_ess <- zones_ess[order(zones_ess$doy_debut), ]
      
      ecrire_geojson(zones_ess,
                     file.path(dir_out, glue("steno_essaimage_{ANNEE}.geojson")))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3. Grille SAFRAN complete (optionnel) — tous attributs, 1 entite par maille
  # ---------------------------------------------------------------------------
  if (grille) {
    pts <- as.points(steno$n_gen, na.rm = TRUE)
    if (nrow(pts) == 0) {
      cat("  [SKIP] grille : aucune maille valide\n")
    } else {
      # Polygones = mailles SAFRAN carrees de 8 km (construites explicitement :
      # st_buffer(endCapStyle="SQUARE") produit des cercles, ne pas l'utiliser)
      xy <- crds(pts)
      d  <- SAFRAN_RES / 2
      geom <- lapply(seq_len(nrow(xy)), function(i) {
        st_polygon(list(cbind(
          c(xy[i,1]-d, xy[i,1]+d, xy[i,1]+d, xy[i,1]-d, xy[i,1]-d),
          c(xy[i,2]-d, xy[i,2]-d, xy[i,2]+d, xy[i,2]+d, xy[i,2]-d)
        )))
      })
      
      ext_val <- function(rst) as.numeric(extract(rst, xy)[, 1])
      
      g <- st_sf(
        doy_essaimage   = as.integer(ext_val(steno$onset)),
        nb_generations  = as.integer(ext_val(steno$n_gen)),
        dj_cumules      = round(ext_val(steno$dj_cum), 1),
        dev_encours_pct = round(ext_val(steno$dev_cur), 1),
        gen_cumulees    = round(ext_val(steno$gen_cont), 2),
        hivernante      = as.integer(ext_val(steno$hiverne)),
        geometry        = st_sfc(geom, crs = 2154)
      )
      g$date_essaimage <- ifelse(
        is.na(g$doy_essaimage), NA_character_,
        format(doy_vers_date(g$doy_essaimage), "%Y-%m-%d")
      )
      g$date_situation <- format(max(steno$dates), "%Y-%m-%d")
      g <- st_transform(g, 4326)
      
      ecrire_geojson(g, file.path(dir_out,
                                  glue("steno_grille_safran_{date_j}.geojson")))
    }
  }
  
  cat(glue("  GeoJSON exportes dans {dir_out}/\n\n"))
}


# =============================================================================
# 10. CARTES NATIONALES
# =============================================================================

carte_france <- function(rst_l93, titre, legende, palette,
                         breaks = NULL, date_label = NULL, fichier = NULL) {
  if (is.null(date_label)) date_label <- format(DATE_FIN_SAFRAN, "%d/%m/%Y")
  rst_wgs <- project(rst_l93, CRS_WGS)
  
  p <- ggplot() +
    geom_spatraster(data = rst_wgs) +
    scale_fill_gradientn(
      colours  = palette,
      values   = if (!is.null(breaks)) rescale(breaks) else NULL,
      na.value = "grey92",
      name     = legende,
      guide    = guide_colorbar(barwidth = 14, barheight = 0.8,
                                title.position = "top")
    ) +
    labs(
      title    = titre,
      subtitle = glue("STENO-DJ (Pineau et al. 2017) | SAFRAN INRAE | {date_label} | {ANNEE}"),
      caption  = "Indicateur d'aide a la decision — non valide terrain | GeoSAS INRAE — CNPF/CRPF AuRA",
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

generer_cartes <- function(steno, dir_out = DIR_OUT) {
  cat(">> Generation des cartes nationales...\n")
  date_j     <- format(max(steno$dates), "%Y%m%d")
  date_label <- format(max(steno$dates), "%d/%m/%Y")
  
  # 1. Date du premier essaimage
  carte_france(
    steno$onset,
    titre   = "Date du premier essaimage — Ips sexdentatus — France",
    legende = glue("Jour de l'annee (1 = 01/01/{ANNEE})"),
    palette = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF", "#FFF3B0"),
    date_label = date_label,
    fichier = file.path(dir_out, glue("steno_essaimage_{date_j}.png"))
  )
  
  # 2. Generations cumulees (indicateur principal)
  carte_france(
    steno$gen_cont,
    titre   = glue("Generations cumulees — Ips sexdentatus — {date_label}"),
    legende = "Generations (entier = generation achevee)",
    palette = c("#F7FCF0", "#BAE4BC", "#7BCCC4", "#F39C12", "#C0392B", "#6E0000"),
    breaks  = c(0, 1, 2, 3, 4, 5),
    date_label = date_label,
    fichier = file.path(dir_out, glue("steno_generations_{date_j}.png"))
  )
  
  # 3. Degres-jours cumules
  carte_france(
    steno$dj_cum,
    titre   = glue("Degres-jours cumules (base {STENO$T_MIN_DEV} C, sous ecorce) — {date_label}"),
    legende = "DJ cumules depuis l'essaimage",
    palette = c("grey95", "#FFF3B0", "#35B779FF", "#31688EFF", "#440154FF"),
    date_label = date_label,
    fichier = file.path(dir_out, glue("steno_dj_{date_j}.png"))
  )
  
  # 4. Avancement de la generation en cours
  carte_france(
    steno$dev_cur,
    titre   = glue("Developpement de la generation en cours — {date_label}"),
    legende = "Avancement (%)",
    palette = c("grey95", "#FFF3B0", "#35B779FF", "#31688EFF", "#440154FF"),
    breaks  = c(0, 25, 50, 75, 100),
    date_label = date_label,
    fichier = file.path(dir_out, glue("steno_dev_{date_j}.png"))
  )
  
  cat(glue("  Cartes exportees dans {dir_out}/\n\n"))
}


# =============================================================================
# 11. SYNTHESE CONSOLE
# =============================================================================

synthese_console <- function(steno) {
  cat("\n", strrep("=", 62), "\n")
  cat("  SYNTHESE STENO-DJ — Ips sexdentatus — France metropolitaine\n")
  cat(strrep("=", 62), "\n")
  
  n_valide <- sum(steno$valide)
  cat(glue("  Date        : {format(max(steno$dates), '%d/%m/%Y')}\n"))
  cat(glue("  Exposition  : {steno$exposition} | K_ECORCE = {steno$k_ecorce}\n"))
  cat(glue("  Cellules    : {n_valide} valides / {ncell(RST_FRANCE)} total\n\n"))
  
  ons <- values(steno$onset)[, 1]
  n_e <- sum(!is.na(ons))
  cat("  --- Essaimage ---\n")
  cat(glue("  Cellules ayant essaime : {n_e} ({round(100*n_e/n_valide,1)}%)\n"))
  if (n_e > 0) {
    q <- quantile(ons, c(0.1, 0.5, 0.9), na.rm = TRUE)
    d <- as.Date(q - 1, origin = glue("{ANNEE}-01-01"))
    cat(glue("  Date d'essaimage (D1 / mediane / D9) : ",
             "{format(d[1],'%d/%m')} / {format(d[2],'%d/%m')} / {format(d[3],'%d/%m')}\n\n"))
  } else cat("\n")
  
  dj <- values(steno$dj_cum)[, 1]
  cat("  --- Degres-jours cumules (sous ecorce) ---\n")
  cat(glue("  Moyenne : {round(mean(dj, na.rm=TRUE))} DJ | ",
           "Max : {round(max(dj, na.rm=TRUE))} DJ\n"))
  cat(glue("  Equivalent : {round(mean(dj, na.rm=TRUE)/steno$par$DJ_GEN, 2)} generation(s) en moyenne\n\n"))
  
  ng <- values(steno$n_gen)[, 1]
  cat("  --- Generations achevees ---\n")
  for (g in 0:steno$par$N_GEN_MAX) {
    n <- sum(ng == g, na.rm = TRUE)
    if (n > 0) cat(glue("  {g} generation(s) : {n} cellules ({round(100*n/n_valide,1)}%)\n"))
  }
  cat("\n")
  
  hv <- values(steno$hiverne)[, 1]
  n_h <- sum(hv == 1, na.rm = TRUE)
  if (n_h > 0)
    cat(glue("  Descendance hivernante (photoperiode) : {n_h} cellules\n\n"))
  
  cat("  /!\\ Indicateur non valide terrain — a croiser avec un piegeage\n")
  cat("      pheromonal local (ipsdienol / cis-verbenol).\n")
  cat(strrep("=", 62), "\n\n")
}


# =============================================================================
# 12. DIAGNOSTIC RAYONNEMENT (a lancer une fois sur ton cache)
# =============================================================================

#' Verifie l'unite reelle de SSI_Q dans le cache et l'impact sur la temperature
#' sous ecorce. A lancer AUSSI pour auditer le script PHENIPS v3.
diagnostiquer_rayonnement <- function(dalles = DALLES[1:min(5, nrow(DALLES)), ]) {
  fichiers <- chemin_cache(dalles$id)
  fichiers <- fichiers[file.exists(fichiers)]
  if (length(fichiers) == 0) { cat("Aucun cache disponible.\n"); return(invisible(NULL)) }
  
  dt <- rbindlist(lapply(fichiers, function(f) {
    x <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    col <- names(x)[grep("ssi", tolower(names(x)))]
    if (length(col) == 0) return(NULL)
    data.table(ssi = x[[col[1]]])
  }), fill = TRUE)
  
  if (nrow(dt) == 0) { cat("SSI_Q introuvable dans le cache.\n"); return(invisible(NULL)) }
  
  v   <- dt$ssi[!is.na(dt$ssi) & dt$ssi > 0]
  med <- median(v)
  conv <- detecter_facteur_rad(med)
  
  cat("\n--- Diagnostic rayonnement SSI_Q ---\n")
  cat(glue("  n = {length(v)} | min={round(min(v),2)} | median={round(med,2)} | max={round(max(v),2)}\n"))
  cat(glue("  Unite deduite : {conv$unite} -> facteur {round(conv$facteur,4)}\n"))
  cat(glue("  Mediane convertie : {round(med * conv$facteur)} Wh/m2\n"))
  cat(glue("  Effet sur la temperature d'ecorce (0.0008518 * I) : ",
           "+{round(0.0008518 * med * conv$facteur, 2)} C\n"))
  cat("  Reference PHENIPS : l'echauffement attendu est de +3 a +8 C en ete.\n")
  cat("  Si tu obtiens < 0.5 C, le facteur de conversion est trop petit\n")
  cat("  (c'est le cas si phenips_france_safran_v3.R utilise WHm2_FACTOR=2.7778\n")
  cat("   sur des valeurs deja exprimees en MJ/m2).\n---\n\n")
  invisible(conv)
}


# =============================================================================
# 13. PIPELINE PRINCIPAL
# =============================================================================

#' Pipeline complet France entiere
#' @param exposition "sunny" (defaut) ou "shaded" (scenario prudent)
#' @param pas_essaimage largeur (jours) des classes de date d'envol en GeoJSON
#' @param grille    TRUE = exporte aussi la grille SAFRAN complete en GeoJSON
#' @param force_dl   TRUE = purge le cache et retelecharge
pipeline_steno <- function(exposition   = "sunny",
                           k_ecorce     = K_ECORCE,
                           pas_essaimage = 15,
                           grille       = TRUE,
                           force_dl     = FALSE) {
  cat("\n", strrep("=", 62), "\n")
  cat("  STENO-DJ — Ips sexdentatus — France metropolitaine — SAFRAN v1\n")
  cat(strrep("=", 62), "\n\n")
  t0 <- proc.time()
  
  if (force_dl) {
    f <- list.files(DIR_CACHE, pattern = "\\.csv$", full.names = TRUE)
    if (length(f) > 0) { file.remove(f); cat(glue("  {length(f)} caches supprimes\n")) }
  }
  recuperer_safran()
  if (length(verifier_dalles(reparer = TRUE)) > 0) recuperer_safran()
  
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  rm(dt); gc()
  
  cat(">> Controle des trous internes...\n")
  verifier_trous(meteo$tmax[[1]])
  cat("\n")
  
  meteo$daylength <- calculer_daylength(RST_FRANCE, as.Date(time(meteo$tmean)))
  
  steno <- modeliser_steno(meteo, exposition = exposition, k_ecorce = k_ecorce)
  rm(meteo); gc()
  
  exporter_rasters(steno)
  exporter_geojson(steno, pas_essaimage = pas_essaimage, grille = grille)
  generer_cartes(steno)
  synthese_console(steno)
  
  cat(glue("  Duree totale : {round((proc.time()-t0)[['elapsed']]/60, 1)} min\n"))
  cat(glue("  Sorties      : {DIR_OUT}/\n\n"))
  
  invisible(steno)
}


# =============================================================================
# 14. LANCEMENT
# =============================================================================

# --- Etape 0 : verifier l'unite du rayonnement dans le cache existant ---
if (FALSE) {
  diagnostiquer_rayonnement()
}

# --- Test sur une seule dalle (Massif Central) ---
if (FALSE) {
  dalle_mc <- DALLES[DALLES$xmin == 600000 & DALLES$ymin == 6200000, ]
  if (nrow(dalle_mc) == 0) dalle_mc <- DALLES[8, ]
  
  telecharger_dalle(dalle_mc,
                    date_deb = as.Date(glue("{ANNEE}-01-01")),
                    date_fin = as.Date(glue("{ANNEE}-07-31")))
  
  dt_t    <- lire_cache_complet(dalles = dalle_mc)
  meteo_t <- assembler_rasters_nationaux(dt_t)
  meteo_t$daylength <- calculer_daylength(RST_FRANCE, as.Date(time(meteo_t$tmean)))
  
  steno_t <- modeliser_steno(meteo_t)
  synthese_console(steno_t)
}

# --- Pipeline complet France ---
# DESACTIVE : la version production passe par R/pipeline_steno.R (charge le
# jeu SAFRAN partage, appelle modeliser_steno() directement).
if (FALSE) {
  # Scenario central : ecorce epaisse de bas de tronc, amortissement 0.2
  steno <- pipeline_steno(exposition = "sunny", k_ecorce = 0.2)

  # Encadrement de l'incertitude sur la temperature sous ecorce :
  #   borne haute = equation de Baier brute (epicea, ecorce fine)
  #   borne basse = temperature de l'air (aucun echauffement radiatif)
  # steno_haut <- pipeline_steno(exposition = "sunny", k_ecorce = 1.0)
  # steno_bas  <- pipeline_steno(exposition = "shaded")
}

# --- Calibration du seuil d'essaimage DJ_VOL ---
if (FALSE) {
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  meteo$daylength <- calculer_daylength(RST_FRANCE, as.Date(time(meteo$tmean)))
  calibrer_seuil_vol(meteo)
}

# --- Reparation des caches tronques ---
if (FALSE) {
  verifier_dalles(reparer = TRUE)
  recuperer_safran()
}

# --- NOTE : calibration terrain de K_ECORCE ---
# Protocole minimal, une saison :
#   - 5 a 10 chablis / grumes de pin en aire de stockage, faces SE et NO
#   - sondes thermocouple inserees entre phloeme et aubier a 1-1,5 m
#     (c'est la zone reellement colonisee par le stenographe)
#   - enregistrement horaire, avril a septembre
#   - regression de (BT_max_observee - Tair_max_SAFRAN) sur SSI_Q
#     -> la pente donne directement K_ECORCE * 0.0008518
# Comparer imperativement grume au sol (forte insolation, ecorce chauffee) et
# arbre sur pied (tronc ombrage par le houppier) : ce sont deux regimes
# thermiques distincts, et le stenographe exploite les deux.

# --- Test unitaire rapide du calcul de degres-jours ---
if (FALSE) {
  # Journee 15-25 C, seuil 10.9 : DJ attendu ~ 9.1 (moyenne 20 - 10.9)
  stopifnot(abs(dj_sinus(15, 25, 10.9, 36) - 9.1) < 0.3)
  # Journee entierement sous le seuil
  stopifnot(dj_sinus(2, 8, 10.9, 36) == 0)
  # Cutoff superieur actif : 30-42 C -> tronque a 36
  stopifnot(dj_sinus(30, 42, 10.9, 36) < (36 - 10.9) + 0.01)
  # Seuil traversant la sinusoide : positif mais < moyenne - seuil
  x <- dj_sinus(5, 20, 10.9, 36)
  stopifnot(x > 0 && x < (12.5 - 10.9) + 5)
  cat("Tests dj_sinus OK\n")
}