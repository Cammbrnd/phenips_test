# =============================================================================
# SAFRAN COMMON — fonctions partagees par PHENIPS, CHAPY, STENO, Monochamus
# -----------------------------------------------------------------------------
# Regroupe tout ce qui etait duplique a l'identique dans les 4 scripts :
# grille de dalles, telechargement, cache CSV, lecture, assemblage rasters,
# daylength. Chaque modele n'a plus qu'a sourcer ce fichier puis definir SES
# propres parametres biologiques et sa fonction modeliser_xxx().
#
# IMPORTANT : on telecharge systematiquement le SUPERSET des 4 variables
# (T_Q, TINF_H_Q, TSUP_H_Q, SSI_Q), meme pour les modeles qui n'utilisent pas
# le rayonnement (Monochamus). Cela permet un cache 100% partage entre les
# 4 pipelines : un seul jeu de requetes API GeoSAS pour tous les modeles.
#
# v2 (controle qualite de couverture) :
#   - dalles de 100 km (et non 200) avec recouvrement de 16 km ;
#   - bbox formatees sans notation scientifique ;
#   - cache renomme (safran_v2_*) pour ne jamais relire un ancien fichier
#     tronque ;
#   - detection des colonnes/lignes de mailles manquantes, re-telechargement
#     par bandes etroites, comblement limite (<= 2 mailles de large), puis
#     ECHEC BRUYANT si des trous persistent (voir SAFRAN_QC_STRICT).
# =============================================================================

library(terra)
library(httr2)
library(glue)
library(data.table)
library(lubridate)

# --- Periode (identique pour tous les modeles) ---
ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))
DATE_FIN_SAFRAN <- Sys.Date() - 5
if (DATE_FIN_SAFRAN < DATE_DEB) DATE_FIN_SAFRAN <- DATE_DEB

# --- API GeoSAS ---
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
GEOSAS_VARS <- "T_Q,TINF_H_Q,TSUP_H_Q,SSI_Q"   # superset partage par les 4 modeles

# --- Grille SAFRAN France metropolitaine (L93) ---
FRANCE_XMIN <- 100000
FRANCE_XMAX <- 1200000
FRANCE_YMIN <- 6050000
FRANCE_YMAX <- 7150000
SAFRAN_RES  <- 8000
DALLE_KM    <- 100                 # 200 avant la v2 : requetes plus petites
DALLE_M     <- DALLE_KM * 1000
MARGE_DALLE <- 2 * SAFRAN_RES      # recouvrement entre dalles voisines (16 km)
CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

# --- Controle qualite ---
# TRUE (defaut) : le job ECHOUE si des mailles manquent ou si des dalles
# restent en erreur. FALSE : simple avertissement (repli d'urgence, a poser
# explicitement dans le workflow : SAFRAN_QC_STRICT: "false").
SAFRAN_QC_STRICT <- !identical(tolower(Sys.getenv("SAFRAN_QC_STRICT", "true")), "false")

# --- Repertoire de cache PARTAGE (un seul, pour les 4 modeles) ---
DIR_CACHE <- "data/cache_safran"
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

#' Signale un probleme de qualite : erreur si strict, avertissement sinon.
#' Le prefixe ::error:: / ::warning:: cree une annotation dans GitHub Actions.
echec_qc <- function(msg) {
  if (SAFRAN_QC_STRICT) {
    message("::error::", msg)
    stop(msg, call. = FALSE)
  }
  message("::warning::", msg, " (SAFRAN_QC_STRICT=false : on continue)")
  invisible(FALSE)
}


# =============================================================================
# GRILLE DE DALLES
# =============================================================================

creer_dalles <- function() {
  # "- 1" : evite une dalle de largeur nulle quand l'emprise est un multiple
  # exact de la taille de dalle (cas de 100 km).
  xs <- seq(FRANCE_XMIN, FRANCE_XMAX - 1, by = DALLE_M)
  ys <- seq(FRANCE_YMIN, FRANCE_YMAX - 1, by = DALLE_M)
  dalles      <- expand.grid(xmin = xs, ymin = ys)
  dalles$xmax <- pmin(dalles$xmin + DALLE_M, FRANCE_XMAX)
  dalles$ymax <- pmin(dalles$ymin + DALLE_M, FRANCE_YMAX)
  dalles$id   <- seq_len(nrow(dalles))
  dalles
}

DALLES <- creer_dalles()

creer_raster_france <- function() {
  rast(
    xmin = FRANCE_XMIN, xmax = FRANCE_XMAX,
    ymin = FRANCE_YMIN, ymax = FRANCE_YMAX,
    resolution = SAFRAN_RES, crs = CRS_L93
  )
}

RST_FRANCE <- creer_raster_france()


# =============================================================================
# TELECHARGEMENT + CACHE
# =============================================================================

chemin_cache <- function(dalle_id) {
  file.path(DIR_CACHE, glue("safran_v2_dalle{sprintf('%03d', dalle_id)}_{ANNEE}.csv"))
}

#' Fichiers de reparation (bandes etroites) : regeneres a chaque run, donc
#' supprimes au debut de recuperer_safran() pour ne jamais relire un ancien.
fichiers_patch <- function() {
  list.files(DIR_CACHE, pattern = glue("^safran_v2_patch_.*_{ANNEE}\\.csv$"),
             full.names = TRUE)
}

csv_cache_valide <- function(f) {
  premieres_lignes <- tryCatch(readLines(f, n = 2, warn = FALSE, encoding = "UTF-8"),
                                error = function(e) "")
  if (length(premieres_lignes) == 0) return(FALSE)
  if (grepl("^\\s*\\{", premieres_lignes[1])) return(FALSE)
  TRUE
}

cache_est_valide <- function(dalle_id) {
  f <- chemin_cache(dalle_id)
  if (!file.exists(f)) return(FALSE)
  tryCatch({
    entete <- names(fread(f, nrows = 0))
    if (length(entete) == 0) return(TRUE)
    col_t <- entete[grep("^time$|^date$|datetime", tolower(entete))][1]
    if (is.na(col_t)) return(TRUE)
    df <- fread(f, select = col_t)
    if (nrow(df) == 0) return(TRUE)
    max_date <- as.Date(max(df[[col_t]], na.rm = TRUE))
    if (is.na(max_date) || is.infinite(as.numeric(max_date))) return(TRUE)
    max_date >= DATE_FIN_SAFRAN
  }, error = function(e) FALSE)
}

#' Telecharge une bbox L93 vers un fichier CSV.
#' Retourne list(statut, code) avec statut parmi :
#'   "ok"     CSV valide ecrit dans `fichier`
#'   "vide"   pas de donnees (HTTP 4xx, corps vide ou non CSV) : normal en mer
#'   "erreur" panne reseau / timeout / HTTP 5xx / 429 apres les relances
#' Un fichier ancien est supprime des que le statut n'est pas "ok" : on ne
#' relit jamais silencieusement une version perimee.
telecharger_bbox <- function(xmin, ymin, xmax, ymax, fichier,
                             vars = GEOSAS_VARS,
                             date_deb = DATE_DEB, date_fin = DATE_FIN_SAFRAN) {
  # sprintf("%.0f") : jamais de notation scientifique ("5e+05") dans l'URL
  bbox    <- sprintf("%.0f,%.0f,%.0f,%.0f", xmin, ymin, xmax, ymax)
  periode <- glue("{date_deb}/{date_fin}")

  req <- request(GEOSAS_URL) |>
    req_url_query(bbox = bbox, crs = GEOSAS_CRS, `parameter-name` = vars,
                  f = "CSV", datetime = periode) |>
    req_timeout(180) |>
    req_error(is_error = \(resp) FALSE) |>      # on lit nous-memes le code HTTP
    req_retry(max_tries = 5, backoff = \(i) 30 * i,
              is_transient = \(resp) resp_status(resp) %in% c(429, 503, 504))

  resultat <- function(statut, code = NA_integer_) {
    if (statut != "ok" && file.exists(fichier)) file.remove(fichier)
    list(statut = statut, code = code)
  }

  resp <- tryCatch(req_perform(req), error = function(e) {
    message("  ERREUR reseau/timeout (", bbox, ") : ", conditionMessage(e))
    NULL
  })
  if (is.null(resp)) return(resultat("erreur"))

  code <- resp_status(resp)
  if (code >= 500 || code == 429) return(resultat("erreur", code))
  if (code != 200)                return(resultat("vide", code))

  corps <- resp_body_raw(resp)
  if (length(corps) < 50) return(resultat("vide", code))
  writeBin(corps, fichier)
  if (!csv_cache_valide(fichier)) {
    message("  reponse non CSV (", bbox, ") : ",
            substr(paste(readLines(fichier, n = 1, warn = FALSE), collapse = ""), 1, 120))
    return(resultat("vide", code))
  }
  list(statut = "ok", code = code)
}

telecharger_dalle <- function(dalle) {
  telecharger_bbox(dalle$xmin - MARGE_DALLE, dalle$ymin - MARGE_DALLE,
                   dalle$xmax + MARGE_DALLE, dalle$ymax + MARGE_DALLE,
                   chemin_cache(dalle$id))
}

recuperer_safran <- function(dalles = DALLES) {
  cat(glue(">> Telechargement SAFRAN partage ({DATE_DEB} -> {DATE_FIN_SAFRAN})...\n"))
  cat(glue("   {nrow(dalles)} dalles de {DALLE_KM} km, recouvrement {MARGE_DALLE/1000} km\n"))
  unlink(fichiers_patch())     # les reparations sont refaites a chaque run
  n <- nrow(dalles)
  n_cache <- 0L; n_dl <- 0L; n_vide <- 0L
  ids_erreur <- integer(0)
  for (i in seq_len(n)) {
    cat(glue("\r  Dalle {i}/{n} | cache={n_cache} dl={n_dl} vides={n_vide} err={length(ids_erreur)}  "))
    d <- dalles[i, ]
    if (cache_est_valide(d$id)) { n_cache <- n_cache + 1L; next }
    res <- telecharger_dalle(d)
    if (res$statut == "ok") n_dl <- n_dl + 1L
    else if (res$statut == "vide") n_vide <- n_vide + 1L
    else ids_erreur <- c(ids_erreur, d$id)
    Sys.sleep(1.0)
  }
  cat("\n")

  # Une seconde chance pour les dalles en erreur (panne passagere de l'API)
  if (length(ids_erreur) > 0) {
    cat(glue("  {length(ids_erreur)} dalle(s) en erreur -> nouvelle tentative dans 60 s\n"))
    Sys.sleep(60)
    restantes <- integer(0)
    for (id in ids_erreur) {
      res <- telecharger_dalle(dalles[dalles$id == id, ])
      if (res$statut == "ok") n_dl <- n_dl + 1L
      else if (res$statut == "vide") n_vide <- n_vide + 1L
      else restantes <- c(restantes, id)
      Sys.sleep(1.0)
    }
    ids_erreur <- restantes
  }

  cat(glue("  Bilan : {n_cache} cache | {n_dl} telecharges | {n_vide} vides (mer/hors domaine) | {length(ids_erreur)} en erreur\n\n"))
  if (length(ids_erreur) > 0) {
    echec_qc(glue("{length(ids_erreur)} dalle(s) SAFRAN en erreur apres relance : ",
                  "{paste(ids_erreur, collapse = ', ')}. Des mailles terrestres manquent probablement."))
  }
  invisible(list(n_cache = n_cache, n_dl = n_dl, n_vide = n_vide, ids_erreur = ids_erreur))
}


# =============================================================================
# LECTURE + ASSEMBLAGE (superset des 4 variables, chaque modele pioche ce qu'il utilise)
# =============================================================================

lire_cache_complet <- function(dalles = DALLES) {
  cat(">> Lecture du cache SAFRAN partage...\n")
  fichiers <- c(chemin_cache(dalles$id), fichiers_patch())
  fichiers <- fichiers[file.exists(fichiers)]
  if (length(fichiers) == 0) stop("Aucun fichier cache. Lancer recuperer_safran() d'abord.")

  valides   <- vapply(fichiers, csv_cache_valide, logical(1))
  corrompus <- fichiers[!valides]
  fichiers  <- fichiers[valides]
  if (length(corrompus) > 0) {
    cat(glue("  !! {length(corrompus)} fichier(s) corrompus -> suppression\n"))
    file.remove(corrompus)
    if (length(fichiers) == 0) stop("Tous les fichiers cache etaient corrompus.")
  }

  temoin <- NULL
  for (f in fichiers) {
    dt_test <- tryCatch(fread(f, nrows = 5), error = function(e) NULL)
    if (!is.null(dt_test) && nrow(dt_test) > 0) { temoin <- dt_test; break }
  }
  if (is.null(temoin)) stop("Aucun fichier cache non-vide trouve.")

  detecter_col <- function(dt, patterns) {
    noms <- tolower(names(dt))
    for (p in patterns) {
      idx <- grep(p, noms, ignore.case = TRUE)
      if (length(idx) > 0) return(names(dt)[idx[1]])
    }
    NULL
  }

  col_date <- detecter_col(temoin, c("^time$", "^date$", "datetime"))
  col_x    <- detecter_col(temoin, c("^x$", "coord_x", "x_l93"))
  col_y    <- detecter_col(temoin, c("^y$", "coord_y", "y_l93"))
  col_tmoy <- detecter_col(temoin, c("t_q", "tmoy", "tmean"))
  col_tmin <- detecter_col(temoin, c("tinf_h_q", "tmin_q", "tmin"))
  col_tmax <- detecter_col(temoin, c("tsup_h_q", "tmax"))
  col_rad  <- detecter_col(temoin, c("ssi_q", "rad", "rayonnement"))

  liste <- lapply(seq_along(fichiers), function(i) {
    if (i %% 5 == 0) cat(glue("\r  Lecture {i}/{length(fichiers)}  "))
    tryCatch({
      dt <- fread(fichiers[i])
      if (nrow(dt) == 0) return(NULL)
      requis <- c(col_date, col_x, col_y, col_tmoy, col_tmin, col_tmax, col_rad)
      if (!all(requis %in% names(dt))) return(NULL)
      setnames(dt, requis, c("date","x","y","t_q","tmin_q","tsup_h_q","ssi_q"),
               skip_absent = FALSE)
      dt[, .(date, x, y, t_q, tmin_q, tsup_h_q, ssi_q)]
    }, error = function(e) NULL)
  })
  cat("\n")
  liste <- Filter(Negate(is.null), liste)
  dt <- rbindlist(liste, use.names = TRUE, fill = TRUE)
  dt[, date := as.Date(date)]
  dt <- dt[date >= DATE_DEB & date <= DATE_FIN_SAFRAN]
  # Les dalles se recouvrent : un meme point peut apparaitre deux fois.
  # unique() garde la premiere occurrence (les dalles precedent les patchs).
  dt <- unique(dt, by = c("date", "x", "y"))

  # Conversion rayonnement J/cm2 -> Wh/m2 (voir diagnostiquer_rayonnement() dans
  # steno_france_safran_v1.R si tu veux revalider cette hypothese d'unite)
  dt[, ssi_q := ssi_q * (10000 / 3600)]

  cat(glue("  {nrow(dt)} lignes | {uniqueN(dt$date)} jours | ",
           "{uniqueN(interaction(dt$x, dt$y))} tuiles SAFRAN\n\n"))
  dt
}


# =============================================================================
# CONTROLE QUALITE DE COUVERTURE
# =============================================================================
# Une maille est "manquante" si elle est vide alors que des mailles valides
# l'encadrent a gauche ET a droite (ou en haut ET en bas) dans un rayon de K
# mailles. Une ligne/colonne n'est signalee que si elle compte au moins
# `seuil` mailles manquantes : les petits golfes et estuaires (quelques
# mailles isolees) ne declenchent donc pas d'alerte.
#
# Fonction en R de base, sans dependance : testable hors pipeline.

#' @param ok matrice logique (lignes = nord -> sud) : TRUE si la maille a des donnees
trouver_trous_matrice <- function(ok, K = 3L, seuil = 8L) {
  nr <- nrow(ok); nc <- ncol(ok)
  # m_dec[i, j] = ok[i - dl, j - dc]  (hors grille : FALSE)
  decaler <- function(dl, dc) {
    out <- matrix(FALSE, nr, nc)
    ii <- seq_len(nr) - dl; jj <- seq_len(nc) - dc
    vi <- ii >= 1 & ii <= nr; vj <- jj >= 1 & jj <= nc
    out[which(vi), which(vj)] <- ok[ii[vi], jj[vj]]
    out
  }
  g <- d <- h <- b <- matrix(FALSE, nr, nc)
  for (k in seq_len(K)) {
    g <- g | decaler(0L,  k)
    d <- d | decaler(0L, -k)
    h <- h | decaler( k, 0L)
    b <- b | decaler(-k, 0L)
  }
  trou_h <- !ok & g & d     # colonne(s) manquante(s) entre deux zones valides
  trou_v <- !ok & h & b     # ligne(s) manquante(s)
  colonnes <- which(colSums(trou_h) >= seuil)
  lignes   <- which(rowSums(trou_v) >= seuil)
  cel <- matrix(FALSE, nr, nc)
  if (length(colonnes)) cel[, colonnes] <- trou_h[, colonnes, drop = FALSE]
  if (length(lignes))   cel[lignes, ]   <- cel[lignes, , drop = FALSE] | trou_v[lignes, , drop = FALSE]
  list(cellules = cel, colonnes = colonnes, lignes = lignes, n_cellules = sum(cel))
}

#' Matrice de presence (TRUE = au moins un point SAFRAN dans la maille).
matrice_presence <- function(dt) {
  pts   <- unique(dt[, .(x, y)])
  cells <- cellFromXY(RST_FRANCE, as.matrix(pts))
  cells <- cells[!is.na(cells)]
  ok <- rep(FALSE, ncell(RST_FRANCE))
  ok[cells] <- TRUE
  matrix(ok, nrow = nrow(RST_FRANCE), ncol = ncol(RST_FRANCE), byrow = TRUE)
}

#' Affiche les trous trouves, avec les coordonnees REELLES des points voisins
#' (permet de distinguer "le serveur n'a pas renvoye ces points" d'un simple
#' decalage de grille).
afficher_trous <- function(diag, dt = NULL) {
  for (j in diag$colonnes) {
    xc <- xFromCol(RST_FRANCE, j)
    li <- which(diag$cellules[, j])
    y  <- yFromRow(RST_FRANCE, li)
    cat(glue("  [TROU] colonne x~{round(xc)} : {length(li)} mailles, y de {round(min(y))} a {round(max(y))}\n"))
    if (!is.null(dt)) {
      proches <- sort(unique(round(dt$x[abs(dt$x - xc) <= 1.5 * SAFRAN_RES])))
      cat("         x reels voisins dans les donnees :",
          paste(head(proches, 12), collapse = ", "), "\n")
    }
  }
  for (i in diag$lignes) {
    yc <- yFromRow(RST_FRANCE, i)
    co <- which(diag$cellules[i, ])
    x  <- xFromCol(RST_FRANCE, co)
    cat(glue("  [TROU] ligne y~{round(yc)} : {length(co)} mailles, x de {round(min(x))} a {round(max(x))}\n"))
  }
}

#' Re-telecharge chaque ligne/colonne manquante par une bande etroite
#' (3 mailles de large, 12 mailles de long) : une requete de geometrie
#' differente de la dalle d'origine.
reparer_trous <- function(diag, taille_bloc = 12L) {
  demi <- 1.5 * SAFRAN_RES
  bandes <- list()
  for (j in diag$colonnes) {
    xc <- xFromCol(RST_FRANCE, j)
    y  <- yFromRow(RST_FRANCE, which(diag$cellules[, j]))
    for (blk in split(y, ceiling(seq_along(y) / taille_bloc))) {
      bandes[[length(bandes) + 1L]] <- c(xc - demi, min(blk) - demi, xc + demi, max(blk) + demi)
    }
  }
  for (i in diag$lignes) {
    yc <- yFromRow(RST_FRANCE, i)
    x  <- xFromCol(RST_FRANCE, which(diag$cellules[i, ]))
    for (blk in split(x, ceiling(seq_along(x) / taille_bloc))) {
      bandes[[length(bandes) + 1L]] <- c(min(blk) - demi, yc - demi, max(blk) + demi, yc + demi)
    }
  }
  cat(glue("  Reparation : {length(bandes)} bande(s) etroite(s) a re-telecharger\n"))
  n_ok <- 0L
  for (k in seq_along(bandes)) {
    b <- bandes[[k]]
    f <- file.path(DIR_CACHE, glue("safran_v2_patch_{sprintf('%03d', k)}_{ANNEE}.csv"))
    res <- telecharger_bbox(b[1], b[2], b[3], b[4], f)
    if (res$statut == "ok") n_ok <- n_ok + 1L
    Sys.sleep(1.0)
  }
  cat(glue("  Reparation : {n_ok}/{length(bandes)} bande(s) recuperee(s)\n"))
  invisible(n_ok)
}

#' Controle la couverture des points SAFRAN lus, tente une reparation, relit.
#' Ne fait JAMAIS echouer le job : la decision finale est prise apres
#' assemblage et comblement (verifier_couverture_finale).
controler_et_reparer_couverture <- function(dt) {
  cat(">> Controle de couverture SAFRAN...\n")
  diag <- trouver_trous_matrice(matrice_presence(dt), K = 3L)
  if (diag$n_cellules == 0L) {
    cat("  Aucune maille manquante detectee.\n\n")
    return(dt)
  }
  cat(glue("  {diag$n_cellules} maille(s) manquante(s) detectee(s) :\n"))
  afficher_trous(diag, dt)
  n_rep <- tryCatch(reparer_trous(diag), error = function(e) {
    message("  Reparation impossible : ", conditionMessage(e)); 0L
  })
  if (n_rep == 0L) { cat("\n"); return(dt) }
  dt2   <- lire_cache_complet()
  diag2 <- trouver_trous_matrice(matrice_presence(dt2), K = 3L)
  cat(glue("  Apres reparation : {diag2$n_cellules} maille(s) manquante(s)\n\n"))
  dt2
}

#' Comble par la moyenne des voisins les colonnes/lignes manquantes de largeur
#' <= largeur_max mailles (champ climatique lisse a 8 km). Retourne la liste
#' meteo corrigee et le nombre de mailles comblees. Les mailles comblees sont
#' comptees et tracees dans le journal : ce n'est pas une donnee mesuree.
combler_trous_raster <- function(meteo, largeur_max = 2L) {
  presence <- app(meteo$tmax, fun = "max", na.rm = TRUE)
  ok   <- !is.na(as.matrix(presence, wide = TRUE))
  diag <- trouver_trous_matrice(ok, K = largeur_max)
  if (diag$n_cellules == 0L) return(list(meteo = meteo, n = 0L))

  cat(glue(">> Comblement de {diag$n_cellules} maille(s) manquante(s) par moyenne des voisins...\n"))
  masque <- rast(presence)
  values(masque) <- ifelse(as.vector(t(diag$cellules)), 1, NA_real_)

  for (nom in c("tmean", "tmin", "tmax", "rad")) {
    r  <- meteo[[nom]]
    f  <- focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    f  <- mask(f, masque)
    r2 <- cover(r, f)
    names(r2) <- names(r)
    time(r2)  <- time(r)
    meteo[[nom]] <- r2
  }
  message("::warning::", diag$n_cellules, " maille(s) SAFRAN absente(s) de la reponse de l'API ",
          "comblee(s) par interpolation (moyenne des voisins).")
  list(meteo = meteo, n = diag$n_cellules)
}

#' Dernier controle sur les rasters assembles : echec si des trous subsistent.
verifier_couverture_finale <- function(meteo) {
  presence <- app(meteo$tmax, fun = "max", na.rm = TRUE)
  diag <- trouver_trous_matrice(!is.na(as.matrix(presence, wide = TRUE)), K = 3L)
  if (diag$n_cellules == 0L) {
    cat("  Couverture finale : OK, aucune maille manquante.\n\n")
    return(invisible(TRUE))
  }
  afficher_trous(diag)
  echec_qc(glue("{diag$n_cellules} maille(s) SAFRAN manquante(s) apres reparation et comblement : ",
                "les cartes auraient des barres blanches."))
}

assembler_rasters_nationaux <- function(dt) {
  cat(">> Assemblage rasters nationaux (partage)...\n")
  dates <- sort(unique(dt$date))
  n_j   <- length(dates)
  vars_info <- list(
    list(col = "t_q",      nom = "tmean"),
    list(col = "tmin_q",   nom = "tmin"),
    list(col = "tsup_h_q", nom = "tmax"),
    list(col = "ssi_q",    nom = "rad")
  )
  rasters <- setNames(vector("list", 4), c("tmean", "tmin", "tmax", "rad"))
  for (vi in vars_info) {
    cat(glue("  Variable {vi$nom}...\n"))
    couches <- lapply(seq_along(dates), function(i) {
      d   <- dates[i]
      sub <- dt[date == d, .(x, y, v = get(vi$col))]
      sub <- sub[!is.na(v)]
      r <- RST_FRANCE
      if (nrow(sub) < 3) {
        values(r) <- NA_real_
      } else {
        vpts <- vect(as.data.frame(sub[, .(x, y, v)]), geom = c("x","y"), crs = CRS_L93)
        r <- rasterize(vpts, RST_FRANCE, field = "v", fun = mean, na.rm = TRUE)
      }
      names(r) <- as.character(d)
      r
    })
    rst <- rast(couches)
    time(rst) <- dates
    rasters[[vi$nom]] <- rst
  }
  cat(glue("  Rasters assembles : {n_j} jours x {ncell(RST_FRANCE)} cellules\n\n"))
  rasters
}

calculer_daylength <- function(rst_ref, dates) {
  cat(">> Calcul daylength (partage)...\n")
  lats <- crds(rst_ref) |> as.data.frame() |>
    vect(geom = c("x","y"), crs = CRS_L93) |> project(CRS_WGS) |>
    crds() |> (\(m) m[, 2])()
  duree <- function(lat_deg, doy) {
    decl <- -asin(0.39779 * cos(pi/180 * (0.98563*(doy+10) + 1.914*sin(pi/180 * 0.98563*(doy-2)))))
    cos_ha <- pmin(pmax(-tan(lat_deg*pi/180) * tan(decl), -1), 1)
    2 * acos(cos_ha) * 12 / pi
  }
  couches <- lapply(seq_along(dates), function(i) {
    r <- rst_ref
    values(r) <- duree(lats, yday(dates[i]))
    names(r) <- as.character(dates[i])
    r
  })
  rst <- rast(couches)
  time(rst) <- dates
  rst
}
