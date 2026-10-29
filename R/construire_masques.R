# =============================================================================
# CONSTRUCTION DES MASQUES D'ESSENCE SUR LA GRILLE SAFRAN (8 km, L93)
# -----------------------------------------------------------------------------
# Source : BD Foret V2 (IGN), couche WFS Geoplateforme
#          LANDCOVER.FORESTINVENTORY.V2:formation_vegetale, champ "essence".
#
#   masque_epicea.csv : mailles qui croisent au moins un peuplement "Sapin, epicea"
#                       -> utilise par PHENIPS (Ips typographus) et CHAPY (Pityogenes)
#   masque_pin.csv    : mailles qui croisent au moins un peuplement de pin
#                       (sylvestre, maritime, laricio/noir, d'Alep, a crochets/cembro,
#                       pins melanges, autre pin)  -> utilise par STENO (Ips sexdentatus)
#
#   masque_chene.csv  : mailles qui croisent au moins un peuplement de "Chenes decidus"
#                       (Quercus robur / petraea / pubescens ; les chenes sempervirents sont
#                       exclus) -> utilise par PHENTHAU (Thaumetopoea processionea)
#
# Colonnes : x, y (centre de maille, L93), part (part de la surface de la maille
# couverte par l'essence, 0-1). Le masque = toutes les lignes du fichier.
#
# ATTENTION : la BD Foret V2 attribue UNE essence par peuplement. "Sapin, epicea"
# ne designe que les peuplements purs (ou presque) de ces essences. L'epicea present
# dans un peuplement classe "Melange de coniferes" ou "Coniferes" n'est pas visible.
# Les masques sont donc prudents (ils peuvent manquer des mailles a epicea minoritaire).
#
# A lancer a la demande (workflow, entree "construire_masques") : la BD Foret change
# peu. Les CSV sont commites dans data/masques/ et lus par les pipelines.
# =============================================================================

suppressPackageStartupMessages({ library(terra); library(httr2); library(jsonlite); library(glue) })
source("R/safran_common.R")

WFS_URL     <- "https://data.geopf.fr/wfs/ows"
WFS_COUCHE  <- "LANDCOVER.FORESTINVENTORY.V2:formation_vegetale"
WFS_PAGE    <- 5000L
TUILE_MAILLES <- 5L                       # tuiles de 5 x 5 mailles = 40 km, alignees sur la grille
DIR_MASQUES <- "data/masques"
dir.create(DIR_MASQUES, showWarnings = FALSE, recursive = TRUE)

# Filtres CQL sur "essence". LIKE evite d'envoyer des accents dans l'URL.
#   "Sapin, epicea"                                    -> 'Sapin,%'
#   "Pin sylvestre", "Pin maritime", "Pin laricio, pin noir", "Pin d'Alep",
#   "Pin a crochets, pin cembro", "Pins melanges", "Pin autre" -> 'Pin%'
#   "Chenes decidus" (accents ^ et ') -> 'Ch%nes d%cidus' (exclut "Chenes sempervirents")
ESSENCES <- list(
  epicea = "essence LIKE 'Sapin,%'",
  pin    = "essence LIKE 'Pin%'",
  chene  = "essence LIKE 'Ch%nes d%cidus'"
)

# Variable d'environnement MASQUES (ex. "chene") : ne construit que ces masques.
# Vide = tous. Les autres CSV et leur entree dans masques_meta.json sont conserves.
NOMS_A_CONSTRUIRE <- trimws(strsplit(Sys.getenv("MASQUES", ""), ",")[[1]])
NOMS_A_CONSTRUIRE <- NOMS_A_CONSTRUIRE[nzchar(NOMS_A_CONSTRUIRE)]
if (length(NOMS_A_CONSTRUIRE) == 0) NOMS_A_CONSTRUIRE <- names(ESSENCES)
if (!all(NOMS_A_CONSTRUIRE %in% names(ESSENCES)))
  stop("MASQUES inconnu(s) : ", paste(setdiff(NOMS_A_CONSTRUIRE, names(ESSENCES)), collapse = ", "), call. = FALSE)

wfs_page <- function(cql, start, avec_propertyname = TRUE) {
  params <- list(SERVICE = "WFS", VERSION = "2.0.0", REQUEST = "GetFeature",
                 TYPENAMES = WFS_COUCHE, OUTPUTFORMAT = "application/json",
                 SRSNAME = "EPSG:2154", COUNT = WFS_PAGE, STARTINDEX = start,
                 SORTBY = "id", CQL_FILTER = cql)
  if (avec_propertyname) params$PROPERTYNAME <- "geom"
  req <- request(WFS_URL) |> req_url_query(!!!params) |> req_timeout(180) |>
    req_error(is_error = function(resp) FALSE) |>
    req_retry(max_tries = 4, backoff = function(i) 5 * i,
              is_transient = function(resp) resp_status(resp) %in% c(429, 500, 502, 503, 504))
  resp <- req_perform(req)
  list(code = resp_status(resp),
       txt  = tryCatch(resp_body_string(resp), error = function(e) ""))
}

compter_retournes <- function(txt) {
  m <- regmatches(txt, regexpr("\"numberReturned\"\\s*:\\s*[0-9]+", txt))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub(".*:\\s*", "", m))
}
compter_total <- function(txt) {
  m <- regmatches(txt, regexpr("\"numberMatched\"\\s*:\\s*[0-9]+", txt))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub(".*:\\s*", "", m))
}

lire_polygones <- function(txt) {
  f <- tempfile(fileext = ".geojson")
  on.exit(unlink(f), add = TRUE)
  writeLines(txt, f, useBytes = TRUE)
  v <- vect(f)
  # GeoServer renvoie des coordonnees L93 ; on force le CRS si GDAL a suppose du WGS84
  if (nrow(v) > 0 && abs(ext(v)[1]) > 360) crs(v) <- CRS_L93
  v
}

# Toutes les pages d'une tuile -> list(part, croise, xy) : un vecteur par maille de la tuile
traiter_tuile <- function(cql_essence, e, avec_pn) {
  r_tuile  <- crop(RST_FRANCE, e)
  part_v   <- numeric(ncell(r_tuile))
  croise_v <- numeric(ncell(r_tuile))
  cql <- sprintf("%s AND BBOX(geom,%.0f,%.0f,%.0f,%.0f,'EPSG:2154')",
                 cql_essence, e[1], e[3], e[2], e[4])
  start <- 0L; n_poly <- 0L
  repeat {
    p <- wfs_page(cql, start, avec_pn$val)
    if (p$code != 200 && avec_pn$val) {          # repli : sans PROPERTYNAME (decision gardee pour la suite)
      message("  [INFO] requete refusee avec PROPERTYNAME (HTTP ", p$code, "), nouvel essai sans")
      avec_pn$val <- FALSE
      p <- wfs_page(cql, start, FALSE)
    }
    if (p$code != 200) {
      msg <- regmatches(p$txt, regexpr("(?<=<ows:ExceptionText>)[^<]+", p$txt, perl = TRUE))
      stop(sprintf("WFS HTTP %s : %s", p$code, if (length(msg)) msg else substr(p$txt, 1, 300)), call. = FALSE)
    }
    n <- compter_retournes(p$txt)
    if (is.na(n) || n == 0L) break
    v <- lire_polygones(p$txt)
    if (nrow(v) > 0) {
      eb <- ext(v)
      if (eb[2] < e[1] - 50000 || eb[1] > e[2] + 50000 || eb[4] < e[3] - 50000 || eb[3] > e[4] + 50000)
        stop(sprintf("Polygones hors de la tuile (axes inverses ou CRS errone ?) : emprise recue %s",
                     paste(round(as.vector(eb)), collapse = " ")), call. = FALSE)
      pr <- values(rasterize(v, r_tuile, field = 1, touches = TRUE))[, 1]
      cv <- values(rasterize(v, r_tuile, cover = TRUE))[, 1]
      pr[is.na(pr)] <- 0; cv[is.na(cv)] <- 0
      croise_v <- pmax(croise_v, pr)
      part_v   <- part_v + cv
    }
    n_poly <- n_poly + n
    start  <- start + n
    total  <- compter_total(p$txt)
    if (!is.na(total) && start >= total) break
    if (start > 500000L) stop("Pagination anormale (plus de 500000 objets pour une tuile)", call. = FALSE)
  }
  idx <- which(croise_v == 1)
  list(n_poly = n_poly,
       d = if (length(idx) > 0) {
         xy <- xyFromCell(r_tuile, idx)
         data.frame(x = xy[, 1], y = xy[, 2], part = round(pmin(part_v[idx], 1), 4))
       } else NULL)
}

construire_masque <- function(nom) {
  cql_essence <- ESSENCES[[nom]]
  cat(glue("\n=== Masque '{nom}' : {cql_essence} ===\n"))
  xs <- seq(FRANCE_XMIN, FRANCE_XMAX - 1, by = TUILE_MAILLES * SAFRAN_RES)
  ys <- seq(FRANCE_YMIN, FRANCE_YMAX - 1, by = TUILE_MAILLES * SAFRAN_RES)
  tuiles <- expand.grid(x0 = xs, y0 = ys)
  cat(glue("  {nrow(tuiles)} tuiles de {TUILE_MAILLES * SAFRAN_RES / 1000} km\n"))

  avec_pn <- new.env(); avec_pn$val <- TRUE
  lignes <- list(); n_err <- 0L; n_poly_tot <- 0L; t0 <- Sys.time()
  for (i in seq_len(nrow(tuiles))) {
    x0 <- tuiles$x0[i]; y0 <- tuiles$y0[i]
    e <- ext(x0, min(x0 + TUILE_MAILLES * SAFRAN_RES, FRANCE_XMAX),
             y0, min(y0 + TUILE_MAILLES * SAFRAN_RES, FRANCE_YMAX))
    res <- NULL
    for (essai in 1:2) {                          # une 2e chance par tuile (erreur reseau ponctuelle)
      res <- tryCatch(traiter_tuile(cql_essence, e, avec_pn),
                      error = function(err) {
                        message("::warning::Tuile ", i, " (", x0, ",", y0, ") essai ", essai, " : ", conditionMessage(err))
                        NULL })
      if (!is.null(res)) break
      Sys.sleep(10)
    }
    if (is.null(res)) {
      n_err <- n_err + 1L
      if (i <= 3 && n_err == i) stop("Les 3 premieres tuiles echouent : requete WFS a revoir (voir ::warning:: ci-dessus)", call. = FALSE)
      next
    }
    n_poly_tot <- n_poly_tot + res$n_poly
    if (!is.null(res$d)) lignes[[length(lignes) + 1]] <- res$d
    if (i %% 50 == 0 || i == nrow(tuiles))
      cat(glue("  tuile {i}/{nrow(tuiles)} | {n_poly_tot} polygones | {sum(vapply(lignes, nrow, 1L))} mailles | ",
               "{round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1)} min\n"))
  }
  if (n_err > 0) stop(glue("{n_err} tuile(s) en erreur : masque '{nom}' incomplet, rien n'est ecrit"), call. = FALSE)

  d <- do.call(rbind, lignes)
  d <- d[order(d$y, d$x), ]
  fichier <- file.path(DIR_MASQUES, glue("masque_{nom}.csv"))
  write.csv(d, fichier, row.names = FALSE)
  cat(glue("  -> {fichier} : {nrow(d)} mailles (sur {ncell(RST_FRANCE)} cellules de la grille)\n"))
  cat(glue("     part de maille couverte : mediane {round(median(d$part), 3)} ; ",
           "mailles avec moins de 1% : {sum(d$part < 0.01)} ; plus de 10% : {sum(d$part >= 0.10)}\n"))
  list(nom = nom, n_mailles = nrow(d), n_polygones = n_poly_tot, filtre = cql_essence,
       date = format(Sys.Date()))
}

res <- lapply(NOMS_A_CONSTRUIRE, construire_masque)

# Construction partielle : on garde les entrees des masques non reconstruits
f_meta <- file.path(DIR_MASQUES, "masques_meta.json")
if (file.exists(f_meta) && !setequal(NOMS_A_CONSTRUIRE, names(ESSENCES))) {
  anciens <- tryCatch(read_json(f_meta)$masques, error = function(e) list())
  anciens <- Filter(function(m) !(m$nom %in% NOMS_A_CONSTRUIRE), anciens)
  res <- c(anciens, res)
}
meta <- list(
  source = "IGN - BD Foret V2 (LANDCOVER.FORESTINVENTORY.V2:formation_vegetale), champ essence",
  grille = list(res_m = SAFRAN_RES, crs = CRS_L93, xmin = FRANCE_XMIN, ymin = FRANCE_YMIN),
  date_construction = format(Sys.Date()),
  masques = res
)
write_json(meta, file.path(DIR_MASQUES, "masques_meta.json"), auto_unbox = TRUE, pretty = TRUE)

# Controles de coherence (mailles connues)
ctl <- list(
  list("epicea", "Vosges (988000, 6781000)",   988000, 6781000, TRUE),
  list("epicea", "Jura (922000, 6626000)",     922000, 6626000, TRUE),
  list("epicea", "Landes (388000, 6352000)",   388000, 6352000, FALSE),
  list("pin",    "Landes (388000, 6352000)",   388000, 6352000, TRUE),
  list("pin",    "Auvergne (692000, 6500000)", 692000, 6500000, NA),
  list("chene",  "Troncais, Allier (674771, 6618932)", 674771, 6618932, TRUE),
  list("chene",  "Massif du Mont-Blanc (999540, 6532918)", 999540, 6532918, FALSE)
)
ctl <- Filter(function(k) k[[1]] %in% NOMS_A_CONSTRUIRE, ctl)
for (k in ctl) {
  d <- read.csv(file.path(DIR_MASQUES, glue("masque_{k[[1]]}.csv")))
  cel <- cellFromXY(RST_FRANCE, cbind(k[[3]], k[[4]]))
  cel_d <- cellFromXY(RST_FRANCE, cbind(d$x, d$y))
  dedans <- cel %in% cel_d
  attendu <- k[[5]]
  statut <- if (is.na(attendu)) "info" else if (dedans == attendu) "OK" else "ECART"
  cat(glue("  controle {k[[1]]} / {k[[2]]} : {ifelse(dedans, 'dans le masque', 'hors masque')} [{statut}]\n"))
  if (statut == "ECART") message("::warning::Controle de coherence en ecart : ", k[[1]], " / ", k[[2]])
}
if (any(vapply(res, function(r) r$n_mailles < 300, logical(1))))
  stop("Masque anormalement petit (< 300 mailles) : verifier la requete WFS", call. = FALSE)
cat("\n>> Masques ecrits dans", DIR_MASQUES, "\n")
