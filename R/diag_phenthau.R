# =============================================================================
# DIAGNOSTIC PHENTHAUproc (processionnaire du chene) -- ISOLE, rien n'est integre
# -----------------------------------------------------------------------------
# Lance par le workflow en mode "diagnostic_phenthau". Ne touche ni aux 4 modeles
# ni aux GeoJSON publies : ecrit seulement dans data/diag/ (commite par le workflow).
#   A. test du package sur ses donnees d'exemple (load_test)
#   B. donnees SAFRAN journalieres sur 2 petites zones, 1 sept 2023 -> 30 sept 2025
#      (via les fonctions de la brique partagee R/safran_common.R, cache a part)
#      - dictionnaire GeoSAS des variables utilisees (T_Q, TINF_H_Q, TSUP_H_Q)
#      - phenthau(dailymeanminmax) pour 2024 et 2025 (+ serie tronquee en cours de saison,
#        + serie SANS le 1er sept N-1 pour mesurer ce que le package exige)
#      - controles de coherence : ordre des stades, NA, tmin <= tmean <= tmax
# Sorties : data/diag/phenthau_diag.txt, phenthau_diag_cellules.csv,
#           phenthau_geosas_dictionnaire.json
# =============================================================================

suppressMessages({
  library(terra); library(PHENTHAUproc); library(jsonlite); library(httr2)
})

DIR_DIAG <- "data/diag"
dir.create(DIR_DIAG, showWarnings = FALSE, recursive = TRUE)
F_RAPPORT <- file.path(DIR_DIAG, "phenthau_diag.txt")
if (file.exists(F_RAPPORT)) file.remove(F_RAPPORT)
rapport <- function(...) {
  l <- paste0(...)
  cat(l, "\n")
  cat(l, "\n", file = F_RAPPORT, append = TRUE)
}
tete <- function(titre) rapport("\n==== ", titre, " ====")

tete("Versions")
rapport("PHENTHAUproc ", as.character(packageVersion("PHENTHAUproc")),
        " | terra ", as.character(packageVersion("terra")),
        " | ", R.version.string)

# --- A. Donnees d'exemple du package -------------------------------------------
tete("A. load_test (donnees d'exemple FVA)")
fva <- PHENTHAUproc::load_test("SpatRaster")
rapport("names(fva) : ", paste(names(fva), collapse = ", "))
rapport("head(time(fva$tmean)) : ", paste(head(terra::time(fva$tmean)), collapse = " "))
res_fva <- suppressMessages(phenthau(fva, params = parameter("dailymeanminmax", year = 2020)))
rapport("names(res) : ", paste(names(res_fva), collapse = ", "))
rapport("nlyr(stages) : ", nlyr(res_fva$stages), " | premier/dernier jour : ",
        paste(range(as.character(terra::time(res_fva$stages))), collapse = " -> "))
rapport("codes stages : ", paste(sort(unique(as.vector(values(res_fva$stages)))), collapse = " "))

# --- B. SAFRAN ------------------------------------------------------------------
tete("B. SAFRAN journalier via la brique partagee")
source("R/safran_common.R")     # telecharger_bbox, lire_cache_complet, assembler_rasters_nationaux

# B1. Dictionnaire GeoSAS : on lit les metadonnees de la collection
url_coll <- sub("/cube$", "", GEOSAS_URL)
dico <- tryCatch({
  resp <- request(url_coll) |> req_url_query(f = "json") |> req_timeout(60) |>
    req_error(is_error = function(r) FALSE) |> req_perform()
  rapport("Collection ", url_coll, " : HTTP ", resp_status(resp))
  j <- fromJSON(resp_body_string(resp), simplifyVector = FALSE)
  pn <- j$parameter_names
  voulus <- c("T_Q", "TINF_H_Q", "TSUP_H_Q", "SSI_Q")
  out <- pn[intersect(voulus, names(pn))]
  rapport("Variables trouvees : ", paste(names(out), collapse = ", "),
          " (", length(pn), " variables dans la collection)")
  for (n in names(out)) {
    p <- out[[n]]
    rapport("  ", n, " : ", paste(unlist(p$description), collapse = " "),
            " | unite=", paste(unlist(p$unit), collapse = " "),
            " | propriete=", paste(unlist(p$observedProperty$label), collapse = " "))
  }
  writeLines(toJSON(list(collection = j$id, titre = j$title, extent = j$extent,
                         parameter_names = out), auto_unbox = TRUE, pretty = TRUE, null = "null"),
             file.path(DIR_DIAG, "phenthau_geosas_dictionnaire.json"))
  out
}, error = function(e) { rapport("Dictionnaire GeoSAS illisible : ", conditionMessage(e)); NULL })

# B2. Deux petites zones (rectangles de ~64 km, L93) : gradient plaine/montagne et zone chaude
ZONES <- list(
  alsace_vosges = c(xc = 1019929, yc = 6786449),
  toulouse      = c(xc = 574008,  yc = 6279127)
)
DEMI <- 32000
SERIE_DEB <- as.Date("2023-09-01")
SERIE_FIN <- as.Date("2025-09-30")
# 3 requetes par zone (taille raisonnable) : sept-dec 2023, 2024, janv-sept 2025
PERIODES <- list(
  c("2023-09-01", "2023-12-31"),
  c("2024-01-01", "2024-12-31"),
  c("2025-01-01", "2025-09-30")
)

# On pilote lire_cache_complet() / assembler_rasters_nationaux() de la brique via ses globales
# (periode + dossier de cache), SANS toucher au cache partage data/cache_safran.
DATE_DEB        <- SERIE_DEB
DATE_FIN_SAFRAN <- SERIE_FIN
ANNEE           <- 2025L

# Premiere date (ou NA) ou une matrice cellules x jours satisfait `cond`
premiere_date <- function(m, dates, cond) {
  ok <- cond(m)
  ok[is.na(ok)] <- FALSE
  idx <- max.col(ok, ties.method = "first")
  as.character(ifelse(rowSums(ok) > 0, as.character(dates[idx]), NA))
}

lancer_phenthau <- function(meteo, annee, fin = NULL, debut = NULL) {
  deb <- if (is.null(debut)) as.Date(paste0(annee - 1, "-09-01")) else debut
  fin <- if (is.null(fin)) as.Date(paste0(annee, "-09-30")) else fin
  d   <- as.Date(terra::time(meteo$tmean))
  k   <- which(d >= deb & d <= fin)
  x <- lapply(meteo[c("tmean", "tmin", "tmax")], function(r) {
    s <- r[[k]]; terra::time(s) <- d[k]; s
  })
  suppressMessages(phenthau(x, params = parameter("dailymeanminmax", year = annee)))
}

lignes <- list()
for (zone in names(ZONES)) {
  tete(paste0("Zone ", zone))
  z <- ZONES[[zone]]
  DIR_CACHE <- file.path(tempdir(), paste0("safran_", zone))
  dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)
  for (k in seq_along(PERIODES)) {
    f <- file.path(DIR_CACHE, sprintf("safran_v2_patch_%03d_%d.csv", k, ANNEE))
    r <- telecharger_bbox(z[["xc"]] - DEMI, z[["yc"]] - DEMI, z[["xc"]] + DEMI, z[["yc"]] + DEMI,
                          f, date_deb = PERIODES[[k]][1], date_fin = PERIODES[[k]][2])
    rapport("  requete ", PERIODES[[k]][1], " -> ", PERIODES[[k]][2], " : ", r$statut,
            " (HTTP ", r$code, ")")
    if (r$statut == "ok") {
      rapport("  en-tete CSV : ", paste(names(data.table::fread(f, nrows = 1)), collapse = ","))
    }
  }
  dt <- tryCatch(lire_cache_complet(dalles = data.frame(id = 999L)), error = function(e) {
    rapport("  lecture impossible : ", conditionMessage(e)); NULL })
  if (is.null(dt)) next
  rapport("  lignes : ", nrow(dt), " | jours : ", length(unique(dt$date)),
          " | mailles : ", nrow(unique(dt[, .(x, y)])),
          " | ", min(dt$date), " -> ", max(dt$date))
  attendus <- as.integer(SERIE_FIN - SERIE_DEB) + 1L
  rapport("  jours attendus : ", attendus, " | manquants : ", attendus - length(unique(dt$date)))
  rapport("  T_Q     (tmean) : ", paste(round(range(dt$t_q, na.rm = TRUE), 1), collapse = " a "), " | NA ", sum(is.na(dt$t_q)))
  rapport("  TINF_H_Q (tmin) : ", paste(round(range(dt$tmin_q, na.rm = TRUE), 1), collapse = " a "), " | NA ", sum(is.na(dt$tmin_q)))
  rapport("  TSUP_H_Q (tmax) : ", paste(round(range(dt$tsup_h_q, na.rm = TRUE), 1), collapse = " a "), " | NA ", sum(is.na(dt$tsup_h_q)))
  incoh <- dt[, mean(tmin_q > t_q + 1e-6 | t_q > tsup_h_q + 1e-6, na.rm = TRUE)]
  rapport("  part des lignes ou tmin <= tmean <= tmax est violee : ", round(100 * incoh, 3), " %")

  meteo <- assembler_rasters_nationaux(dt)       # grille nationale, NA hors zone (comme apres masquage)
  n_cel <- sum(!is.na(values(meteo$tmean[[1]])))
  rapport("  rasters : ", nlyr(meteo$tmean), " jours, ", n_cel, " mailles avec donnees")

  cas <- list(
    list(nom = "2024", annee = 2024L, fin = NULL, debut = NULL),
    list(nom = "2025", annee = 2025L, fin = NULL, debut = NULL),
    list(nom = "2025_tronque_15juin", annee = 2025L, fin = as.Date("2025-06-15"), debut = NULL),
    list(nom = "2025_sans_sept_N-1", annee = 2025L, fin = NULL, debut = as.Date("2025-01-01"))
  )
  for (cs in cas) {
    tete(paste0(zone, " / cas ", cs$nom))
    tt <- system.time(
      res <- tryCatch(lancer_phenthau(meteo, cs$annee, cs$fin, cs$debut),
                      error = function(e) conditionMessage(e), warning = function(w) w)
    )
    if (is.character(res)) { rapport("  ERREUR : ", res); next }
    if (inherits(res, "warning")) {
      rapport("  AVERTISSEMENT : ", conditionMessage(res))
      res <- tryCatch(suppressWarnings(lancer_phenthau(meteo, cs$annee, cs$fin, cs$debut)),
                      error = function(e) conditionMessage(e))
      if (is.character(res)) { rapport("  ERREUR apres avertissement : ", res); next }
    }
    rapport("  temps : ", round(tt[3], 1), " s | sorties : ", paste(names(res), collapse = ", "))
    dts <- as.Date(terra::time(res$stages))
    rapport("  stages : ", length(dts), " jours, ", min(dts), " -> ", max(dts))

    avec <- which(!is.na(values(meteo$tmean[[1]])[, 1]))      # numeros de mailles avec donnees
    S <- values(res$stages)[avec, , drop = FALSE]
    rapport("  mailles avec donnees : ", length(avec),
            " | mailles avec NA dans stages : ", sum(rowSums(is.na(S)) > 0),
            " | hors donnees non-NA : ", sum(!is.na(values(res$stages)[-avec, 1])))
    # codes monotones dans le temps (un stade ne revient pas en arriere)
    recul <- sum(apply(S, 1, function(v) any(diff(v[!is.na(v)]) < 0)))
    rapport("  mailles ou un stade recule dans le temps : ", recul)

    prem <- sapply(1:8, function(s) premiere_date(S, dts, function(m) m >= s))
    prem <- matrix(prem, nrow = length(avec))
    colnames(prem) <- paste0("stade_", 1:8)
    ordre_ko <- sum(apply(prem, 1, function(v) {
      v <- v[!is.na(v)]; length(v) > 1 && any(diff(as.numeric(as.Date(v))) < 0) }))
    rapport("  mailles ou l'ordre des dates de stades est incoherent : ", ordre_ko)
    rapport("  mailles sans stade 8 (adulte) a la fin de la serie : ", sum(is.na(prem[, "stade_8"])))

    ecl <- sapply(c("custers", "wagenhoff", "meurisse"), function(m)
      premiere_date(values(res[[m]])[avec, , drop = FALSE], dts, function(x) x > 0))
    bud <- premiere_date(values(res$budswelling)[avec, , drop = FALSE], dts, function(x) x > 0)
    leaf <- premiere_date(values(res$leafunfolding)[avec, , drop = FALSE], dts, function(x) x > 0)
    ppa1 <- premiere_date(values(res$ppa_biocide)[avec, , drop = FALSE], dts, function(x) x == 1)
    ppa2 <- premiere_date(values(res$ppa_biocide)[avec, , drop = FALSE], dts, function(x) x == 2)
    mort <- values(res$mortality)[avec, 1]

    tab <- data.frame(
      zone = zone, cas = cs$nom, annee = cs$annee,
      x = xFromCell(meteo$tmean, avec), y = yFromCell(meteo$tmean, avec),
      debourrement_budswelling = bud, debourrement_leafunfolding = leaf,
      eclosion_custers = ecl[, "custers"], eclosion_wagenhoff = ecl[, "wagenhoff"],
      eclosion_meurisse = ecl[, "meurisse"],
      prem, ppa_debut = ppa1, ppa_fin = ppa2, mortalite_pct = mort,
      stringsAsFactors = FALSE
    )
    lignes[[length(lignes) + 1L]] <- tab

    med <- function(v) if (all(is.na(v))) "NA" else as.character(as.Date(median(as.numeric(as.Date(v)), na.rm = TRUE)))
    rapport("  dates MEDIANES sur la zone :")
    rapport("    debourrement (budswelling) ", med(bud), " | feuilles (leafunfolding) ", med(leaf))
    rapport("    eclosion custers ", med(ecl[, "custers"]), " | wagenhoff ", med(ecl[, "wagenhoff"]),
            " | meurisse ", med(ecl[, "meurisse"]))
    rapport("    stades 1..8 : ", paste(sprintf("%d=%s", 1:8, apply(prem, 2, med)), collapse = " | "))
    rapport("    mortalite L1 (%) : min ", min(mort, na.rm = TRUE), " | mediane ", median(mort, na.rm = TRUE),
            " | max ", max(mort, na.rm = TRUE))
    rapport("    nb de mailles par mortalite >= 50 % : ", sum(mort >= 50, na.rm = TRUE))
  }
  rm(meteo, dt); invisible(gc())
}

if (length(lignes) > 0) {
  tout <- do.call(rbind, lignes)
  utils::write.csv(tout, file.path(DIR_DIAG, "phenthau_diag_cellules.csv"), row.names = FALSE)
  rapport("\nCSV : ", nrow(tout), " lignes -> phenthau_diag_cellules.csv")
}
tete("Fin du diagnostic PHENTHAU")
