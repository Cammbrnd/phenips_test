# =============================================================================
# VALIDATION de l'extraction historique par maille (annee EN COURS uniquement)
# -----------------------------------------------------------------------------
# Compare les tables data/historique_tmp/<modele>_<annee>.csv aux GeoJSON publies par
# les pipelines de production (memes donnees SAFRAN, memes masques) : chaque centre de
# maille doit tomber dans un polygone de meme classe. Resultats en annotations ::notice::.
# Ignore les annees passees (pas de GeoJSON de reference).
# =============================================================================

suppressPackageStartupMessages(library(sf))
annee <- suppressWarnings(as.integer(Sys.getenv("ANNEE_SAFRAN", "")))
if (is.na(annee) || annee != as.integer(format(Sys.Date(), "%Y"))) {
  cat("Validation ignoree : pas une annee en cours.\n"); quit(status = 0)
}

pts <- function(d) st_transform(st_as_sf(d[, c("x", "y")], coords = c("x", "y"), crs = 2154), 4326)
lire <- function(f) tryCatch(utils::read.csv(file.path("data/historique_tmp", sprintf("%s_%d.csv", f, annee))),
                             error = function(e) NULL)
texte_vers_num <- function(x) {
  x <- trimws(as.character(x)); x[x == "-"] <- "0"
  suppl <- ifelse(grepl("s$", x), 0.5, 0)
  suppressWarnings(as.numeric(sub("s$", "", x))) + suppl
}
comparer <- function(nom, d, col, geojson, champ, transf = identity, test = function(a, b) a == b) {
  if (is.null(d) || !col %in% names(d)) { message("::warning::VALIDATION ", nom, " : colonne absente"); return(invisible()) }
  if (!file.exists(geojson)) { message("::warning::VALIDATION ", nom, " : ", geojson, " absent"); return(invisible()) }
  g <- st_read(geojson, quiet = TRUE)
  j <- st_join(pts(d), g[champ])
  ref <- transf(j[[champ]][seq_len(nrow(d))])
  ok  <- !is.na(ref) & !is.na(d[[col]])
  eg  <- test(d[[col]][ok], ref[ok])
  message(sprintf("::notice::VALIDATION %s : %d / %d mailles concordantes (%.1f %%) ; %d mailles hors polygone ou sans valeur",
                  nom, sum(eg), sum(ok), 100 * mean(eg), sum(!ok)))
  if (any(!eg)) {
    t <- head(sort(table(paste(d[[col]][ok][!eg], "->", ref[ok][!eg]), decreasing = TRUE), 6))
    message("::notice::VALIDATION ", nom, " ecarts (table -> geojson) : ", paste(names(t), as.integer(t), sep = " x", collapse = " | "))
  }
}
doy_iso <- function(x) as.integer(format(as.Date(as.character(x)), "%j"))

ph <- lire("phenips")
for (v in c("sunny", "shaded")) {
  comparer(paste("PHENIPS", v, "generations"), ph, paste0("gen_", v),
           sprintf("data/phenips_generations_%s_latest.geojson", v), "generation", texte_vers_num)
  comparer(paste("PHENIPS", v, "premier envol"), ph, paste0("doy_", v),
           sprintf("data/phenips_onset_%s_latest.geojson", v), "date_iso", doy_iso)
}
ch <- lire("chapy")
comparer("CHAPY generations", ch, "gen", "data/chapy_generations_latest.geojson", "generation", as.numeric)
comparer("CHAPY premier envol", ch, "doy_envol", "data/chapy_onset_latest.geojson", "date_iso", doy_iso)
st <- lire("steno")
comparer("STENO moyen generations", st, "gen_moyen", "data/steno_generations_latest.geojson", "nb_generations", as.numeric)
# essaimage STENO : polygones par tranche de 15 jours (doy_debut .. doy_fin)
if (!is.null(st) && file.exists("data/steno_essaimage_latest.geojson")) {
  g <- st_read("data/steno_essaimage_latest.geojson", quiet = TRUE)
  j <- st_join(pts(st), g[c("doy_debut", "doy_fin")])
  doy <- st$doy_moyen
  ok  <- !is.na(doy) & !is.na(j$doy_debut)
  dedans <- doy[ok] >= j$doy_debut[ok] & doy[ok] <= j$doy_fin[ok]
  message(sprintf("::notice::VALIDATION STENO moyen essaimage : %d / %d mailles dans la bonne tranche (%.1f %%)",
                  sum(dedans), sum(ok), 100 * mean(dedans)))
}
