# =============================================================================
# DIAGNOSTIC BD FORET (WFS IGN Geoplateforme) : couche, champs, codes d'essences
# -----------------------------------------------------------------------------
# Lance par le workflow en mode "diagnostic_bdforet". Ne modifie rien : ecrit des
# fichiers texte dans data/diag/, que le workflow commite pour qu'on puisse les lire.
#   1. liste les couches WFS "FORESTINVENTORY" (versions de la BD Foret)
#   2. pour chacune : champs (DescribeFeatureType) + 3 objets d'exemple
#   3. echantillon de 18 zones de 20 x 20 km (L93) : table (essence, tfv) avec effectifs
# =============================================================================

library(httr2)
library(jsonlite)

DIR_DIAG <- "data/diag"
dir.create(DIR_DIAG, showWarnings = FALSE, recursive = TRUE)
WFS <- "https://data.geopf.fr/wfs/ows"

journal <- character(0)
note <- function(...) { l <- paste0(...); cat(l, "\n"); journal <<- c(journal, l) }

wfs_get <- function(params, timeout = 120) {
  req <- request(WFS) |> req_url_query(!!!params) |> req_timeout(timeout) |>
    req_error(is_error = function(resp) FALSE) |> req_retry(max_tries = 2)
  resp <- req_perform(req)
  list(code = resp_status(resp), txt = tryCatch(resp_body_string(resp), error = function(e) ""))
}

# --- 1. Couches ---------------------------------------------------------------
cap <- wfs_get(list(SERVICE = "WFS", VERSION = "2.0.0", REQUEST = "GetCapabilities"))
note("GetCapabilities : HTTP ", cap$code, " | ", nchar(cap$txt), " caracteres")
noms <- unique(regmatches(cap$txt, gregexpr("(?<=<Name>)[^<]+(?=</Name>)|(?<=:Name>)[^<]+(?=</)", cap$txt, perl = TRUE))[[1]])
foret <- grep("FOREST|FORET", noms, ignore.case = TRUE, value = TRUE)
note("Couches contenant FOREST/FORET (", length(foret), ") :")
for (n in foret) note("   ", n)
writeLines(noms, file.path(DIR_DIAG, "bdforet_toutes_couches.txt"))

cibles <- grep("FORESTINVENTORY", foret, value = TRUE)
if (length(cibles) == 0) cibles <- foret

# --- 2. Champs et exemples ------------------------------------------------------
champs_non_geom <- list()
for (tn in cibles) {
  note("\n=== ", tn, " ===")
  d <- wfs_get(list(SERVICE = "WFS", VERSION = "2.0.0", REQUEST = "DescribeFeatureType", TYPENAMES = tn))
  elts <- regmatches(d$txt, gregexpr("<xsd:element[^>]*name=\"[^\"]+\"[^>]*>", d$txt))[[1]]
  nom_c  <- sub(".*name=\"([^\"]+)\".*", "\\1", elts)
  type_c <- ifelse(grepl("type=\"", elts), sub(".*type=\"([^\"]+)\".*", "\\1", elts), "?")
  note("Champs : ", paste0(nom_c, " (", type_c, ")", collapse = ", "))
  champs_non_geom[[tn]] <- nom_c[!grepl("^gml:", type_c)]
  ex <- wfs_get(list(SERVICE = "WFS", VERSION = "2.0.0", REQUEST = "GetFeature", TYPENAMES = tn,
                     COUNT = 3, OUTPUTFORMAT = "application/json", SRSNAME = "EPSG:2154"))
  note("Exemple : HTTP ", ex$code)
  if (ex$code == 200) {
    p <- tryCatch(fromJSON(ex$txt)$features$properties, error = function(e) NULL)
    if (!is.null(p)) { note(paste(capture.output(print(p)), collapse = "\n")) }
  } else note(substr(ex$txt, 1, 300))
}

# --- 3. Echantillon d'essences ----------------------------------------------------------
zones <- list("Vosges" = c(988000, 6781000), "Jura" = c(922000, 6626000), "Alpes_N" = c(966000, 6472000),
  "Alpes_S" = c(979000, 6362000), "Auvergne" = c(692000, 6500000), "Morvan" = c(776000, 6667000),
  "Landes" = c(388000, 6352000), "Pyrenees" = c(455000, 6205000), "Cevennes" = c(748000, 6345000),
  "Ardennes" = c(815000, 6968000), "Bretagne" = c(255000, 6806000), "Limousin" = c(599000, 6512000),
  "Vosges_nord" = c(1015000, 6875000), "Corse" = c(1204000, 6142000), "Sologne" = c(625000, 6712000),
  "Provence" = c(933000, 6315000), "Normandie" = c(481000, 6882000), "Vercors" = c(897000, 6437000))

tn <- cibles[grep("V2", cibles)[1]]
if (is.na(tn)) tn <- cibles[1]
note("\n=== Echantillon sur : ", tn, " ===")
pn <- paste(setdiff(champs_non_geom[[tn]], c("gid", "id", "fid")), collapse = ",")
note("propertyName utilise : ", pn)
tout <- list()
t0 <- Sys.time()
for (z in names(zones)) {
  xy <- zones[[z]]
  bb <- sprintf("%d,%d,%d,%d,urn:ogc:def:crs:EPSG::2154", xy[1] - 10000, xy[2] - 10000, xy[1] + 10000, xy[2] + 10000)
  r <- wfs_get(list(SERVICE = "WFS", VERSION = "2.0.0", REQUEST = "GetFeature", TYPENAMES = tn,
                    COUNT = 4000, OUTPUTFORMAT = "application/json", SRSNAME = "EPSG:2154",
                    BBOX = bb, PROPERTYNAME = pn))
  if (r$code != 200) { note("  ", z, " : HTTP ", r$code, " ", substr(r$txt, 1, 200)); next }
  p <- tryCatch(fromJSON(r$txt)$features$properties, error = function(e) NULL)
  if (is.null(p) || nrow(p) == 0) { note("  ", z, " : vide"); next }
  p$zone <- z
  tout[[z]] <- p
  note("  ", z, " : ", nrow(p), " objets")
}
note("Duree echantillon : ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")

if (length(tout) > 0) {
  d <- do.call(rbind, lapply(tout, function(x) x[, intersect(names(x), names(tout[[1]])), drop = FALSE]))
  write.csv(d, file.path(DIR_DIAG, "bdforet_echantillon.csv"), row.names = FALSE)
  cols_cat <- setdiff(names(d), c("zone", "geometry"))
  for (cc in cols_cat) {
    if (is.character(d[[cc]]) && length(unique(d[[cc]])) < 200) {
      tb <- sort(table(d[[cc]]), decreasing = TRUE)
      note("\n--- valeurs de '", cc, "' (", length(tb), ") ---")
      note(paste(sprintf("%6d  %s", as.integer(tb), names(tb)), collapse = "\n"))
    }
  }
}

writeLines(journal, file.path(DIR_DIAG, "bdforet_diagnostic.txt"))
cat("\n>> Diagnostic ecrit dans", DIR_DIAG, "\n")
