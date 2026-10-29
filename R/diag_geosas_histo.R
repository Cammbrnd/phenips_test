# =============================================================================
# DIAGNOSTIC GEOSAS pour l'historique SAFRAN (1960 -> aujourd'hui)
# -----------------------------------------------------------------------------
# Repond a trois questions avant de lancer 66 annees de calcul :
#   1. de quelle date a quelle date la collection safran-isba est-elle disponible ?
#   2. combien de temps et de Mo pour UNE dalle sur UNE annee (1958, 1960, 1990, 2020...) ?
#   3. une requete sur plusieurs annees d'un coup est-elle acceptee ?
# Ecrit data/diag/geosas_histo.txt (commite par le workflow).
# =============================================================================

suppressPackageStartupMessages({ library(httr2); library(jsonlite) })
source("R/safran_common.R")

DIR_DIAG <- "data/diag"
dir.create(DIR_DIAG, showWarnings = FALSE, recursive = TRUE)
journal <- character(0)
note <- function(...) { l <- paste0(...); cat(l, "\n"); journal <<- c(journal, l) }

# --- 1. Metadonnees de la collection ----------------------------------------------------
base <- sub("/cube$", "", GEOSAS_URL)
r <- request(base) |> req_url_query(f = "json") |> req_timeout(60) |>
  req_error(is_error = function(x) FALSE) |> req_perform()
note("Collection : ", base, " -> HTTP ", resp_status(r))
if (resp_status(r) == 200) {
  meta <- tryCatch(fromJSON(resp_body_string(r), simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(meta)) {
    note("Titre : ", meta$title)
    ext <- meta$extent
    note("Etendue temporelle : ", paste(capture.output(str(ext$temporal)), collapse = " | "))
    note("Etendue spatiale   : ", paste(capture.output(str(ext$spatial)), collapse = " | "))
    note("Parametres : ", paste(names(meta$parameter_names), collapse = ", "))
    note("Cles du document : ", paste(names(meta), collapse = ", "))
  }
}

# --- 2. Une dalle (Vosges), plusieurs periodes -------------------------------------------
bbox <- "884000,6734000,1016000,6866000"          # dalle 100 km + marge de 16 km
tester <- function(periode) {
  f <- tempfile(fileext = ".csv"); on.exit(unlink(f), add = TRUE)
  t0 <- Sys.time()
  req <- request(GEOSAS_URL) |>
    req_url_query(bbox = bbox, crs = GEOSAS_CRS, `parameter-name` = GEOSAS_VARS,
                  f = "CSV", datetime = periode) |>
    req_timeout(900) |> req_error(is_error = function(x) FALSE)
  resp <- tryCatch(req_perform(req), error = function(e) { note("  ", periode, " : ERREUR ", conditionMessage(e)); NULL })
  if (is.null(resp)) return(invisible())
  duree <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  corps <- resp_body_raw(resp)
  note("  ", periode, " : HTTP ", resp_status(resp), " | ", duree, " s | ", round(length(corps) / 1e6, 2), " Mo")
  if (resp_status(resp) != 200 || length(corps) < 100) {
    note("     debut de la reponse : ", substr(rawToChar(corps[seq_len(min(300, length(corps)))]), 1, 300)); return(invisible())
  }
  writeBin(corps, f)
  d <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(d)) { note("     lecture CSV impossible"); return(invisible()) }
  col_t <- names(d)[grep("^time$|^date$|datetime", tolower(names(d)))][1]
  dates <- as.Date(substr(as.character(d[[col_t]]), 1, 10))
  note("     lignes : ", nrow(d), " | colonnes : ", paste(names(d), collapse = ","))
  note("     dates : ", min(dates, na.rm = TRUE), " -> ", max(dates, na.rm = TRUE),
       " (", length(unique(dates)), " jours)")
  for (v in intersect(c("T_Q", "TINF_H_Q", "TSUP_H_Q", "SSI_Q"), names(d)))
    note("     ", v, " : NA ", sum(is.na(d[[v]])), " | mediane ", round(median(d[[v]], na.rm = TRUE), 1),
         " | min ", round(min(d[[v]], na.rm = TRUE), 1), " | max ", round(max(d[[v]], na.rm = TRUE), 1))
  invisible()
}

note("\n=== Une annee, une dalle ===")
for (an in c(1958, 1959, 1960, 1975, 1990, 2005, 2020)) tester(paste0(an, "-01-01/", an, "-12-31"))
note("\n=== Plusieurs annees d'un coup ===")
tester("1960-01-01/1964-12-31")
tester("1960-01-01/1969-12-31")

writeLines(journal, file.path(DIR_DIAG, "geosas_histo.txt"))
cat("\n>> Diagnostic ecrit dans", DIR_DIAG, "\n")
