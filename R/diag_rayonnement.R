# =============================================================================
# DIAGNOSTIC RAYONNEMENT (SSI_Q) : verifie l'hypothese d'unite de safran_common.R
# -----------------------------------------------------------------------------
# safran_common.R convertit SSI_Q en Wh/m2 avec le facteur 10000/3600, c'est-a-dire
# qu'il suppose des J/cm2. Ce script lit le cache SAFRAN BRUT (avant conversion) et
# affiche des statistiques, avec un verdict sur l'unite. Lance par le workflow en
# mode "diagnostic_rayonnement" : ne telecharge rien, ne modifie rien.
#
# Ordres de grandeur du cumul quotidien de rayonnement global en France :
#   hiver ~ 200-400 J/cm2 | ete ~ 2000-2600 J/cm2   (annuel median ~ 1000-1400)
#   en MJ/m2 : hiver ~ 2-4 | ete ~ 20-26            (annuel median ~ 10-14)
#   en Wh/m2 : hiver ~ 550-1100 | ete ~ 5500-7200
# =============================================================================

library(data.table)

dir_cache <- "data/cache_safran"
fichiers  <- list.files(dir_cache, pattern = "\\.csv$", full.names = TRUE)
cat(">> Fichiers de cache :", length(fichiers), "\n")
if (length(fichiers) == 0) stop("Cache SAFRAN vide : lancer d'abord le workflow complet.")

# Echantillon reparti sur la France (au plus 12 fichiers non vides)
fichiers <- fichiers[file.size(fichiers) > 1000]
pas      <- max(1L, floor(length(fichiers) / 12))
fichiers <- fichiers[seq(1L, length(fichiers), by = pas)]
cat(">> Echantillon :", length(fichiers), "fichiers\n")

lire <- function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  noms <- tolower(names(dt))
  c_rad  <- names(dt)[grep("ssi_q|rayonnement|^rad$", noms)[1]]
  c_date <- names(dt)[grep("^time$|^date$|datetime", noms)[1]]
  if (is.na(c_rad) || is.na(c_date)) return(NULL)
  data.table(date = as.Date(dt[[c_date]]), ssi = as.numeric(dt[[c_rad]]), col = c_rad)
}
dt <- rbindlist(Filter(Negate(is.null), lapply(fichiers, lire)))
if (nrow(dt) == 0) stop("Aucune colonne de rayonnement lisible dans le cache.")

cat(">> Colonne lue :", unique(dt$col), "| lignes :", nrow(dt),
    "| du", format(min(dt$date)), "au", format(max(dt$date)), "\n\n")
dt[, mois := as.integer(format(date, "%m"))]

q <- function(x) round(quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE), 1)
glob <- q(dt$ssi)
cat("SSI_Q BRUT (tel que renvoye par GeoSAS) : p5 / mediane / p95 =",
    paste(glob, collapse = " / "), "\n")
mens <- dt[, .(mediane = round(median(ssi, na.rm = TRUE), 1)), by = mois][order(mois)]
cat("Mediane brute par mois :", paste0(mens$mois, ":", mens$mediane, collapse = "  "), "\n\n")

med <- as.numeric(glob[2])
if (is.na(med) || med <= 0) {
  verdict <- "INDETERMINE (mediane nulle ou NA)"; unite <- NA; facteur <- NA
} else if (med < 60) {
  unite <- "MJ/m2";  facteur <- 1e6 / 3600
} else if (med < 600) {
  unite <- "W/m2 (moyenne journaliere)"; facteur <- 24
} else {
  unite <- "J/cm2";  facteur <- 10000 / 3600
}
facteur_code <- 10000 / 3600
if (!is.na(facteur)) {
  verdict <- if (abs(facteur - facteur_code) < 1e-6) {
    "OK : l'unite observee est compatible avec le facteur du code (J/cm2 -> Wh/m2)"
  } else {
    sprintf("ECART : unite observee %s, facteur correct %.3f, facteur du code %.3f (x%.0f)",
            unite, facteur, facteur_code, facteur / facteur_code)
  }
}
cat("Unite deduite de la magnitude :", unite, "\n")
cat("Rayonnement CONVERTI par le code (Wh/m2) : mediane =",
    round(med * facteur_code), "| ete (juin-aout) =",
    round(median(dt[mois %in% 6:8, ssi], na.rm = TRUE) * facteur_code),
    "| hiver (dec-fev) =",
    round(median(dt[mois %in% c(12, 1, 2), ssi], na.rm = TRUE) * facteur_code), "\n")
cat("   attendu : ete ~ 5500-7200, hiver ~ 550-1100\n\n")
cat("VERDICT :", verdict, "\n")
message("::notice::SSI_Q brut mediane=", med, " ; unite deduite=", unite, " ; ", verdict)
