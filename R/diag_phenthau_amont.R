# Bilan du fetch SAFRAN amont (mode workflow "test_amont") : ecrit data/diag/phenthau_amont_bilan.txt
suppressMessages(library(terra))
dir.create("data/diag", showWarnings = FALSE, recursive = TRUE)
m <- lapply(readRDS("data/cache_meteo/meteo_amont.rds"), terra::unwrap)
d <- as.Date(terra::time(m$tmean))
l <- c(
  paste0("variables : ", paste(names(m), collapse = ", ")),
  paste0("jours : ", length(d), " | ", min(d), " -> ", max(d)),
  paste0("jours manquants dans la serie : ", length(setdiff(as.character(seq(min(d), max(d), by = "day")), as.character(d)))),
  paste0("grille : ", nrow(m$tmean), " x ", ncol(m$tmean), " | mailles avec donnees (derniere couche) : ",
         sum(!is.na(values(m$tmean[[nlyr(m$tmean)]])))),
  paste0("tmean : ", paste(round(range(values(m$tmean), na.rm = TRUE), 1), collapse = " a "),
         " | tmin : ", paste(round(range(values(m$tmin), na.rm = TRUE), 1), collapse = " a "),
         " | tmax : ", paste(round(range(values(m$tmax), na.rm = TRUE), 1), collapse = " a ")),
  paste0("taille meteo_amont.rds : ", round(file.size("data/cache_meteo/meteo_amont.rds") / 1e6, 1), " Mo")
)
writeLines(l, "data/diag/phenthau_amont_bilan.txt")
cat(paste(l, collapse = "\n"), "\n")
