# HISTORIQUE STENO-DJ : une annee (ANNEE_SAFRAN), 3 scenarios (pire / moyen / moins_pire), mailles a pin.
source("R/safran_common.R")
source("R/steno_functions.R")
source("R/historique_common.R")

meteo <- charger_meteo_masquee("pin")
cel   <- cellules_valides(meteo)
res   <- data.frame(annee = ANNEE_H, x = cel$x, y = cel$y)

for (sc in names(SCENARIOS)) {
  steno <- modeliser_steno(meteo, exposition = "shaded", par = modifyList(STENO, SCENARIOS[[sc]]))
  gen <- terra::values(steno$n_gen)[cel$id, 1]
  doy <- terra::values(steno$onset)[cel$id, 1]
  res[[paste0("gen_", sc)]] <- gen
  res[[paste0("doy_", sc)]] <- doy
  resume_colonne(sprintf("STENO %s %d generations", sc, ANNEE_H), gen)
  message("::notice::STENO ", sc, " ", ANNEE_H, " : premier envol jour ",
          paste(range(doy, na.rm = TRUE), collapse = " a "), " | mailles sans envol : ", sum(is.na(doy)))
  rm(steno); invisible(gc())
}
ecrire_table_annee(res, "steno")
