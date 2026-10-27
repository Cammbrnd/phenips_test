# HISTORIQUE PHENIPS-Clim : une annee (ANNEE_SAFRAN), 2 scenarios d'exposition, mailles a epicea.
source("R/safran_common.R")
source("R/phenips_functions.R")
source("R/historique_common.R")

meteo <- charger_meteo_masquee("epicea")
cel   <- cellules_valides(meteo)
res   <- data.frame(annee = ANNEE_H, x = cel$x, y = cel$y)

for (v in names(VARIANTES_PHENIPS)) {
  param <- VARIANTES_PHENIPS[[v]]
  pheno <- modeliser_phenips(meteo, exposure = param$exposure, scenario = param$scenario)
  rst   <- get_gen_rst(pheno, max(pheno$dates))
  gen   <- if (is.null(rst)) rep(NA_real_, nrow(cel)) else
             normaliser_generation(etiquettes_raster(rst))[cel$id]
  res[[paste0("gen_", v)]] <- gen
  res[[paste0("doy_", v)]] <- doy_premier(pheno$onset, pheno$dates)[cel$id]
  resume_colonne(sprintf("PHENIPS %s %d generations", v, ANNEE_H), gen)
  message("::notice::PHENIPS ", v, " ", ANNEE_H, " : premier envol jour ",
          paste(range(res[[paste0("doy_", v)]], na.rm = TRUE), collapse = " a "),
          " | mailles sans envol : ", sum(is.na(res[[paste0("doy_", v)]])))
  rm(pheno); invisible(gc())
}
ecrire_table_annee(res, "phenips")
