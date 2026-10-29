# HISTORIQUE CHAPY : une annee (ANNEE_SAFRAN), mailles a epicea.
source("R/safran_common.R")
source("R/chapy_functions.R")
source("R/historique_common.R")

meteo <- charger_meteo_masquee("epicea")
cel   <- cellules_valides(meteo)
res   <- data.frame(annee = ANNEE_H, x = cel$x, y = cel$y)

pheno <- modeliser_chapy(meteo, mode = "max")
rst   <- get_gen_rst_chapy(pheno, date_fin_chapy(pheno))
gen   <- if (is.null(rst)) rep(NA_real_, nrow(cel)) else
           normaliser_generation(etiquettes_raster(rst))[cel$id]
res$gen      <- gen
res$doy_envol <- terra::values(raster_premier_envol_doy(pheno)$rst)[cel$id, 1]
resume_colonne(sprintf("CHAPY %d generations", ANNEE_H), gen)
message("::notice::CHAPY ", ANNEE_H, " : premier envol jour ",
        paste(range(res$doy_envol, na.rm = TRUE), collapse = " a "),
        " | mailles sans envol : ", sum(is.na(res$doy_envol)))
ecrire_table_annee(res, "chapy")
