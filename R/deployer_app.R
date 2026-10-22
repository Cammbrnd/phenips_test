# =============================================================================
# DEPLOIEMENT DE L'APP SHINY SUR shinyapps.io (lance par .github/workflows/deploy-app.yml)
# -----------------------------------------------------------------------------
# Identifiants lus dans les secrets GitHub (Settings > Secrets and variables > Actions) :
#   SHINYAPPS_NAME    nom du compte shinyapps.io (ex. dsf-fr)
#   SHINYAPPS_TOKEN   jeton  (shinyapps.io > Account > Tokens > Show)
#   SHINYAPPS_SECRET  secret (idem, bouton "Show secret")
# Les versions de packages envoyees sont celles du runner (CRAN a jour), donc plus
# de probleme de package local trop ancien.
# =============================================================================

nom    <- Sys.getenv("SHINYAPPS_NAME")
token  <- Sys.getenv("SHINYAPPS_TOKEN")
secret <- Sys.getenv("SHINYAPPS_SECRET")
app    <- Sys.getenv("SHINYAPPS_APP", "phenips")

if (!nzchar(nom) || !nzchar(token) || !nzchar(secret)) {
  message("::error::Secrets manquants : SHINYAPPS_NAME, SHINYAPPS_TOKEN et SHINYAPPS_SECRET doivent etre definis dans les secrets du depot.")
  quit(status = 1)
}

library(rsconnect)
setAccountInfo(name = nom, token = token, secret = secret)

cat("Packages utilises par l'app (versions du runner) :\n")
for (p in c("shiny", "leaflet", "sf", "glue", "httr2", "jsonlite", "digest"))
  cat(sprintf("  %-9s %s\n", p, as.character(packageVersion(p))))

deployApp(appDir = "app", appName = app, account = nom, server = "shinyapps.io",
          forceUpdate = TRUE, launch.browser = FALSE)
cat("\n>> App deployee : https://", nom, ".shinyapps.io/", app, "/\n", sep = "")
