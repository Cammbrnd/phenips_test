# =============================================================================
# App Shiny PHENIPS-Clim — visualisation seule
# -----------------------------------------------------------------------------
# Ne fait AUCUN calcul : lit les deux GeoJSON publies par GitHub Actions
# et les affiche sur une carte leaflet. Prevu pour shinyapps.io (plan gratuit).
#
# A FAIRE avant deploiement : remplacer <TON_USER>/<TON_REPO> ci-dessous
# par le nom de ton depot GitHub (branche par defaut = main, a adapter si besoin).
# =============================================================================

library(shiny)
library(leaflet)
library(sf)
library(glue)

UTILISATEUR_GH <- "Cammbrnd"
REPO_GH        <- "phenips_test"
BRANCHE_GH     <- "main"

url_geojson <- function(nom_fichier) {
  glue("https://raw.githubusercontent.com/{UTILISATEUR_GH}/{REPO_GH}/{BRANCHE_GH}/data/{nom_fichier}")
}

URL_GEN   <- url_geojson("phenips_generations_latest.geojson")
URL_ONSET <- url_geojson("phenips_onset_latest.geojson")

ui <- fluidPage(
  titlePanel("PHENIPS-Clim — France metropolitaine"),
  div(
    style = "margin-bottom: 10px; color: #666; font-size: 13px;",
    textOutput("statut", inline = TRUE)
  ),
  leafletOutput("carte", height = "82vh")
)

server <- function(input, output, session) {

  lire_geojson_securise <- function(url) {
    tryCatch(sf::st_read(url, quiet = TRUE), error = function(e) NULL)
  }

  gen   <- reactive({ lire_geojson_securise(URL_GEN) })
  onset <- reactive({ lire_geojson_securise(URL_ONSET) })

  output$statut <- renderText({
    n_gen   <- if (is.null(gen()))   0 else nrow(gen())
    n_onset <- if (is.null(onset())) 0 else nrow(onset())
    glue("{n_gen} polygones generations | {n_onset} polygones onset")
  })

  output$carte <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = 2.5, lat = 46.5, zoom = 6)
  })

  observe({
    g <- gen()
    if (is.null(g) || nrow(g) == 0) return(invisible(NULL))
    leafletProxy("carte") |>
      clearGroup("generations") |>
      addPolygons(
        data        = g,
        fillColor   = ~couleur,
        color       = "#333333",
        weight      = 0.5,
        fillOpacity = 0.65,
        group       = "Generations",
        label       = ~label
      )
  })

  observe({
    o <- onset()
    if (is.null(o) || nrow(o) == 0) return(invisible(NULL))
    leafletProxy("carte") |>
      clearGroup("onset") |>
      addPolygons(
        data        = o,
        color       = "#333333",
        weight      = 0.5,
        fillOpacity = 0.5,
        group       = "Onset (1er essaimage)",
        label       = ~date_onset
      )
  })

  observe({
    leafletProxy("carte") |>
      addLayersControl(
        overlayGroups = c("Generations", "Onset (1er essaimage)"),
        options       = layersControlOptions(collapsed = FALSE)
      )
  })
}

shinyApp(ui, server)
