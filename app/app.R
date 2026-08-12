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
library(httr2)

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
    tryCatch({
      # On telecharge via httr2 (pile TLS distincte de celle de GDAL) puis on
      # lit le fichier en local. Evite l'erreur "schannel: CertGetCertificateChain
      # trust error" que GDAL peut declencher sur certains postes Windows
      # (proxy/antivirus qui intercepte le SSL) quand on lui passe l'URL directement.
      fichier_tmp <- tempfile(fileext = ".geojson")
      request(url) |>
        req_perform(path = fichier_tmp)
      sf::st_read(fichier_tmp, quiet = TRUE)
    }, error = function(e) {
      message("Erreur lecture GeoJSON (", url, "): ", conditionMessage(e))
      NULL
    })
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

    # Palette continue sur la date (date_iso -> Date), du plus precoce (clair)
    # au plus tardif (fonce). "YlOrRd" va du jaune (tot) au rouge fonce (tard).
    dates_num <- as.Date(o$date_iso)
    pal_onset <- colorNumeric(
      palette = "YlOrRd",
      domain  = dates_num,
      na.color = "#cccccc"
    )

    leafletProxy("carte") |>
      clearGroup("onset") |>
      removeControl("legende_onset") |>
      addPolygons(
        data        = o,
        fillColor   = ~pal_onset(as.Date(date_iso)),
        color       = "#333333",
        weight      = 0.5,
        fillOpacity = 0.75,
        group       = "Onset (1er essaimage)",
        label       = ~date_onset
      ) |>
      addLegend(
        position = "bottomright",
        pal      = pal_onset,
        values   = dates_num,
        title    = "1er essaimage",
        labFormat = labelFormat(date = TRUE),
        layerId  = "legende_onset",
        group    = "Onset (1er essaimage)"
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
